;+
; NAME:
;   xy2traceset
;
; PURPOSE:
;   Convert from an array of x,y positions to a trace set
;
; CALLING SEQUENCE:
;   xy2traceset, xpos, ypos, tset, [ invvar=, func=func, ncoeff=ncoeff, $
;    xmin=xmin, xmax=xmax, maxiter=maxiter, inputfunc=inputfunc, $
;    inmask=inmask, outmask=outmask, yfit=yfit, inputans=inputans, $
;    _EXTRA=EXTRA ]
;
; INPUTS:
;   xpos       - X positions corresponding to YPOS as an [nx,Ntrace] array
;   ypos       - Y centers as an [nx,ntrace] array
;
; OPTIONAL KEYWORDS:
;   invvar     - Inverse variance for weighted fits.
;   func       - Function for trace set; options are:
;                'poly'
;                'legendre'
;                'chebyshev'
;                'chebyshev_split'
;                Default to 'legendre'
;   ncoeff     - Number of coefficients in fit; default to 3
;   xmin       - Explicitly set XMIN for trace set rather than using minimum
;                in XPOS
;   xmax       - Explicitly set XMAX for trace set rather than using maximum
;                in XPOS
;   maxiter    - Maximum number of rejection iterations; set to 0 for no
;                rejection; default to 10.
;   inmask     - Mask set to 1 for good points and 0 for rejected points;
;                same dimensions as XPOS, YPOS.  Points rejected by INMASK
;                are always rejected from the fits (the rejection is "sticky"),
;                and will also be marked as rejected in OUTMASK.
;   inputans   - ???
;   inputfunc  - An array which matches the size of ypos, which is multiplied
;                  to the normal function before SVD decomposition
;   silent     - Set to suppress print and splog outputs
;   EXTRA      - Keywords passed to either the function FUNC, or DJS_REJECT().
;                Note that keywords like MAXREJ relate to each individual trace.
;
; OUTPUTS:
;   tset       - Structure containing trace set
;
; OPTIONAL OUTPUTS:
;   outmask    - Mask set to 1 for good points and 0 for rejected points;
;                same dimensions as XPOS, YPOS.
;   yfit       - Fit values at each XPOS.
;
; COMMENTS:
;   The fits are done to one trace at a time, where each trace is treated
;   completely independently.
;
;   Note that both MAXDEV and MAXSIG can be set for applying both rejection
;   schemes at once.
;
;   Additional keywords can be passed to the fitting functions with _EXTRA.
;   By not setting any of these rejection keywords, no rejection is performed.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_reject()
;   fchebyshev()
;   fchebyshev_split()
;   flegendre()
;   fpoly()
;   func_fit()
;
; REVISION HISTORY:
;   19-May-1999  Written by David Schlegel, Princeton.
;   04-Aug-1999  Added chebyshev option (DJS).
;   02-Sep-2000  Modify to use rejection schemes in DJS_REJECT() (DJS).
;   07-Dec-2000  Added /silent keyword (DPF)
;   10-Jul-2001  Add polynomial option
;-
;------------------------------------------------------------------------------
pro xy2traceset, xpos, ypos, tset, invvar=invvar, func=func, ncoeff=ncoeff, $
 xmin=xmin, xmax=xmax, maxiter=maxiter, inputfunc=inputfunc, $
 inmask=inmask, outmask=outmask, yfit=yfit, inputans=inputans, $
 silent=silent, _EXTRA=EXTRA

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xy2traceset, xpos, ypos, tset, [invvar=, func=, ncoeff=, $'
      print, ' xmin=, xmax=, maxiter=, inmask=, outmask=, yfit=, inputans=, _EXTRA= ]'
      return
   endif

   if (NOT keyword_set(func)) then func = 'legendre'
   if (n_elements(maxiter) EQ 0) then maxiter = 10

   ndim = size(ypos, /n_dim)
   dims = size(ypos, /dim)

   if (ndim EQ 1) then begin
      nx = dims[0]
      ntrace = 1
   endif else if (ndim EQ 2) then begin
      nx = dims[0]
      ntrace = dims[1]
   endif else begin
      message, 'XPOS contains invalid number of dimensions'
   endelse

   case func of
     'poly':     function_name = 'poly'
     'legendre': function_name = 'flegendre'
     'chebyshev': function_name = 'chebyshev'
     'chebyshev_split': function_name = 'chebyshev_split'
     else: message, 'Unknown function' + func
   endcase

   if (NOT keyword_set(ncoeff)) then ncoeff = 3

   tset = $
   { func    :    func              , $
     xmin    :    0.0               , $
     xmax    :    0.0               , $
     coeff   :    dblarr(ncoeff, ntrace) $
   }

   if (size(xmin, /tname) NE 'UNDEFINED') then tset.xmin = xmin $
    else tset.xmin = min(xpos)
   if (size(xmax, /tname) NE 'UNDEFINED') then tset.xmax = xmax $
    else tset.xmax = max(xpos)
   xmid = 0.5 * (tset.xmin + tset.xmax)
   xrange = tset.xmax - tset.xmin

   outmask = bytarr(nx, ntrace)

   yfit = ypos*0.0
   if (NOT keyword_set(inputans)) then curans = fltarr(ncoeff)

   ; Header for Burles counter
   if not keyword_set(silent) then print, ''
   if not keyword_set(silent) then print, ' TRACE# NPOINTS NREJECT'

   ;----------
   ; Loop over each trace

   for itrace=0, ntrace-1 do begin

      xnorm = 2.0 * (xpos[*,itrace] - xmid) / xrange ; X positions renormalized

      if (keyword_set(inputans)) then curans = inputans[*,itrace] 
      if (keyword_set(inmask)) then curinmask = inmask[*,itrace]

      ;----------
      ; Rejection iteration loop

      iiter = 0
      qdone = 0
      ycurfit = 0
      if (keyword_set(invvar)) then tempivar = invvar[*,itrace] $
      else tempivar = bytarr(nx) + 1
      curoutmask = tempivar GT 0

      while (NOT keyword_set(qdone) AND iiter LE maxiter) do begin

         if keyword_set(inputfunc) then $
           res = func_fit(xnorm, ypos[*,itrace], ncoeff, $
             invvar=tempivar*curoutmask, $
             function_name=function_name, yfit=ycurfit, inputans=curans, $
             inputfunc=inputfunc[*,itrace], _EXTRA=EXTRA) $
         else res = func_fit(xnorm, ypos[*,itrace], ncoeff, $
             invvar=tempivar*curoutmask, $
          function_name=function_name, yfit=ycurfit, inputans=curans, $
          _EXTRA=EXTRA)

         curinmask = curoutmask

         qdone = djs_reject(ypos[*,itrace], ycurfit, invvar=tempivar, $
          inmask=curinmask, outmask=curoutmask, _EXTRA=EXTRA)

         iiter = iiter + 1
      endwhile
 
      yfit[*,itrace] = ycurfit 
      tset.coeff[*,itrace] = res
      outmask[*,itrace] = curoutmask

      ; Burles counter of row number...
      junk = where(curoutmask EQ 0, nreject)
      if not keyword_set(silent) then $
        print, format='(i7,i8,i8,a1,$)', itrace, nx, nreject, string(13b)
   endfor

   junk = where(outmask EQ 0, nreject)
   if not keyword_set(silent) then splog, 'Total rejected = ', nreject

   return
end
;------------------------------------------------------------------------------
