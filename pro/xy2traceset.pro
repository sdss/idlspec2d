;+
; NAME:
;   xy2traceset
;
; PURPOSE:
;   Convert from an array of x,y positions to a trace set
;
; CALLING SEQUENCE:
;   xy2traceset, xpos, ypos, tset, [ func=func, ncoeff=ncoeff, $
;    xmin=xmin, xmax=xmax, maxdev=maxdev, maxsig=maxsig, maxiter=maxiter, $
;    singlerej=singlerej, xmask=xmask, yfit=yfit, inputans=inputans, $
;    invvar=invvar, _EXTRA=KeywordsForFuncFit ]
;
; INPUTS:
;   xpos       - X positions corresponding to YPOS as an [nx,Ntrace] array
;   ypos       - Y centers as an [nx,ntrace] array
;
; OPTIONAL KEYWORDS:
;   func       - Function for trace set; options are:
;                'legendre'
;                'chebyshev'
;                Default to 'legendre'
;   ncoeff     - Number of coefficients in fit; default to 3
;   ncoeff     - Inverse variance for weighted func_fit
;   xmin       - Explicitly set XMIN for trace set rather than using minimum
;                in XPOS
;   xmax       - Explicitly set XMAX for trace set rather than using maximum
;                in XPOS
;   maxdev     - Maximum deviation in the fit to YPOS; set to 0 for no reject;
;                default to 0.
;   maxsig     - Maximum deviation in the fit to YPOS in terms of the 1-sigma
;                dispersion of the residuals; set to 0 for no reject;
;                default to 0.
;   maxiter    - Maximum number of rejection iterations; set to 0 for no
;                rejection; default to 10 if either MAXDEV or MAXSIG are set.
;                Rejection iterations continues until convergence
;                (actually, until the number of rejected pixels is unchanged).
;   singlerej  - If set, then reject at most one deviant point per iteration,
;                rejecting the worst one each time.  In this case, MAXITER
;                represents the maximum number of points that can be rejected
;                per trace.
;   inputans   - ???
;
; OUTPUTS:
;   tset       - Structure containing trace set
;
; OPTIONAL OUTPUTS:
;   xmask      - Mask set to 1 for good points and 0 for rejected points;
;                same dimensions as XPOS, YPOS.
;   yfit       - Fit values at each XPOS.
;
; COMMENTS:
;   Note that both MAXDEV and MAXSIG can be set for applying both rejection
;   schemes at once.
;
;   The HALFINTWO keyword can be passed to FCHEBYSHEV by this procedure.
;
; EXAMPLES:
;
; BUGS:
;   Should probably change default to no rejection.
;
; PROCEDURES CALLED:
;   fchebyshev()
;   flegendre()
;   func_fit()
;
; REVISION HISTORY:
;   19-May-1999  Written by David Schlegel, Princeton.
;   04-Aug-1999  Added chebyshev option (DJS).
;-
;------------------------------------------------------------------------------
pro xy2traceset, xpos, ypos, tset, func=func, ncoeff=ncoeff, $
 xmin=xmin, xmax=xmax, maxdev=maxdev, maxsig=maxsig, maxiter=maxiter, $
 singlerej=singlerej, xmask=xmask, yfit=yfit, inputans=inputans, $
 invvar=invvar, _EXTRA=KeywordsForFuncFit

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xy2traceset, xpos, ypos, tset, [ func=, ncoeff=, $
      print, ' xmin=, xmax=, maxdev=, maxsig=, maxiter=, /singlerej, yfit=, xmask= ]'
      return
   endif

   if (NOT keyword_set(func)) then func = 'legendre'
   if (N_elements(maxdev) EQ 0) then maxdev = 0
   if (N_elements(maxsig) EQ 0) then maxsig = 0
   if (N_elements(maxiter) EQ 0) then begin
      if (maxdev NE 0 or maxsig NE 0) then numiter = 10 $
       else numiter = 0
   endif else begin
      numiter = maxiter
   endelse

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

   if (func EQ 'legendre' OR func EQ 'chebyshev') then begin
      if (func EQ 'legendre') then function_name = 'flegendre'
      if (func EQ 'chebyshev') then function_name = 'fchebyshev'

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

      xmask = bytarr(nx, ntrace)

      yfit = ypos*0.0
      if (NOT keyword_set(inputans)) then curans = fltarr(ncoeff)

      totalreject = 0
      ; Header for Schlegel counter
      print, ''
      print, ' TRACE# NPOINTS NREJECT'

      for itrace=0, ntrace-1 do begin
;         res = svdfit(2.0*(xpos[*,i]-xmid)/xrange, ypos[*,itrace], ncoeff, $
;          /double, function_name=function_name, singular=singular)

         xnorm = 2.0*(xpos[*,itrace]-xmid) / xrange ; X positions renormalized
         nreject = 0
         good = lonarr(nx) + 1

         if (keyword_set(inputans)) then curans = inputans[*,itrace] 

         ; Rejection iteration loop

         iiter = 0
         ngood = nx
         nglast = nx+1 ; Set to anything other than NGOOD for 1st iteration
         while ( (iiter EQ 0 AND numiter EQ 0) $
              OR (ngood NE nglast AND iiter LE numiter AND ngood GE 1) $
          ) do begin

            if (iiter EQ 0) then begin
               qgood = bytarr(nx) + 1
               igood = lindgen(nx)
            endif else begin
               nglast = ngood
               if (keyword_set(singlerej)) then begin
                  ydiff = ycurfit[igood] - ypos[igood,itrace]
                  if (keyword_set(invvar)) then $
                    invsig = sqrt(invvar[igood,itrace] > 0) $
                  else invsig = (fltarr(ngood) + 1.0)/ stddev(ydiff)
                  worstdiff = max(abs(ydiff), iworst)
                  if (keyword_set(maxdev) AND worstdiff GT maxdev) then $
                     qgood[igood[iworst]] = 0 $
                  else begin
                    worstsig = max(abs(ydiff)*invsig, iworst)
                    if (keyword_set(maxsig) AND worstsig GT maxsig) then $
                     qgood[igood[iworst]] = 0 
                  endelse
               endif else begin
                  ydiff = ycurfit - ypos[*,itrace]
                  if (keyword_set(invvar)) then $
                    invsig = sqrt(invvar[*,itrace] > 0) $
                  else invsig = (fltarr(nx) + 1.0)/ stddev(ydiff)
                  qgood = bytarr(nx) + 1
                  if (keyword_set(maxdev)) then $
                   qgood = qgood AND (abs(ydiff) LT maxdev)
                  if (keyword_set(maxsig)) then $
                   qgood = qgood AND (abs(ydiff)*invsig LT maxsig)
               endelse
               igood = where(qgood, ngood)
            endelse

            nreject = nx - ngood
            if (keyword_set(invvar)) then $
              tempivar  = invvar[*,itrace]*qgood $
            else tempivar = qgood

            res = func_fit(xnorm, ypos[*,itrace], ncoeff, $
             invvar=tempivar, $
             function_name=function_name, yfit=ycurfit, inputans=curans, $
             _EXTRA=KeywordsForFuncFit)

;            if (func EQ 'legendre') then $
;             ycurfit = flegendre(xnorm, ncoeff) # res
;            if (func EQ 'chebyshev') then $
;             ycurfit = fchebyshev(xnorm, ncoeff, _EXTRA=KeywordsForFuncFit) # res

            iiter = iiter + 1
         endwhile
  
         yfit[*,itrace] = ycurfit 
         tset.coeff[*,itrace] = res

         xmask[*,itrace] = qgood

         ; Schlegel counter of row number...
         print, format='(i7,i8,i8,a1,$)', itrace, nx, nreject, string(13b)
         totalreject = totalreject + nreject
      endfor

      splog, 'Total rejected: ', totalreject
      print, ''

   endif else begin
      message, 'Unknown function' + func
   endelse

   return
end
;------------------------------------------------------------------------------
