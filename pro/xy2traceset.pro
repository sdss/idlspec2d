;+
; NAME:
;   xy2traceset
;
; PURPOSE:
;   Convert from an array of x,y positions to a trace set
;
; CALLING SEQUENCE:
;   xy2traceset, xpos, ypos, tset, [func=func, ncoeff=ncoeff]
;
; INPUTS:
;   xpos       - X positions corresponding to YPOS as an [nx,Ntrace] array
;   ypos       - Y centers as an [nx,nTrace] array
;
; OPTIONAL KEYWORDS:
;   func       - Function for trace set; options are:
;                'legendre'
;                'chebyshev'
;                Default to 'legendre'
;   ncoeff     - Number of coefficients in fit; default to 3
;   xmin       - Explicitly set XMIN for trace set rather than using minimum
;                in XPOS
;   xmax       - Explicitly set XMAX for trace set rather than using maximum
;                in XPOS
;   maxdev     - Maximum deviation of X in pixels; default to rejecting any
;                XPOS positions that deviate by more than 1.0 pixels from
;                the fit.  Rejection iterations continues until convergence
;                (actually, until the number of rejected pixels is unchanged).
;
; OUTPUTS:
;   tset       - Structure containing trace set
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   flegendre()
;   fchebyshev()
;
; REVISION HISTORY:
;   19-May-1999  Written by David Schlegel, Princeton.
;   04-Aug-1999  Added chebyshev option (DJS).
;-
;------------------------------------------------------------------------------
pro xy2traceset, xpos, ypos, tset, func=func, ncoeff=ncoeff, $
 xmin=xmin, xmax=xmax, maxdev=maxdev

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xy2traceset, xpos, ypos, tset, [func=func, ncoeff=ncoeff, $
      print, ' xmin=xmin, xmax=xmax]'
      return
   endif

   if (NOT keyword_set(func)) then func = 'legendre'
   if (NOT keyword_set(maxdev)) then maxdev = 1.0

   ndim = size(ypos, /n_dim)
   dims = size(ypos, /dim)

   if (ndim EQ 1) then begin
      nx = dims[0]
      nTrace = 1
   endif else if (ndim EQ 2) then begin
      nx = dims[0]
      nTrace = dims[1]
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
        coeff   :    fltarr(ncoeff, nTrace) $
      }

      if (size(xmin, /tname) NE 'UNDEFINED') then tset.xmin = xmin $
       else tset.xmin = min(xpos)
      if (size(xmax, /tname) NE 'UNDEFINED') then tset.xmax = xmax $
       else tset.xmax = max(xpos)
      xmid = 0.5 * (tset.xmin + tset.xmax)
      xrange = tset.xmax - tset.xmin

      for i=0, nTrace-1 do begin
;         res = svdfit(2.0*(xpos[*,i]-xmid)/xrange, ypos[*,i], ncoeff, $
;          /double, function_name=function_name, singular=singular)

         xnorm = 2.0*(xpos[*,i]-xmid) / xrange ; X positions renormalized
         nreject = 1
         totalreject = 0
         good = lonarr(nx) + 1

         ; Rejection iteration loop

         igood = lindgen(nx)
         ngood = nx
         nglast = nx+1 ; Set to anything other than NGOOD for 1st iteration
         while (ngood NE nglast) do begin
            res = func_fit(xnorm[igood], ypos[igood,i], ncoeff, $
              function_name=function_name)

            if (func EQ 'legendre') then $
               yfit = flegendre(xnorm[igood], ncoeff) # res
            if (func EQ 'chebyshev') then $
               yfit = fchebyshev(xnorm[igood], ncoeff) # res

            nglast = ngood
            igood = where( abs(yfit-ypos[igood,i]) LT maxdev, ngood)
         endwhile    
   
         tset.coeff[*,i] = res
         if (totalreject GT 0) then $
          print, 'Rejected ', totalreject, ' pixels on trace ', i

         ; Burles counter of row number...
         print, format='($, ".",i4.4,a5)', i, string([8b,8b,8b,8b,8b])
      endfor

   endif else begin
      error, 'Unknown function' + func
   endelse

   return
end
;------------------------------------------------------------------------------
