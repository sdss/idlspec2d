;+
; NAME:
;   xy2traceset
;
; PURPOSE:
;   Convert from an array of x,y positions to a trace set
;
; CALLING SEQUENCE:
;   xy2traceset, xpos, ypos, tset, [ func=func, ncoeff=ncoeff, $
;    xmin=xmin, xmax=xmax, maxdev=maxdev, maxiter=maxiter, $
;    singlerej=singlerej, xmask=xmask ]
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
;   xmin       - Explicitly set XMIN for trace set rather than using minimum
;                in XPOS
;   xmax       - Explicitly set XMAX for trace set rather than using maximum
;                in XPOS
;   maxdev     - Maximum deviation in the fit to YPOS; default to rejecting any
;                points that deviate by more than 1.0 from the fit.
;   maxiter    - Maximum number of rejection iterations; default to 10;
;                set to 0 for no rejection.
;                Rejection iterations continues until convergence
;                (actually, until the number of rejected pixels is unchanged).
;   singlerej  - If set, then reject at most one deviant point per iteration,
;                rejecting the worst one each time.  In this case, MAXITER
;                represents the maximum number of points that can be rejected
;                per trace.
;
; OUTPUTS:
;   tset       - Structure containing trace set
;
; OPTIONAL OUTPUTS:
;   xmask      - Mask set to 1 for good points and 0 for rejected points;
;                same dimensions as XPOS, YPOS.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Should probably change default to no rejection.
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
 xmin=xmin, xmax=xmax, maxdev=maxdev, maxiter=maxiter, $
 singlerej=singlerej, xmask=xmask

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xy2traceset, xpos, ypos, tset, [ func=, ncoeff=, $
      print, ' xmin=, xmax=, maxdev=, maxiter=, /singlerej, xmask= ]'
      return
   endif

   if (NOT keyword_set(func)) then func = 'legendre'
   if (NOT keyword_set(maxdev)) then maxdev = 1.0
   if (N_elements(maxiter) EQ 0) then maxiter = 10

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

      for itrace=0, ntrace-1 do begin
;         res = svdfit(2.0*(xpos[*,i]-xmid)/xrange, ypos[*,itrace], ncoeff, $
;          /double, function_name=function_name, singular=singular)

         xnorm = 2.0*(xpos[*,itrace]-xmid) / xrange ; X positions renormalized
         nreject = 1
         totalreject = 0
         good = lonarr(nx) + 1

         ; Rejection iteration loop

         iiter = 0
         ngood = nx
         nglast = nx+1 ; Set to anything other than NGOOD for 1st iteration
         while ( (iiter EQ 0 AND maxiter EQ 0) $
              OR (ngood NE nglast AND iiter LE maxiter AND ngood GT 1) $
          ) do begin

            if (iiter EQ 0) then begin
               qgood = bytarr(nx) + 1
               igood = lindgen(nx)
            endif else begin
               nglast = ngood
               if (keyword_set(singlerej)) then begin
                  worstdiff = max(abs(yfit[igood]-ypos[igood,itrace]), iworst)
                  if (worstdiff GT maxdev) then begin
                     qgood[igood[iworst]] = 0
                     igood = where(qgood, ngood)
                  endif
               endif else begin
                  qgood = abs(yfit-ypos[*,itrace]) LT maxdev
                  igood = where(qgood, ngood)
               endelse
            endelse

            totalreject = nx - ngood
            res = func_fit(xnorm[igood], ypos[igood,itrace], ncoeff, $
              function_name=function_name)

            if (func EQ 'legendre') then $
               yfit = flegendre(xnorm, ncoeff) # res
            if (func EQ 'chebyshev') then $
               yfit = fchebyshev(xnorm, ncoeff) # res

            iiter = iiter + 1
         endwhile
   
         tset.coeff[*,itrace] = res
         if (totalreject GT 0) then $
          print, 'Rejected ', totalreject, ' of ', nx, $
           ' points on trace ', itrace

         xmask[*,itrace] = qgood

         ; Burles counter of row number...
         print, format='($, ".",i4.4,a5)', itrace, string([8b,8b,8b,8b,8b])
      endfor

   endif else begin
      message, 'Unknown function' + func
   endelse

   return
end
;------------------------------------------------------------------------------
