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
pro xy2traceset, xpos, ypos, tset, func=func, ncoeff=ncoeff

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xy2traceset, xpos, ypos, tset, [func=func, ncoeff=ncoeff]'
      return
   endif

   if (NOT keyword_set(func)) then func = 'legendre'

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

      tset.xmin = min(xpos)
      tset.xmax = max(xpos)
      xmid = 0.5 * (tset.xmin + tset.xmax)
      xrange = tset.xmax - tset.xmin

      for i=0, nTrace-1 do begin
;         res = svdfit(2.0*(xpos[*,i]-xmid)/xrange, ypos[*,i], ncoeff, $
;          /double, /legendre, singular=singular)
         res = svdfit(2.0*(xpos[*,i]-xmid)/xrange, ypos[*,i], ncoeff, $
          /double, function_name=function_name, singular=singular)
         tset.coeff[*,i] = res
      endfor

   endif else begin
      error, 'Unknown function' + func
   endelse

   return
end
;------------------------------------------------------------------------------
