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
;   xpos       - X centers as an [ny,nTrace] array
;   ypos       - Y positions corresponding to XPOS as an [ny,Ntrace] array
;
; OPTIONAL KEYWORDS:
;   func       - Function for trace set; options are:
;                'legendre'
;   ncoeff     - Number of coefficients in fit
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
;
; REVISION HISTORY:
;   19-May-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro xy2traceset, xpos, ypos, tset, func=func, ncoeff=ncoeff

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xy2traceset, xpos, ypos, tset, [func=func, ncoeff=ncoeff]'
      return
   endif

   if (NOT keyword_set(func)) then func = 'legendre'

   nTrace = (size(xpos))[2]
   ny = (size(xpos))[3]

   case func of
   'legendre': begin
      if (NOT keyword_set(ncoeff)) then ncoeff = 3

      tset = $
      { func    :    'legendre'        , $
        ymin    :    0.0               , $
        ymax    :    0.0               , $
        coeff   :    fltarr(ncoeff, nTrace) $
      }

      tset.ymin = min(ypos)
      tset.ymax = max(ypos)
      ymid = 0.5 * (tset.ymin + tset.ymax)
      yrange = tset.ymax - tset.ymin

      for i=0, nTrace-1 do begin
;         res = svdfit(2.0*(ypos[*,i]-ymid)/yrange, xpos[*,i], ncoeff, $
;          /double, /legendre, singular=singular)
         res = svdfit(2.0*(ypos[*,i]-ymid)/yrange, xpos[*,i], ncoeff, $
          /double, function_name='flegendre', singular=singular)
         tset.coeff[*,i] = res
      endfor

      end
   else: error, 'Unknown function' + func
   endcase

   return
end
;------------------------------------------------------------------------------
