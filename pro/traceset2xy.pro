;+
; NAME:
;   traceset2xy
;
; PURPOSE:
;   Convert from a trace set to an array of x,y positions
;
; CALLING SEQUENCE:
;   traceset2xy, tset, xpos, ypos
;
; INPUTS:
;   tset       - Structure containing trace set
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   xpos       - X centers as an [ny,nTrace] array
;   ypos       - Y positions corresponding to XPOS as an [ny,Ntrace] array
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
pro traceset2xy, tset, xpos, ypos

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - traceset2xy, tset, xpos, ypos'
      return
   endif

   case tset.func of
   'legendre': begin
      ncoeff = (size(tset.coeff))[1]
      nTrace = (size(tset.coeff))[2]
      ny = long(tset.ymax - tset.ymin + 1)

      ymid = 0.5 * (tset.ymin + tset.ymax)
      yrange = tset.ymax - tset.ymin

      ypos = djs_laxisgen([ny, nTrace], iaxis=0)
      xpos = fltarr(ny, nTrace)
      ypos = 2.0 * (ypos - ymid) / yrange

      for iTrace=0, nTrace-1 do begin
         xpos[*,iTrace] = flegendre(ypos[*,iTrace], ncoeff) $
          # tset.coeff[*,iTrace]
;         for iy=0, ny-1 do xpos[iy,iTrace] = $
;          transpose(tset.coeff[*,iTrace]) # svdleg(ypos[iy,iTrace], ncoeff)
      endfor

      end
   else: error, 'Unknown function' + func
   endcase

   return
end
;------------------------------------------------------------------------------
