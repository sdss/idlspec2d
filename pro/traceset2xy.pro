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
      ndim = size(tset.coeff, /n_dim)
      dims = size(tset.coeff, /dim)

      if (ndim EQ 1) then begin
         ncoeff = dims[0]
         nTrace = 1
      endif else if (ndim EQ 2) then begin
         ncoeff = dims[0]
         nTrace = dims[1]
      endif else begin
         message, 'TSET.COEFF contains invalid number of dimensions'
      endelse

      ny = long(tset.ymax - tset.ymin + 1)

      ymid = 0.5 * (tset.ymin + tset.ymax)
      yrange = tset.ymax - tset.ymin
      ypos = djs_laxisgen([ny, nTrace], iaxis=0) + tset.ymin

      xpos = fltarr(ny, nTrace)
      yvec = 2.0 * findgen(ny)/(ny-1) - 1.0
      legarr = flegendre(yvec, ncoeff)
      for iTrace=0, nTrace-1 do begin
         xpos[*,iTrace] = legarr # tset.coeff[*,iTrace]
      endfor

      end
   else: error, 'Unknown function' + func
   endcase

   return
end
;------------------------------------------------------------------------------
