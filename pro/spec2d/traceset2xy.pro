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
;   xpos       - X positions corresponding to YPOS as an [nx,Ntrace] array
;   ypos       - Y centers as an [nx,nTrace] array
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

      nx = long(tset.xmax - tset.xmin + 1)

      xmid = 0.5 * (tset.xmin + tset.xmax)
      xrange = tset.xmax - tset.xmin
      xpos = djs_laxisgen([nx, nTrace], iaxis=0) + tset.xmin

      ypos = fltarr(nx, nTrace)
      xvec = 2.0 * findgen(nx)/(nx-1) - 1.0
      legarr = flegendre(xvec, ncoeff)
      for iTrace=0, nTrace-1 do begin
         ypos[*,iTrace] = legarr # tset.coeff[*,iTrace]
      endfor

      end
   else: error, 'Unknown function' + func
   endcase

   return
end
;------------------------------------------------------------------------------
