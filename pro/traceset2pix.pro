;+
; NAME:
;   traceset2pix
;
; PURPOSE:
;   Use a traceset to convert lambda to pix
;
; CALLING SEQUENCE:
;   pix = traceset2pix(invset, lambda)
;
; INPUTS:
;   invset     - Structure containing trace set
;   lambda     - wavelengths at which to evaluate
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
;   djs_laxisgen()
;   flegendre()
;
; REVISION HISTORY:
;   19-May-1999  Written by David Schlegel, Princeton.
;   14-Oct-1999  D. Finkbeiner
;-
;------------------------------------------------------------------------------
function traceset2pix, invset, lambda

   ; Need 3 parameters
   if (N_params() LT 2) then begin
      print, 'Syntax - traceset2xy, invset, xpos, ypos'
      return, -1
   endif

   if (invset.func EQ 'legendre' OR invset.func EQ 'chebyshev') then begin

      ndim = size(invset.coeff, /n_dim)
      dims = size(invset.coeff, /dim)

      if (ndim EQ 1) then begin
         ncoeff = dims[0]
         nTrace = 1
      endif else if (ndim EQ 2) then begin
         ncoeff = dims[0]
         nTrace = dims[1]
      endif else begin
         message, 'INVSET.COEFF contains invalid number of dimensions'
      endelse

      nx = long(invset.xmax - invset.xmin + 1)


	nrows = 320
        nlines = n_elements(lambda)
	xnorm = (2.0d*lambda - (invset.xmin + invset.xmax))/ $
	          (invset.xmax - invset.xmin)

	ysky = dblarr(nrows,nlines)
	pix = dblarr(nrows,nlines)
	for i=0,nrows-1 do begin
	  ysky[i,*] = float(i)
	  if (invset.func EQ 'legendre') then $
              pix[i,*] = flegendre(xnorm,ncoeff) # invset.coeff[*,i]
	  if (invset.func EQ 'chebyshev') then $
              pix[i,*] = fchebyshev(xnorm,ncoeff) # invset.coeff[*,i]
 
	endfor	  

   endif else begin
      error, 'Unknown function' + func
   endelse

   return, pix
end
;------------------------------------------------------------------------------
