;+
; NAME:
;	ushift
; PURPOSE: (one line)
;	Shift data using a damped sinc onto slightly non-linear x values
;	
; DESCRIPTION:
;
;	This function will shift an array of data pointed to by x and
;	extending for n points.  The amount of the shift is given by shift.
;	The result of the operation is placed at xp.  A shift that is within
;	0.0001 of a whole number is treated to be that of the whole number.  If
;	the shift is by an integral number of pixels then the shift involves
;	reindexing the data, no interpolation is done.  If the shift is some
;	non-integral amount then the data is resampled using a damped sinc
;	function.
;
;	The sense of the shift is as follows: think of the array plotted on a
;	fixed scale.  A shift of 1 corresponds to shifting the data by one data
;	point to the right relative to the fixed scale, ie. x[3]=xp[4].
;
;	The data will fall off one end or another of the output vector as a
;	result of the shift.  However, this is not as significant as the edge
;	effect, the convolution is not complete for any data point within 10
;	points of the edge, so those points cannot be trusted.  The missing
;	points in the convolution are assumed to be equal to the end points.
;
; CATEGORY:
;       Numerical
; CALLING SEQUENCE:
;	yp = ushift(y, invset, xp)
; INPUTS:
;	y     - Input data array to be shifted.
;	invset - Structure which contains coefficients to map xin to pixels
;       xin    - The x-axis values for the requested yp
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	Return value is the shifted array.
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
;	The input and output arrays cannot be the same.
; PROCEDURE:
; MODIFICATION HISTORY:
;	Adapted from Zodiac routine: shiftc/sshift
;	  Marc W. Buie, Lowell Observatory, 1992 October 2
;-

function ushift,data,invset,xin, function_name=function_name


   EPS     = 1.0e-5 ; Smallest fractional shift allowed.
   DAMPFAC = 3.65   ; Damping factor for gaussian.
   NS      = 21     ; Number of points in the sinc convolution kernal.

   if (NOT keyword_set(function_name)) then function_name='legendre'

   npts = n_elements(xin)
   ndata   = n_elements(data)
   yp   = fltarr(npts)
   buffer = 300
   pixels = fltarr(npts+2*buffer)
   pad = fltarr(ndata+2*buffer)
   ncoeff = (size(invset))[1]-2
   


   xnorm = (2.0*xin - (invset[1]+invset[0]))/(invset[1]-invset[0])
   if (function_name EQ 'legendre') then $
     pixels[buffer:npts+buffer-1] = $
       flegendre(xnorm, ncoeff) # invset[2:*]
   if (function_name EQ 'chebyshev') then $
     pixels[buffer:npts+buffer-1] = $
       fchebyshev(xnorm, ncoeff) # invset[2:*]

;
;	Put in some unit step buffer space at the beginning and end.
;
   pixels[0:buffer] = pixels[buffer] - reverse(findgen(buffer+1))
   pixels[npts+buffer-1:*] = pixels[npts+buffer-1] + findgen(buffer+1)

   pad[buffer:ndata+buffer-1] = data
   pad[0:buffer-1] = data[0]
   pad[ndata+buffer:*] = data[ndata-1]


; Convolve the sinc array with the input data.  This is the shift.
   for point=0,npts-1 do begin
     curpixel = pixels[point+buffer]
     spot = fix(curpixel)
     fshift = curpixel - spot
; Do the fractional shift first (if necessary).
     if ( ( abs(fshift) gt EPS ) and ( abs(fshift) lt 1.0-EPS) ) then begin

       first = point - NS/2 + buffer
       last = point + NS/2 + buffer
       y = spot - pixels[first:last]
       py = !pi * y
       sinc = exp( -y^2/DAMPFAC^2 ) * sin(py)/py
       sinc = rotate(sinc,2)

       lobe = (indgen(NS) - 10) + spot + buffer
       vals = pad[lobe]

       yp[point] = total(sinc*vals)/total(sinc)


; No fractional shift, just copy the data.
     endif else begin

      yp[point] = pad[spot+buffer]

     endelse
   endfor

   return,yp

end
