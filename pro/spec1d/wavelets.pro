function wavelets, data, nscale

     if (NOT keyword_set(nscale)) then nscale = 7

     ;
     ;  Return the array of wavelet deconvolutions
     ;  The sum should exactly equal the original data
     ;

     ; kernel is gaussian with sigma = 1 pixel
 
     nkernel = 7
     nkernelodd =  nkernel/2 
     halfkernel = (errorf((findgen(nkernelodd)+1.5)/sqrt(2.0)) - $
                    errorf((findgen(nkernelodd)+0.5)/sqrt(2.0)))/2.0

     kernel = [reverse(halfkernel),errorf(0.5/sqrt(2.0)),halfkernel]

     ;  Above doesn't work, as kernel needs to expand with each scale

     npix = n_elements(data)

     wavelet = fltarr(npix,nscale+1)

     smooth = data 
     for i=0,nscale - 1 do begin
        sigma = 2.0^i
        nkernelodd = fix(3.0*sigma)
        halfkernel = errorf((findgen(nkernelodd+1)+0.5)/sqrt(2.0 * sigma))
        halfkernel = (halfkernel[1:nkernelodd] - halfkernel[0:nkernelodd-1])/2.0
        kernel = [reverse(halfkernel),errorf(0.5/sqrt(2.0*sigma)),halfkernel]

        temp = convol(smooth,kernel,/EDGE_TRUNCATE)
        wavelet[*,i] = smooth - temp
        smooth = temp
     endfor

     wavelet[*,nscale] = smooth 
     return, wavelet
end 
