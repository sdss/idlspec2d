
;
;	Pad input array with mean with at least npix
;
;

function fftpad, data, direction, _EXTRA=KeywordsForFft

     if (N_elements(direction) EQ 0) then direction = -1

     npix = n_elements(data)
     nextpower = fix(alog(2*npix) / alog(2.0)) + 1

     if (nextpower GT 31) then begin
       print, 'too many elements'
       return, -1
     endif

     nfinalpix = 2L^nextpower
     paddata = make_array(nfinalpix,type=size(data,/type), value = mean(data))

     first = (nfinalpix - npix)/2
     paddata[first:first+npix-1] = data


     ;
     ;	Here we need to apodize the paddata to get rid of ringing...
     ;	Skipping for now (laziness)


     return, fft(paddata, direction, _Extra=KeywordsForFFT)
end

