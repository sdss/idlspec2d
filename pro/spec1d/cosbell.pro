;Return a cosbell window function for subsequent applicaiton to data
;and variance

function cosbell, npix, fraction

  if N_params() LT 2 then return, -1
  
  window = fltarr(npix)+1.0

  nfirst = lindgen(fix(npix * fraction))
  nlast = npix - lindgen(fix(npix * fraction)) - 1
  
  window[nfirst] = 0.5 * (1.0 - cos(!Pi*nfirst/(n_elements(nfirst)-1)))
  window[nlast] = 0.5 * (1.0 - cos(!Pi*nfirst/(n_elements(nfirst)-1)))
  
  return, window
end       
        

