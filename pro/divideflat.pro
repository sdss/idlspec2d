pro divideflat,  flux, fluxivar, fflat, fibermask

    npix = (size(flux))[1] 
    nTrace = (size(flux))[2] 

    if((size(fluxivar))[1] NE npix OR (size(fluxivar))[2] NE nTrace) then $
        message, 'flux and fluxivar are not the same size'

    if((size(fflat))[1] NE npix OR (size(fflat))[2] NE nTrace) then $
        message, 'flux and fflat are not the same size'

    if((size(fibermask))[1] NE nTrace) then $
        message, 'flux and fibermask have different number of fibers'
    
    for i=0,nTrace-1 do begin
      if (fibermask[i]) then begin
;
;	Find where flat fields might be trouble
;
        badflat = where(fflat[*,i] LE 0.0)
        goodflat = where(fflat[*,i] GT 0.0)
        if (goodflat[0] EQ -1) then begin
          print, ' No good flat field points in trace ', i
          fluxivar[*,i] = 0.0
          flux[*,i] = 0.0
        endif else begin
          if (badflat[0] NE -1) then begin
            fluxivar[badflat,i] = 0.0
            flux[badflat,i] = 0.0
          endif
;
;	Only divide good fflat pixels
;
          flux[goodflat,i] = flux[goodflat,i] / fflat[goodflat,i]
          fluxivar[goodflat,i] = fluxivar[goodflat,i] * fflat[goodflat,i]^2
        endelse
      endif else begin
;
;	Set both to zero for now
;
        fluxivar[*,i] = 0.0
        flux[*,i] = 0.0
      endelse
    endfor
end
          

