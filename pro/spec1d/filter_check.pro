;+
; NAME:
;   filter_check
;
; PURPOSE:
;   
;   Fit for the filter curves by comparing photometry and spectrophotometry
;
; CALLING SEQUENCE:
;   
;   filter_check,spallfile,binboundsfile
;
; INPUTS:
;   spallfile      - spectro info
;
; OPTIONAL INPUTS:
;
;   binboundsfile  - bin boundaries 
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   binRfile   - file containing output stuff
;
; COMMENTS:
;
; BUGS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; DATA FILES:
;
; REVISION HISTORY:
;   05-APr-2000  Written by M. Blanton, Fermiland
;-
;------------------------------------------------------------------------------
pro filter_check, spselect, flux,invvar,loglam, binRbase, $ 
binbounds=binbounds, nosubtract=nosubtract

  ;--------
  ; Solve in each band
  nselect=n_elements(spselect)
  for band=1, 3 do begin

    ;--------
    ; Create photocounts
    photomags=spselect[*].psfCounts[band] $
        +(spselect[*].fiberCounts[2]-spselect[*].psfCounts[2])
    photocounts=10.^(-0.4*photomags)-spselect[*].counts_spectro*0.79

    ;--------
    ; Read in bin boundaries, bin spectra if desired
    if(keyword_set(binbounds)) then begin
      ;filename=binboundsbase+'.'+strtrim(string(i),2)+'.dat'
       ;openr,11,filename
       ;readf,nbins
       ;binbounds=lonarr(nbins+1)
       ;readf,binbounds
       ;close,11
       nbins=n_elements(binbounds)-1
       binflux=fltarr(nbins+1-keyword_set(nosubtract),nselect)
       bin_spectra,flux,invvar,binbounds,binflux=binflux
    endif else begin
       if(keyword_set(nosubtract)) then return
       nbins=0
       binflux=fltarr(1,nselect)
    endelse 

    ;--------
    ; Flux convolved with bins
    if(not keyword_set(nosubtract)) then $
       binflux[nbins,*]=spselect[*].counts_spectro[band]

    ;-------
    ; Iterate
    indx=lindgen(n_elements(spselect))
    for iter=0, 1 do begin

      ;--------
      ; Find solution 
      filter_solve,binflux[*,indx],photocounts[indx],binR

      ;-------
      ; Sigma clip
      resid=transpose(binflux##binR)-photocounts
      rms=djsig(resid[indx])
      help,rms
      mean=djs_mean(resid[indx])
      help,mean
      indx=where((resid/rms)^2 lt 100.)

    endfor
    print, binR
erase & plot, 10.0^(0.5*(loglam[binbounds[0:nbins-1]]+loglam[binbounds[1:nbins]])),binR,psym=10
stop
erase & plot, spselect.psfcounts[1]-spselect.psfcounts[2], resid, psym=1
djs_oplot, photocounts[indx], resid[indx], psym=1, color='red'

stop

    ;--------
    ; Output solution 
    filename=binRbase+'.'+strtrim(string(band),2)+'.dat'
    openw,11,filename
    printf,11,binR
    close,11

  endfor 

end
;------------------------------------------------------------------------------
