;+
; NAME:
;   sdss_spec_block
; PURPOSE:
;   take a set of sdss spectra and output as a big block
; CALLING SEQUENCE:
;   sdss_spec_block
; INPUTS:
;   plate - list of plates
;   fiberid - list of fibers
;   mjd - list of mjds
; OPTIONAL INPUTS:
; KEYWORDS:
;   masksky     mask out sky line regions in all spectra (only do this if
;                 the spectra span a wide redshift range)
; OUTPUTS:
;   avlum       on output, average spectrum (lambda L_lambda units)
; DEPENDENCIES:
;   idlutils
;   idlspec2d
; BUGS:
;   Input and output wavelength grids MUST be logarithmically spaced.
;   User has to worry about units.
; REVISION HISTORY:
;   2001-11-17  written - Hogg (NYU)
;-
pro sdss_spec_block, plate, fiberid, mjd, block_flux=block_flux, $
                     block_ivar=block_ivar, block_lambda=block_lambda

masksky=1

; set defaults
nspec= n_elements(plate)
if(n_elements(avloglam) eq 0) then $
  avloglam=double(alog10(3500.)+(alog10(9500.)-alog10(3500.))* $
                  (dindgen(4000)+0.5)/4000.)
nlam= n_elements(avloglam)
if(n_elements(avlum) eq 0) then avlum=dblarr(nlam)
startavlum= avlum

; group by plate and get all spectra
platelist=plate[uniq(plate,sort(plate))]
avlum= 0d
isp=0L
for iplate=0L, n_elements(platelist)-1L do begin
    splog,'iplate= '+string(iplate)

    plate_indx=where(plate eq platelist[iplate],plate_count)
    readspec, plate[plate_indx],fiberid[plate_indx],mjd=mjd[plate_indx], $
      flux=flux,flerr=flerr,invvar=invvar,andmask=andmask,ormask=ormask, $
      loglam=loglam, zans=zans
    flux= double(flux)
    invvar= double(invvar)
    for i=0L,plate_count-1 do begin
        splog,'i= '+string(i)

; read spectrum and inverse variance
        if keyword_set(masksky) then $
          invvar[*,i]= skymask(invvar[*,i], andmask[*,i], ormask[*,i])
        ;flux[*,i]= flux[*,i]*(10.0^loglam[*,i])
        ;invvar[*,i]= invvar[*,i]/(10.0^loglam[*,i])
        
; shift wavelength scale to rest-frame
        restloglam= loglam[*,i]-alog10(1.d + zans[i].z)

; slam onto new grid (MUST BE LOGARITHMIC SPACING)
        lum=fltarr(n_elements(restloglam))
        combine1fiber, restloglam,flux[*,i],invvar[*,i], $
          newloglam=avloglam, newflux=tmp1, newivar=tmp2, maxiter=0
        lum= double(tmp1)
        lum_ivar= double(tmp2)

; increment vectors
        if(n_elements(block_flux) eq 0) then begin
            nl=n_elements(lum)
            block_flux=fltarr(nl,nspec)
            block_ivar=fltarr(nl,nspec)
            block_lambda=fltarr(nl)
        endif 
        block_flux[*,isp]=lum
        block_ivar[*,isp]=lum_ivar
        block_lambda[*]=10.^(avloglam)
        isp=isp+1L
    endfor
endfor
    
return
end
