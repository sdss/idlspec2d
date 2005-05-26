;+
; NAME:
;   sdss_spec_block
; PURPOSE:
;   take a set of sdss spectra and output as a big block
; CALLING SEQUENCE:
;   sdss_spec_block
; INPUTS:
;   plate - [N] list of plates
;   fiberid - [N] list of fibers
;   mjd - [N] list of mjds
; OPTIONAL INPUTS:
;   avloglam - [nl] wavelength grid
;   vdisp - try to get everything to this velocity dispersion (in km/s)
; OPTIONAL KEYWORDS:
;   deextinct - apply corrections for galactic extinction
; OUTPUTS:
;   block_flux - [nl, N] list of fluxes (as output by readspec)
;   block_ivar - [nl, N] list of inverse variances
;   block_lamba - [nl] output wavelengths (angstroms)
; DEPENDENCIES:
;   idlutils
;   idlspec2d
; BUGS:
;   Input and output wavelength grids MUST be logarithmically spaced.
;   ivar is unchanged when vdisp is used
;   revision history wrong
; COMMENTS:
;   Units are in 10^{-17} ergs/cm^2/s/A 
;   Bolometric flux is kept constant in deredshifting step
;    (ie. lambda -> lambda/(1+z), flux -> flux*(1+z), ivar ->
;    ivar/(1+z)^2)
; REVISION HISTORY:
;   2001-11-17  written - Hogg (NYU)
;-
pro sdss_spec_block, plate, fiberid, mjd, block_flux=out_block_flux, $
                     block_ivar=out_block_ivar, $
                     block_lambda=out_block_lambda, $
                     avloglam=avloglam, deextinct=deextinct, vdisp=vdisp, $
                     minerror=minerror

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
        
; shift wavelength scale to rest-frame
        restloglam= loglam[*,i]-alog10(1.d + zans[i].z)

; de-extinct if desired
        if(keyword_set(deextinct)) then begin
            wave=10.^(loglam[*,i])
            alam=ext_ccm(wave)
            glactc, zans[i].plug_ra, zans[i].plug_dec, 2000., gl, gb, 1, /deg
            ebv=dust_getval(gl, gb)
            extvoebv=3.1
            ext=ebv*alam*extvoebv
            flux[*,i]=flux[*,i]*10.^(0.4*ext)
        endif

; slam onto new grid (MUST BE LOGARITHMIC SPACING)
        lum=fltarr(n_elements(restloglam))
        combine1fiber, restloglam,flux[*,i]*(1.+zans[i].z), $
          invvar[*,i]/(1.+zans[i].z)^2, newloglam=avloglam, newflux=tmp1, $
          newivar=tmp2, maxiter=0
        lum= double(tmp1)
        lum_ivar= double(tmp2)

; smooth if desired
        if(keyword_set(vdisp)) then begin
            adddisp2=vdisp^2-zans[i].vdisp^2
            if(adddisp2 gt 0.) then begin
;;                 sigma=sqrt(adddisp2)/(2.99792e+5*alog(10.))
                lum=k_smooth(avloglam,lum,sqrt(adddisp2))
                lum_ivar_new=k_smooth(avloglam,lum_ivar,sqrt(adddisp2))
                iz=where(lum_ivar eq 0, nz)
                if(nz gt 0) then $
                  lum_ivar_new[iz]=0.
                lum_ivar=lum_ivar_new
            endif
        endif

;; apply minimum errors
        if(keyword_set(minerror)) then begin
            inz=where(lum_ivar gt 0., nnz)
            if(nnz gt 0) then begin
                fracerr=1./(sqrt(lum_ivar[inz])*lum[inz])
                fracerr=sqrt(fracerr^2+minerror^2)
                lum_ivar[inz]=1./(fracerr*lum[inz])^2
            endif
        endif

; increment vectors
        nl=n_elements(lum)
        if(n_elements(block_flux) ne nl*nspec OR $
           n_elements(block_ivar) ne nl*nspec OR $
           n_elements(block_lambda) ne nl) then begin
            block_flux=fltarr(nl,nspec)
            block_ivar=fltarr(nl,nspec)
            block_lambda=fltarr(nl)
        endif 
        block_flux[*,plate_indx[i]]=lum
        block_ivar[*,plate_indx[i]]=lum_ivar
        block_lambda[*]=10.^(avloglam)
    endfor
endfor

out_block_flux=block_flux
out_block_ivar=block_ivar
out_block_lambda=block_lambda
    
return
end
