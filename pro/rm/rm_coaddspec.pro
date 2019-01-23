;+
; NAME:
;   rm_coadd1spec
; PURPOSE:
;   Coadd one object
;
; INPUTS:
;   inloglam       - Wavelengths in log10-Angstroms [NPIX,NSPEC]
;   objflux        - Flux [NPIX,NSPEC]
;
; REQUIRED KEYWORDS:
;   newloglam      - Wavelengths for output evaluation, also in log10-Angstroms
;                    [NNEWPIX]
;
; OPTIONAL INPUTS:
;   objivar        - Inverse variance [NPIX,NSPEC]
;   finalmask      - Pixel mask [NPIX,NSPEC]
;   indisp         - Dispersion values [NPIX,NSPEC]
;   skyflux        - Sky flux vectors [NPIX,NSPEC]
;   binsz          - Bin separation for INLOGLAM; if not set, then default
;                    to INLOGLAM[1]-INLOGLAM[0].
;   nord           - Order of spline fit; default to 3.
;   bkptbin        - Break point binning; default to 1.2 * BINSZ.
;   maxsep         - Maximum separation between input wavelengths.  The spline
;                    fit is split into pieces, with the breaks wherever this
;                    spacing is exceeded.  Default to 2.0 * BINSZ.
;   _EXTRA         - Keywords for DJS_REJECT().
;   verbose        - If set, then output messages about bad break points and
;                    masked data points.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   finalmask      - Modified from its input by setting the COMBINEREJ bit
;                    for deviant pixels in individual spectra that have
;                    been rejected.
;   objivar        - Modified from itts input by setting to zero wherever
;                    the COMBINEREJ bit has been set in FINALMASK.
;   newflux        - Resampled flux [NNEWPIX].
;   newivar        - Resampled inverse variance [NNEWPIX].
;   andmask        - Resampled mask. For each mask bit, set that bit only if
;                    every input spectrum at this wavelength has that bit set
;                    (e.g., this is a logical AND) [NNEWPIX].
;   ormask         - Resampled mask. For each mask bit, set that bit if any
;                    of the input spectra at this wavelength has that bit set
;                    (e.g., this is a logical OR) [NNEWPIX].
;   newdisp        - Resampled dispersion values [NNEWPIX].
;   newsky         - Resampled sky flux [NNEWPIX].
;
;------------------------------
pro rm_coadd1spec, inloglam, objflux, objivar, $
    inandmask=inandmask, inormask=inormask, indisp=indisp, skyflux=skyflux, $
    newloglam=newloglam, newflux=newflux, newivar=newivar, $
    andmask=andmask, ormask=ormask, newdisp=newdisp, newsky=newsky, $
    nord=nord, binsz=binsz

   if not keyword_set(binsz) then binsz=1d-4

   rm_combine1fiber, inloglam, objflux, objivar, $
          inandmask=inandmask, inormask=inormask, $
          indisp=indisp, skyflux=skyflux, $
          newloglam=newloglam, newflux=bestflux, newivar=bestivar, $
          andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
          newsky=bestsky, $
          nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep
          ; do not reject outliers during bspline-fit
          ; since individual epochs could vary in flux significantly
          ; maxiter=50, upper=3.0, lower=3.0, maxrej=1

   newflux = bestflux & newivar = bestivar
   andmask = bestandmask & ormask = bestormask
   newdisp = bestdispersion & newsky = bestsky

end
;
;-------------------------------
;+
; NAME:
;   rm_coaddspec
;
; PURPOSE:
;   Coadd all epochs for the objects on the RM plate, and output a single spPlate file
;
; OPTIONAL INPUTS:
;   specdir_prefix   - Path prefix of the individual spPlate* files, e.g., 'recalib/'
;                      To use the pipeline reduction, set specdir_prefix=''
;                      To use the WH sky-sub improved reduction, set specdir_prefix='wh_skysub/'
;                        and set the outfile name properly
;   mjd_coadd        - If set, only coadd those epochs specified by mjd
;   id_coadd         - Indices of coadded object in fibermap; default is all 1000 objects	
;
; REVISION HISTORY:
;   08-May-2014  Changed HUD0,1,4,6 outputs from double to float to match standard SDSS format
;------------------------------

pro rm_coaddspec, specdir_prefix=specdir_prefix, mjd_coadd=mjd_coadd, id_coadd=id_coadd, $
    outfile=outfile

   forward_function rm_combine1fiber

   topdir = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') + '/'


   ; Default is to use the new xyfit reduction       
   If not keyword_set(specdir_prefix) then specdir_prefix = 'wh_skysub/' ;'recalib/'

   ; Read in the master file
   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   ; If specified, only coadd these objects
   if n_elements(id_coadd) gt 0 then fibermap = fibermap[id_coadd]
   
   platearr = fibermap.plate
   fiberarr = fibermap.fiberid
   mjdarr = fibermap.mjd
   ind = where(platearr[*,0] gt 0, nnn)
   platearr = platearr[ind,0]
   fiberarr = fiberarr[ind,*]
   mjdarr = mjdarr[ind,0]
   ; If specified, only coadd these epochs
   if n_elements(mjd_coadd) gt 0 then begin
      flag = lonarr(nnn)
      for j=0L,n_elements(mjd_coadd) - 1L do begin
         indd = where(mjdarr eq mjd_coadd[j])
         if indd[0] ne -1 then flag[indd] = 1L
      endfor
      ind_mjd = where(flag eq 1)
      platearr = platearr[ind_mjd]
      fiberarr = fiberarr[ind_mjd,*]
      mjdarr = mjdarr[ind_mjd]
   endif
   nepoch = n_elements(platearr)
   nfiber = n_elements(fibermap)
   
   ; Set up the common wavelength grid of the final coadd
   wavemin=3.5523 & wavemax=4.0171
   binsz = 1.d-4
   spotmin = 0L
   spotmax = round((wavemax - wavemin)/binsz)
   nfinalpix = spotmax - spotmin + 1L
   newloglam = dindgen(nfinalpix) * binsz + wavemin
   
   wavearr = dblarr(nfinalpix,nfiber,nepoch)
   fluxarr = dblarr(nfinalpix,nfiber,nepoch)
   ivararr = dblarr(nfinalpix,nfiber,nepoch)
   skyarr = dblarr(nfinalpix,nfiber,nepoch)
   disparr = dblarr(nfinalpix,nfiber,nepoch)
   andmaskarr = lonarr(nfinalpix,nfiber,nepoch)
   ormaskarr = lonarr(nfinalpix,nfiber,nepoch)
   
   ; These are the final spPlate images [NFINALPIX,NFIBER]
   finalflux = dblarr(nfinalpix,nfiber)
   finalivar = dblarr(nfinalpix,nfiber)
   finalsky = dblarr(nfinalpix,nfiber)
   finaldisp = dblarr(nfinalpix,nfiber)
   finalandmask = lonarr(nfinalpix,nfiber)
   finalormask = lonarr(nfinalpix,nfiber)
   
   ; Get the plugmap structure from the first epoch
   rm_readspec,platearr[0],fiberarr[0,*],mjd=mjdarr[0],plugmap=plugmap1,/silent
   
   ; Populate each epoch
   for i=0L, nepoch - 1L do begin
   
      path = topdir + string(platearr[i],format='(i4.4)') + '/' + specdir_prefix
      rm_readspec,platearr[i],fiberarr[i,*],mjd=mjdarr[i],wave=wave1,flux=flux1, $
       invvar=ivar1,sky=sky1,disp=disp1,andmask=andmask1,ormask=ormask1,path=path,$
       /silent
      npix = (size(wave1))[1]
      if npix gt nfinalpix then npix = nfinalpix
      
      ; Note the wavelength array could differ slightly for different epochs
      wavearr[0:npix-1,*,i] = alog10(wave1)
      ; Set the extra wave array gracefully in increasing order
      if npix lt nfinalpix then begin
         for jj=0L,nfiber-1 do wavearr[npix:nfinalpix-1,jj,i] $
            = wavearr[npix-1,jj,i] + binsz*(1.+findgen(nfinalpix-npix))
      endif 

      fluxarr[0:npix-1,*,i] = flux1
      ivararr[0:npix-1,*,i] = ivar1
      skyarr[0:npix-1,*,i] = sky1
      disparr[0:npix-1,*,i] = disp1
      andmaskarr[0:npix-1,*,i]=andmask1
      ormaskarr[0:npix-1,*,i]=ormask1
   endfor
   
   ; Combine each object
   for i=0L, nfiber - 1L do begin
   
      inandmask = reform(andmaskarr[*,i,*])
      inormask = reform(ormaskarr[*,i,*])
      inloglam = reform(wavearr[*,i,*])
      objflux = reform(fluxarr[*,i,*])
      objivar = reform(ivararr[*,i,*])
      indisp = reform(disparr[*,i,*])
      skyflux = reform(skyarr[*,i,*])
      
      ;--------------------------------------------------------------
      ; Set the pixel ivar=0 if the inormask/inandmask is really bad
      ; so that these pixels won't be coadded
      ; What bits to use ???
      ; badmask=1 for bad pixel
      badmask = make_array(size=size(objivar),/long)
      ;badmask = badmask OR ((inormask AND pixelmask_bits('BADSKYCHI')) NE 0)
      objivar = objivar * (1 - badmask)
      ;--------------------------------------------------------------
   
      rm_coadd1spec, inloglam, objflux, objivar, $
       inandmask=inandmask, inormask=inormask, indisp=indisp, skyflux=skyflux, $
       newloglam=newloglam, newflux=newflux, newivar=newivar, $
       andmask=andmask, ormask=ormask, newdisp=newdisp, newsky=newsky, $
       nord=nord, binsz=binsz
   
      finalflux[*,i] = newflux
      finalivar[*,i] = newivar
      finalsky[*,i] = newsky
      finaldisp[*,i] = newdisp
      finalandmask[*,i] = andmask
      finalormask[*,i] = ormask
      
      splog, 'Finished coadding object: ', i+1
   
      ; message, 'stop'
 
   endfor
   
   ; Output the final spPlate file
   spfile1 = topdir+string(platearr[0],format='(i4.4)') + '/' + specdir_prefix $
      + 'spPlate-' + string(platearr[0],format='(i4.4)')+ '-' $
      + string(mjdarr[0],format='(i5.5)') + '.fits'
   tmp = mrdfits(spfile1,0,bighdr,/silent)
   tmp = mrdfits(spfile1,1,hdrfloat,/silent)
   tmp = mrdfits(spfile1,2,hdrlong,/silent)
   tmp = mrdfits(spfile1,3,hdrlong,/silent)
   tmp = mrdfits(spfile1,4,hdrfloat,/silent)
   tmp = mrdfits(spfile1,5,hdrplug,/silent)
   tmp = mrdfits(spfile1,6,hdrsky,/silent)
   totexptime = 0.
   for i=0L,nepoch-1 do begin
      spfile1 = topdir+string(platearr[i],format='(i4.4)') + '/' + specdir_prefix $
        + 'spPlate-' + string(platearr[i],format='(i4.4)')+ '-' $
        + string(mjdarr[i],format='(i5.5)') + '.fits'
      tmp = mrdfits(spfile1,0,tmphdr,/silent)
      thisexptime = fxpar(tmphdr,'EXPTIME')
      totexptime = totexptime + thisexptime
   endfor

   maxmjd = max(mjdarr)
   If not keyword_set(outfile) then $
    outfile = topdir + '0000/wh_skysub/spPlate-0000-'+string(maxmjd,format='(i5.5)')+'.fits'
   
   ; HDU #0 is flux
   sxaddpar, bighdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
   sxaddpar, bighdr, 'NAXIS1', nfinalpix
   sxaddpar, bighdr, 'NAXIS2', nfiber
   sxaddpar, bighdr, 'MJD', maxmjd   ; 99999
   mjdlist = mjdarr
   mjdlist = mjdlist[uniq(mjdlist, sort(mjdlist))]
   mjdlist = strtrim(strcompress(string(mjdlist,format='(99a)')),2)
   sxaddpar, bighdr, 'MJDLIST', mjdlist, after='MJD'
   sxaddpar, bighdr, 'COEFF0', wavemin, $
    ' Central wavelength (log10) of first pixel'
   sxaddpar, bighdr, 'COEFF1', binsz, ' Log10 dispersion per pixel'
   sxaddpar, bighdr, 'NEXP', nepoch, $
    ' Number of coadded epochs', before='EXPTIME'
   sxaddpar, bighdr, 'EXPTIME', totexptime, $
    ' Minimum of total exposure times for all cameras'
   spawn, 'uname -n', uname
   sxaddpar, bighdr, 'UNAME', uname[0]
   mwrfits, float(finalflux), outfile, bighdr, /create

   ; HDU #1 is inverse variance
   sxaddpar, hdrfloat, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
   sxaddpar, hdrfloat, 'EXTNAME', 'IVAR', ' Inverse variance'
   mwrfits, float(finalivar), outfile, hdrfloat

   ; HDU #2 is AND-pixelmask
   sxaddpar, hdrlong, 'EXTNAME', 'ANDMASK', ' AND Mask'
   mwrfits, finalandmask, outfile, hdrlong

   ; HDU #3 is OR-pixelmask
   sxaddpar, hdrlong, 'EXTNAME', 'ORMASK', ' OR Mask'
   mwrfits, finalormask, outfile, hdrlong

   ; HDU #4 is dispersion map
   sxaddpar, hdrfloat, 'BUNIT', 'pixels'
   sxaddpar, hdrfloat, 'EXTNAME', 'WAVEDISP', ' Wavelength dispersion'
   mwrfits, float(finaldisp), outfile, hdrfloat

   ; HDU #5 is plugmap
   ; Use the plugmap from the first epoch
   finalplugmap = plugmap1
   ; Reassign fiberid
   finalplugmap.fiberid = lindgen(nfiber)+1
   sxaddpar, hdrplug, 'EXTNAME', 'PLUGMAP', ' Plugmap structure'
   mwrfits, finalplugmap, outfile, hdrplug

   ; HDU #6 is the sky
   sxaddpar, hdrsky, 'EXTNAME', 'SKY', ' Subtracted sky flux'
   mwrfits, float(finalsky), outfile, hdrsky


end

