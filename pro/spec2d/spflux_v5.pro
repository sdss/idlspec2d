
;   objname    - These must be spFrame file names all from either spectro-1
;                or spectro-2, but not both!
;   adderr     - Additional error to add to the formal errors, as a
;                fraction of the flux; default to 0.03 (3 per cent).

;------------------------------------------------------------------------------
; Read the Kurucz models at the specified log-lambda and resolution
; ISELECT - If set, then only return these model numbers
function spflux_read_kurucz, loglam, dispimg, iselect=iselect1, $
 kindx_return=kindx_return

   common com_spflux_kurucz, kfile, kflux, kindx, kloglam, nmodel, $
    gridsig, gridlam, gridflux

   ;----------
   ; Read the high-resolution Kurucz models

   if (NOT keyword_set(kfile)) then begin
      ; Read Christy's file with the high-resolution Kurucz models
; ARE THESE FILES ALREADY CONVOLVED WITH THE SDSS RESPONSE, IN WHICH
; CASE I'M DOING THIS TWICE!!!???
      splog, 'Reading Kurucz models'
      kfile = filepath('kurucz_stds_v5.fit', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
      kflux = mrdfits(kfile, 0, hdr, /silent)  ; flux
      kindx = mrdfits(kfile, 1, /silent)
      dims = size(kflux, /dimens)
      npix = dims[0]
      nmodel = dims[1]
      kdlam = sxpar(hdr,'CD1_1')
      kloglam = dindgen(npix) * kdlam + sxpar(hdr,'CRVAL1')

      ; Compute the spectro-photo fluxes for these spectra
      ; by just using these high-resolution redshift=0 spectra.
      ; Formally, I should re-compute these fluxes after redshifting
      ; and convolving with the dispersion, but I'm choosing to
      ; ignore those tiny effects.
;      wavevec = 10.d0^kloglam
;      flambda2fnu = wavevec^2 / 2.99792e18
;      fthru = filter_thru(kflux * rebin(flambda2fnu,npix,nmodel), $
;       waveimg=wavevec, /toair)
;      mag = - 2.5 * alog10( fthru * 10^((48.6 - 2.5*17.)/2.5) )
;      kindx.mag = transpose(mag)

      ; Construct a grid of these models at various possible resolutions
      splog, 'Convolving Kurucz models with dispersions'
      gridsig = 0.70 + findgen(50) * 0.02 ; span range [0.7,1.7]
      nkpix = 7
      nres = n_elements(gridsig)
      gridlam = kloglam ; Keep the same wavelength mapping
      gridflux = fltarr(npix, nres, nmodel)
      for ires=0L, nres-1 do begin
         kern = exp(-0.5 * (findgen(nkpix*2+1) - nkpix)^2 / (gridsig[ires])^2)
         kern = kern / total(kern)
         for imodel=0L, nmodel-1 do begin
            gridflux[*,ires,imodel] = convol(kflux[*,imodel], kern, /center)
         endfor
      endfor
   endif

   ;----------
   ; Return just with KINDX if LOGLAM,DISPIMG are not set

   if (n_elements(iselect1) GT 0) then iselect = iselect1 $
    else iselect = lindgen(nmodel)
   nselect = n_elements(iselect)
   if (arg_present(kindx_return)) then kindx_return = kindx[iselect]
   if (NOT keyword_set(loglam)) then return, 0

   ;----------
   ; Map the wavelength values and dispersion values onto pixel numbers
   ; on the grid of models.

   xpos = (loglam - gridlam[0]) / (gridlam[1] - gridlam[0])
   xpos = (xpos > 0) < (n_elements(gridlam)-1)
   ypos = (dispimg - gridsig[0]) / (gridsig[1] - gridsig[0])
   ypos = (ypos > 0) < (n_elements(gridsig)-1)

   ;----------
   ; Interpolate within these models

   ndim = size(loglam, /n_dimen)
   dims = size(loglam, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then nobj = 1 $
    else nobj = dims[1]
   modflux = fltarr(npix, nobj, nselect)
   for iobj=0L, nobj-1 do begin
      for j=0L, nselect-1 do begin
         imodel = iselect[j]
         modflux[*,iobj,j] = interpolate(gridflux[*,*,imodel], $
          xpos[*,iobj], ypos[*,iobj])
      endfor
   endfor

   return, modflux
end
;------------------------------------------------------------------------------
; Create a mask of 1's and 0's, where wavelenths that should not be used
; for fluxing (like near stellar features) are masked.

function spflux_masklines, loglam

   rejwave = [ $
    [3830.0 + [-8,8]] , $  ; ? (H-7 is at 3835 Ang)
    [3889.0 + [-8,8]] , $  ; H-6
    [3933.7 + [-8,8]] , $  ; Ca_k
    [3968.5 + [-8,8]] , $  ; Ca_H (and H-5 at 3970. Ang)
    [4101.7 + [-8,8]] , $  ; H-delta
    [[4290., 4320.]]  , $  ; G-band
    [4340.5 + [-8,8]] , $  ; H-gamma
    [4861.3 + [-10,10]] , $  ; H-beta
    [6562.8 + [-12,12]] ] ; H-alpha
   airtovac, rejwave
   nreject = n_elements(rejwave) / 2

   ; Mask =1 for good points, =0 for bad points
   mask = bytarr(size(loglam,/dimens)) + 1B
   for i=0L, nreject-1 do begin
      mask = mask AND (loglam LT alog10(rejwave[0,i]) $
       OR loglam GT alog10(rejwave[1,i]))
   endfor

   return, mask
end
;------------------------------------------------------------------------------
; Divide the spectrum by a median-filtered spectrum.
; The median-filtered version is computed ignoring stellar absorp. features.

function spflux_medianfilt, loglam, objflux, objivar, mask=mask, width=width, $
 newivar=newivar, _EXTRA=KeywordsForMedian

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then nspec = 1 $
    else nspec = dims[1]

; Should also put the telluric features here...???
   ; Specify where the lines are to mask out
   linewave = [3830.0, 3889.0, $
   ;         H-delta, Ca_k,   Ca_H,   G-band,     
              4101.7, 3933.7, 3968.5, 4300, 4310, $
   ;         H-gamma H-beta, Mgb, H-alpha
              4340.5, 4861.3, 5153.0, 6562.8]
   linewidth = 16.0 ; in Angstroms

   ;----------
   ; Loop over each spectrum

   medflux = 0 * objflux
   if (arg_present(objivar)) then newivar = 0 * objivar
   for ispec=0L, nspec-1 do begin

      ; Mask =1 for good points, =0 for bad points
      qgood = spflux_masklines(loglam[*,ispec])

      ; Median-filter, but skipping masked points
      igood = where(qgood, ngood)
      thisback = fltarr(dims[0])
      if (ngood GT 1) then begin
         thisback[igood] = djs_median(objflux[igood,ispec], width=width, $
          _EXTRA=KeywordsForMedian)
      endif
      thisback = djs_maskinterp(thisback, (qgood EQ 0), /const)

      ; Force the ends of the background to be the same as the spectrum,
      ; which will force the ratio of the two to be unity.
      hwidth = ceil((width-1)/2.)
      thisback[0:hwidth] = objflux[0:hwidth,ispec]
      thisback[npix-1-hwidth:npix-1] = objflux[npix-1-hwidth:npix-1,ispec]

      medflux[*,ispec] = objflux[*,ispec] / thisback
      if (arg_present(objivar)) then $
      newivar[*,ispec] = objivar[*,ispec] * thisback^2
   endfor

   return, medflux
end
;------------------------------------------------------------------------------
function spflux_bestmodel, loglam, objflux, objivar, dispimg, kindx=kindx1

   filtsz = 99 ; ???
   cspeed = 2.99792458e5

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   if (ndim EQ 1) then nspec = 1 $
    else nspec = dims[1]

   ;----------
   ; Median-filter the object fluxes

   medflux = spflux_medianfilt(loglam, objflux, objivar, mask=(objivar NE 0), $
    width=filtsz, /reflect, newivar=medivar)
   sqivar = sqrt(medivar)

   ;----------
   ; Load the Kurucz models into memory

   junk = spflux_read_kurucz(kindx=kindx)
   nmodel = n_elements(kindx)

; NEED TO MASK OUT 5577, etc !!!???
   ;----------
   ; Fit the redshift just by using a canonical model

   ifud = where(kindx.teff EQ 6000 AND kindx.g EQ 4 AND kindx.feh EQ -1.5)
   if (ifud[0] EQ -1) then $
    message, 'Could not find fiducial model!'
   nshift = 20
   logshift = (-nshift/2. + findgen(nshift)) * 1.d-4
   chivec = fltarr(nshift)
   for ishift=0L, nshift-1 do begin
      modflux = spflux_read_kurucz(loglam-logshift[ishift], $
       dispimg, iselect=ifud)
      ; Median-filter this model
      medmodel = spflux_medianfilt(loglam, modflux, $
       width=filtsz, /reflect)
      for ispec=0L, nspec-1 do begin
         chivec[ishift] = chivec[ishift] + computechi2(medflux[*,ispec], $
          sqivar[*,ispec], medmodel[*,ispec])
      endfor
   endfor
   zshift = (10.d^logshift - 1) ; Convert log-lambda shift to redshift
   zpeak = find_nminima(chivec, zshift, errcode=errcode)
   splog, 'Best-fit velocity for std star = ', zpeak * cspeed, ' km/s'
   if (errcode NE 0) then $
    splog, 'Warning: Error code ', errcode, ' fitting std star'

   ;----------
   ; Generate the Kurucz models at the specified wavelengths + dispersions,
   ; using the best-fit redshift

   modflux = spflux_read_kurucz(loglam-alog10(1.+zpeak), dispimg)

   ;----------
   ; Loop through each model, computing the best chi^2
   ; as the sum of the best-fit chi^2 to each of the several spectra
   ; for this same object.
   ; We do this after a median-filtering of both the spectra + the models.

   chiarr = fltarr(nmodel,nspec)
   chivec = fltarr(nmodel)
   for imodel=0L, nmodel-1 do begin
      ; Median-filter this model
      medmodel = spflux_medianfilt(loglam, modflux[*,*,imodel], $
       width=filtsz, /reflect)

      for ispec=0L, nspec-1 do begin
         chiarr[imodel,ispec] = computechi2(medflux[*,ispec], $
          sqivar[*,ispec], medmodel[*,ispec])
      endfor
      chivec[imodel] = total(chiarr[imodel,*])
   endfor

   ; Return the best-fit model
   minchi2 = min(chivec, ibest)
   dof = total(sqivar NE 0)
   bestflux = modflux[*,*,ibest]
   kindx1 = create_struct(kindx[ibest], 'IMODEL', ibest, 'Z', zpeak)

   ;----------
   ; Plot the filtered object spectrum, overplotting the best-fit Kurucz model
   ; Only plot the first spectrum -- this assumes it is the blue CCD ???

   djs_plot, [3840., 4120.], [0.0, 1.4], /xstyle, /ystyle, /nodata, $
    xtitle='Wavelength [Ang]', ytitle='Normalized Flux'
   djs_oplot, 10^loglam[*,0], medflux[*,0]
   djs_oplot, 10^loglam[*,0], medmodel[*,0], color='red'
   xyouts, 3860, 0.2, kindx1.model, charsize=1.5
   djs_xyouts, 4000, 0.2, $
    string(minchi2/dof, format='("\chi^2/DOF=",f5.2)'), charsize=1.5
   djs_xyouts, 3860, 0.1, string(kindx1.feh, kindx1.teff, kindx1.g, $
    zpeak*cspeed, $
    format='("Fe/H=", f4.1, "  T_{eff}=", f6.0, "  g=", f3.1, "  cz=",f5.0)'), $
    charsize=1.5

;set_plot,'x' ; ???
;imodel = ibest & ispec = 0
;medmodel = spflux_medianfilt(loglam, modflux[*,*,imodel], $
; width=filtsz, /reflect)
;foo=computechi2(medflux[*,ispec],sqivar[*,ispec],medmodel[*,ispec],yfit=yfit)
;stop

   return, bestflux
end
;------------------------------------------------------------------------------
pro spframe_read, filename, indx, objflux=objflux, objivar=objivar, $
 mask=mask, wset=wset, loglam=loglam, dispset=dispset, dispimg=dispimg, $
 plugmap=plugmap, skyflux=skyflux, hdr=hdr, adderr=adderr

   qtrim = n_elements(indx) GT 0

   if (NOT keyword_set(filename)) then $
    message, 'Must specify FILENAME'

   if (arg_present(hdr)) then hdr = headfits(filename)
   if (arg_present(objflux) $
    OR (arg_present(objivar) AND keyword_set(adderr))) then begin
      objflux = mrdfits(filename, 0, /silent)
      if (qtrim) then objflux = objflux[*,indx]
   endif
   if (arg_present(objivar)) then begin
      objivar = mrdfits(filename, 1, /silent)
      if (qtrim) then objivar = objivar[*,indx]
      if (keyword_set(adderr)) then begin
         gmask = objivar NE 0 ; =1 for good points
         objivar = 1.0 / ( 1.0/(objivar + (1-gmask)) $
          + (adderr * (objflux>0))^2 ) * gmask
      endif
   endif
   if (arg_present(mask)) then begin
      mask = mrdfits(filename, 2, /silent)
      if (qtrim) then mask = mask[*,indx]
   endif
   if (arg_present(wset) OR arg_present(loglam)) then begin
      wset = mrdfits(filename, 3, /silent)
      if (qtrim) then wset = traceset_trim(wset, indx)
      if (arg_present(loglam)) then traceset2xy, wset, xtmp, loglam
      xtmp = 0
   endif
   if (arg_present(dispset) OR arg_present(dispimg)) then begin
      dispset = mrdfits(filename, 4, /silent)
      if (qtrim) then dispset = traceset_trim(dispset, indx)
      if (arg_present(dispimg)) then traceset2xy, dispset, xtmp, dispimg
      xtmp = 0
   endif
   if (arg_present(plugmap)) then begin
      plugmap = mrdfits(filename, 5, /silent)
      if (qtrim) then plugmap = plugmap[indx]
   endif
   if (arg_present(skyflux)) then begin
      skyflux = mrdfits(filename, 6, /silent)
      if (qtrim) then skyflux = skyflux[*,indx]
   endif

   return
end

;------------------------------------------------------------------------------
function spflux_pca, loglam, mratio, mrativar, nkeep=nkeep, $
 newloglam=newloglam, acoeff=acoeff, outmask=outmask

   ; We'll want to reform all the arrays into NPIX * NSPECTRA,
   ; since they may be passed as NPIX * NSPECTRA1 * NSPECTRA2.
   dims = size(loglam, /dimens)
   npix = dims[0]
   nspec = n_elements(loglam) / npix

   ; Get the full wavelength mapping but on a uniform grid in log-wavelength
   minlog1 = min(loglam, max=maxlog1)
   newloglam = wavevector(minlog1, maxlog1)

   ; Reject pixels near stellar lines
   mask = spflux_masklines(loglam)

   pcaflux = pca_solve(reform(mratio,npix,nspec), $
    reform(mrativar*mask,npix,nspec), reform(loglam,npix,nspec), $
    nkeep=nkeep, nreturn=nkeep, newloglam=newloglam, $
    acoeff=acoeff, outmask=outmask)

   return, pcaflux
end

;------------------------------------------------------------------------------
pro spflux_v5, objname, adderr=adderr, combinedir=combinedir

   if (n_elements(adderr) EQ 0) then adderr = 0.03
   nfile = n_elements(objname)
   nkeep = 1 ; ???

   ;----------
   ; Get the list of spectrograph ID and camera names

   camname = strarr(nfile)
   expnum = lonarr(nfile)
   for ifile=0, nfile-1 do begin
      spframe_read, objname[ifile], hdr=hdr
      camname[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
      expnum[ifile] = sxpar(hdr, 'EXPOSURE')
   endfor

   ;----------
   ; Figure out which objects are F stars.
   ; Assume that the plug map is the same for all exposures.

   spframe_read, objname[0], plugmap=plugmap, hdr=hdr
   objtype = strtrim(plugmap.objtype,2)
   iphoto = where(objtype EQ 'SPECTROPHOTO_STD' OR objtype EQ 'REDDEN_STD', $
    nphoto)
   if (nphoto EQ 0) then begin
      splog, 'WARNING: No spectro-photo stars!'
      return
   endif

   ;----------
   ; Replace the magnitudes for the F stars with the PSF fluxes
   ; from the calibObj files !!!???

   plateid = sxpar(hdr, 'PLATEID')
   tsobj = plug2tsobj(plateid, plugmap=plugmap)
   plugmap[iphoto].mag = 22.5 - 2.5 * alog10(tsobj[iphoto].psfflux)

   ;----------
   ; Read the raw F-star spectra

   npix = 2048
   loglam = fltarr(npix, nfile, nphoto)
   objflux = fltarr(npix, nfile, nphoto)
   objivar = fltarr(npix, nfile, nphoto)
   dispimg = fltarr(npix, nfile, nphoto)
   for ifile=0L, nfile-1 do begin
      spframe_read, objname[ifile], iphoto, wset=wset1, loglam=loglam1, $
       objflux=objflux1, objivar=objivar1, dispimg=dispimg1, adderr=adderr

      ; Make a map of the size of each pixel in delta-(log10-Angstroms),
      ; and re-normalize the flux to ADU/(dloglam).
      ; Use the default value of 1.d-4 for DLOGLAM.
      correct_dlam, objflux1, objivar1, wset1, dlam=dloglam

      loglam[*,ifile,*] = loglam1
      objflux[*,ifile,*] = objflux1
      objivar[*,ifile,*] = objivar1
      dispimg[*,ifile,*] = dispimg1
   endfor

   ;----------
   ; Read the dust maps

   euler, plugmap[iphoto].ra, plugmap[iphoto].dec, ll, bb, 1
   ebv = dust_getval(ll, bb, /interp)

   ;----------
   ; For each star, find the best-fit model.

   modflux = 0 * objflux
   for ip=0L, nphoto-1 do begin
      thismodel = spflux_bestmodel(loglam[*,*,ip], objflux[*,*,ip], $
       objivar[*,*,ip], dispimg[*,*,ip], kindx=thisindx)

      ; The returned models are redshifted, but not fluxed or
      ; reddening-corrected.  Do that now...
      extcurve = ext_odonnell(10.^loglam[*,*,ip], 3.1)
      thismodel = thismodel * 10.^(-extcurve * 3.1 * ebv[ip] / 2.5)

      ; Now integrate the apparent magnitude for this spectrum,
      ; over the wavelength range [3006,10960] Ang.
      ; The units of FTHRU are such that m = -2.5*alog10(FTHRU) + (48.6-2.5*17)
      tmploglam = 3.4780d0 + lindgen(5620) * 1.d-4
      tmpdispimg = 0 * tmploglam + 1.0 ; arbitrarily select this resolution
      tmpflux = spflux_read_kurucz(tmploglam, tmpdispimg, $
       iselect=thisindx.imodel)
      wavevec = 10.d0^tmploglam
      flambda2fnu = wavevec^2 / 2.99792e18
      fthru = filter_thru(tmpflux * flambda2fnu, waveimg=wavevec, /toair)
      thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)

      thismodel = thismodel $
;       * 10.^((thisindx.mag[2] - plugmap[iphoto[ip]].mag[2])/2.5)
       * 10.^((thismag[2] - plugmap[iphoto[ip]].mag[2])/2.5)

      modflux[*,*,ip] = thismodel
      if (ip EQ 0) then kindx = replicate(thisindx, nphoto)
      kindx[ip] = thisindx

;set_plot,'x' ; ???
;ii = where(loglam[*,*,ip] GT alog10(3700.) $
; AND loglam[*,*,ip] LT alog10(4150),ct)
;if (ct GT 0) then begin
; splot,10^(loglam[*,*,ip])[ii], $
;  (objflux[*,*,ip])[ii]/(median(objflux[*,*,ip]))[ii], ps=3
; soplot,10^(loglam[*,*,ip])[ii], $
;  (modflux[*,*,ip])[ii]/(median(modflux[*,*,ip]))[ii],color='red', ps=3
;endif
   endfor

   ;----------
   ; The MRATIO vectors are the "raw" flux-calib vectors for each exposure+CCD

   mratio = objflux / modflux
   mrativar = objivar * modflux^2

   ;----------
   ; Keep track of which F stars are good

   qfinal = bytarr(nphoto) + 1B

   ;----------
   ; Loop over each exposure, and compute the PCA fit to MRATIO
   ; using outlier-rejection.
   ; Iterate, rejecting entire stars if they are terrible fits.

   iblue = where(strmatch(camname,'b*'), nblue)
   ired = where(strmatch(camname,'r*'), nred)
   qdone = 0L
   while (NOT qdone) do begin
      ifinal = where(qfinal,nfinal) ; This is the list of the good stars

      pca_b = spflux_pca(loglam[*,iblue,ifinal], $
       mratio[*,iblue,ifinal], mrativar[*,iblue,ifinal], nkeep=nkeep, $
       newloglam=loglam_b, acoeff=acoeff_b, outmask=mask_b)

      pca_r = spflux_pca(loglam[*,ired,ifinal], $
       mratio[*,ired,ifinal], mrativar[*,ired,ifinal], nkeep=nkeep, $
       newloglam=loglam_r, acoeff=acoeff_r, outmask=mask_r)

      ; Find which star has the most pixels rejected, and reject
      ; that star if it's bad enough
      fracgood_b = total(mask_b,nkeep) / (size(mask_b,/dimens))[0]
      fracgood_r = total(mask_r,nkeep) / (size(mask_r,/dimens))[0]
      fracgood = 0.5 * (fracgood_b + fracgood_r)
      minfrac = min(fracgood, iworst)
      if (minfrac LT 0.80) then begin
         if (nfinal LE nphoto/2.) then begin
            splog, 'WARNING: Already rejected ', nphoto-nfinal, ' of ', $
             nphoto, ' std stars'
            qdone = 1B
         endif else begin
            splog, 'Rejecting std star in fiber = ', iphoto[ifinal[iworst]] $
             + 320 * strmatch(camname[0],'*2') + 1
            ifinal[iworst] = 0B
         endelse
      endif else begin
         qdone = 1B ; No other stars to reject
      endelse
   endwhile
; SHOULD ALSO REJECT BASED UPON CHI^2 W.R.T. KURUCZ MODELS!!!???
   splog, 'Rejected ', nphoto-nfinal, ' of ', nphoto, ' std stars'

   ;----------
   ; Select break points for spline.
   ; These are the same values we have always used.
   
   bkptfile = filepath('blue.bkpts', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, bkptfile, allbkpts_b
   bkptfile = filepath('red.bkpts', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, bkptfile, allbkpts_r
   allbkpts_b = alog10(allbkpts_b) ; Convert to log-10 Angstroms
   allbkpts_r = alog10(allbkpts_r) ; Convert to log-10 Angstroms

   ;----------
   ; Construct the final (B-splined) flux-calibration vectors

   acoeff_b = reform(acoeff_b,nkeep,nblue,nfinal)
   acoeff_r = reform(acoeff_r,nkeep,nred,nfinal)
   for ifile=0, nfile-1 do begin
      ; Is this a blue CCD?
      ii = where(ifile EQ iblue, ct)
      if (ct EQ 1) then begin
         thisloglam = loglam_b
         thispca = pca_b
         thiscoeff = acoeff_b[*,ii,*]
         thisbkpts = allbkpts_b
      endif

      ; Is this a blue CCD?
      ii = where(ifile EQ ired, ct)
      if (ct EQ 1) then begin
         thisloglam = loglam_r
         thispca = pca_r
         thiscoeff = acoeff_r[*,ii,*]
         thisbkpts = allbkpts_r
      endif

      ; Construct a straight average of all the final F star
      ; reconstructed PCA fluxes.
      fluxcalib = 0 * thisloglam
      for j=0, nfinal-1 do $
       fluxcalib = fluxcalib + thispca # thiscoeff[*,0,j] / nfinal
; ???
;splot, 10^thisloglam, fluxcalib
;for j=0, nfinal-1 do $
; soplot, 10^thisloglam, thispca # thiscoeff[*,0,j], color=j+1

      ; Fit a B-spline to this flux-calibration vector, after
      ; masking the usual stellar absorp. features
      mask = spflux_masklines(thisloglam)
      indx = where(mask)
      logmin = min(thisloglam[indx], max=logmax)

      ; The following selects break points exactly like we always have,
      ; hard-wired to particular values.
;      ibk = where(thisbkpts GE logmin AND thisbkpts LE logmax, nbk)
;      if (nbk LT 4) then $
;       message, 'Error selecting break points'
;      bkpt = [logmin, thisbkpts[ibk[1:nbk-2]], logmax]

      ; The following generates a number of uniformly-spaced break points
      bkpt = 0
      fullbkpts = bspline_bkpts(thisloglam[indx], nord=4, nbkpts=30, $
       bkpt=bkpt, /silent)

      calibset = bspline_iterfit(thisloglam[indx], fluxcalib[indx], $
       nord=4, bkpt=bkpt, lower=3, upper=3, $
       maxrej=ceil(0.05*n_elements(indx)))
;foo = bspline_valu(thisloglam, calibset)
;splot, 10^thisloglam, fluxcalib
;soplot, 10^thisloglam, foo, color='green'

      ;----------
      ; Make plots

      xplot = thisloglam[sort(thisloglam)]
      yplot = bspline_valu(xplot, calibset)
      xrange = [10.^logmin, 10.^logmax]
      yrange = minmax(yplot)
      djs_plot, xrange, yrange, xrange=xrange, /xsytle, /nodata, $
       xtitle='Wavelength [Ang]', ytitle='Counts/(10^{-17}erg/cm^2/s/Ang'
      for j=0, nfinal-1 do $
       djs_oplot, 10.^thisloglam, thispca # thiscoeff[*,0,j], psym=3
      djs_oplot, 10.^xplot, yplot, color='red'

      ;----------
      ; Create header cards describing the fit range

      hdr = ['']
      sxaddpar, hdr, 'WAVEMIN', 10.^logmin
      sxaddpar, hdr, 'WAVEMAX', 10.^logmax

      ; Write the output file
      calibfile = djs_filepath(string(camname[ifile], expnum[ifile], $
       format='("spFluxcalib-", a2, "-", i8.8, ".fits")'), $
       root_dir=combinedir)
      mwrfits, 0, calibfile, hdr, /create
      mwrfits, calibset, calibfile
   endfor

   return
end
;------------------------------------------------------------------------------
