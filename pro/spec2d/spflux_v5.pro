
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
    [5893.0 + [-12,12]] , $  ; Mg
    [6562.8 + [-12,12]] , $ ; H-alpha
    [8500.8 + [-14,14]] , $ ; ?
    [8544.6 + [-14,14]] , $ ; ?
    [8665.0 + [-14,14]] , $ ; ?
    [8753.3 + [-14,14]] , $ ; ?
    [8866.1 + [-14,14]] , $ ; ?
    [9017.5 + [-14,14]] , $ ; ?
    [9232.0 + [-14,14]] ]   ; ?
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

; NEED TO MASK OUT 5577, TELLURIC BANDS, etc !!!???
   ;----------
   ; Fit the redshift just by using a canonical model

   ifud = where(kindx.teff EQ 6000 AND kindx.g EQ 4 AND kindx.feh EQ -1.5)
   if (ifud[0] EQ -1) then $
    message, 'Could not find fiducial model!'
   nshift = 20
   logshift = (-nshift/2. + findgen(nshift)) * 1.d-4
   chivec = fltarr(nshift)
   for ishift=0L, nshift-1 do begin
; VERIFY THAT THE DISPIMG BELOW IS IN EXACTLY THE SAME UNITS
; OF /10^-4 dloglam and not /pix !!!???
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
function spflux_bspline, loglam, mratio, mrativar, outmask=outmask, $
 _EXTRA=KeywordsForBkpts

   isort = sort(loglam)
   nord = 4

   ; The following generates break points spaced every Nth good value
   bkpt = 0
   fullbkpts = bspline_bkpts(loglam[isort], nord=4, _EXTRA=KeywordsForBkpts, $
    bkpt=bkpt, /silent)

   outmask = 0
   sset = bspline_iterfit(loglam[isort], mratio[isort], $
    invvar=mrativar[isort], nord=nord, bkpt=bkpt, lower=3, upper=3, $
    maxrej=ceil(0.05*n_elements(indx)), outmask=outmask)

   return, sset
end

;------------------------------------------------------------------------------
function spflux_mratio_flatten, loglam1, mratio1, mrativar1, pres=pres

   ;--------
   ; Re-form the input data arrays from multi-dimensional to N x M

   ndim = size(loglam1, /n_dimen)
   dims = size(loglam1, /dimens)
   npix = dims[0]
   nobj = n_elements(loglam1) / npix
   loglam = reform(loglam1, npix, nobj)
   mratio = reform(mratio1, npix, nobj)
   mrativar = reform(mrativar1, npix, nobj)

   ;--------
   ; Re-bin the spectra to the same spacing

   minlog1 = min(loglam, max=maxlog1)
   newloglam = wavevector(minlog1, maxlog1)
   nnewpix = n_elements(newloglam)

   newratio = fltarr(nnewpix, nobj)
   newivar = fltarr(nnewpix, nobj)

   for iobj=0L, nobj-1 do begin
      isort = sort(loglam[*,iobj])
      combine1fiber, loglam[isort,iobj], mratio[isort,iobj], $
       mrativar[isort,iobj], $
       newloglam=newloglam, newflux=newratio1, newivar=newivar1
      newratio[*,iobj] = newratio1
      newivar[*,iobj] = newivar1
   endfor

   ;--------
   ; Compute the straight weighted mean at each wavelength

   denom = total(newivar, 2) ; avoid divide-by-zeros
   meanratio = total(newratio * newivar, 2) / (denom + (denom EQ 0))

   ibadpix = where(meanratio LE 0, nbadpix)
   if (nbadpix GT 0) then newivar[ibadpix,*] = 0

   ;--------
   ; Now for each object, compute the polynomial fit of it relative to the mean

   npoly = 3
   flatarr = fltarr(npix, nobj)
   pres = fltarr(npoly, nobj)
   for iobj=0L, nobj-1 do begin
      ii = where(newivar[*,iobj] GT 0)
      thisloglam = newloglam[ii]
      thisratio = newratio[ii,iobj] / meanratio[ii]
      thisivar = newivar[ii,iobj] * meanratio[ii]^2
      pres1 = poly_fit(thisloglam, thisratio, npoly-1, $
       measure_errors=1./sqrt(thisivar))
      flatarr[*,iobj] = poly(loglam[*,iobj], pres1)
      pres[*,iobj] = reform(pres1, npoly)
   endfor

   pres = reform(pres, [npoly, dims[1:ndim-1]])
   return, reform(flatarr, dims)
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
       objflux=objflux1, objivar=objivar1, dispimg=dispimg1, $
       mask=mask1, adderr=adderr

      ; Make a map of the size of each pixel in delta-(log10-Angstroms),
      ; and re-normalize the flux to ADU/(dloglam).
      ; Use the default value of 1.d-4 for DLOGLAM.
      correct_dlam, objflux1, objivar1, wset1, dlam=dloglam

      loglam[*,ifile,*] = loglam1
      objflux[*,ifile,*] = objflux1
      objivar[*,ifile,*] = skymask(objivar1, mask1, mask1)
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
   ; Keep track of which F stars are good

   qfinal = bytarr(nphoto) + 1B

   iblue = where(strmatch(camname,'b*'), nblue)
   ired = where(strmatch(camname,'r*'), nred)

   ;----------
   ; Loop over each exposure, and compute the PCA fit to MRATIO
   ; using outlier-rejection.
   ; Iterate, rejecting entire stars if they are terrible fits.

   qdone = 0L
   while (NOT qdone) do begin
      ifinal = where(qfinal,nfinal) ; This is the list of the good stars

      ;----------
      ; The MRATIO vectors are the "raw" flux-calib vectors for each expos+CCD

      mratio = objflux / modflux
      mrativar = objivar * modflux^2

      ; Ignore regions around the stellar features
      mrativar = mrativar * spflux_masklines(loglam)

      ;----------
      ; For each camera (blue or red), divide-out a low-order polynomial from
      ; MRATIO each star to get them all to the same mean flux levels.
      ; This takes out gross throughput differences between exposures.
      ; Also, it will remove the ~5% large-scale spectrophotometry errors
      ; between invidual stars, both from spectrograph throughput variations
      ; and from slight mis-typing of the stars.

      flatarr_b = spflux_mratio_flatten(loglam[*,iblue,ifinal], $
       mratio[*,iblue,ifinal], mrativar[*,iblue,ifinal], pres=pres_b)
      mratio[*,iblue,ifinal] = mratio[*,iblue,ifinal] / flatarr_b
      mrativar[*,iblue,ifinal] = mrativar[*,iblue,ifinal] * flatarr_b^2

      flatarr_r = spflux_mratio_flatten(loglam[*,ired,ifinal], $
       mratio[*,ired,ifinal], mrativar[*,ired,ifinal], pres=pres_r)
      mratio[*,ired,ifinal] = mratio[*,ired,ifinal] / flatarr_r
      mrativar[*,ired,ifinal] = mrativar[*,ired,ifinal] * flatarr_r^2

      ;----------
      ; Do the B-spline fits

      everyn = nblue * nfinal * 3
      ii = where(mrativar[*,ired,ifinal] GT 0)
      sset_b = spflux_bspline((loglam[*,iblue,ifinal])[ii], $
       (mratio[*,iblue,ifinal])[ii], (mrativar[*,iblue,ifinal])[ii], $
       everyn=everyn, outmask=mask_b)
;set_plot,'x'
;foo=bspline_valu(loglam[*,iblue,ifinal],sset_b)
;splot,10^loglam[*,iblue,ifinal],mratio[*,iblue,ifinal],ps=3,yr=[0,20]
;ii=sort(loglam[*,iblue,ifinal])
;soplot,10^(loglam[*,iblue,ifinal])[ii],foo[ii]

      everyn = nred * nfinal * 1.5
      ii = where(mrativar[*,ired,ifinal] GT 0)
      sset_r = spflux_bspline((loglam[*,ired,ifinal])[ii], $
       (mratio[*,ired,ifinal])[ii], (mrativar[*,ired,ifinal])[ii], $
       everyn=everyn, outmask=mask_r)
;set_plot,'x'
;jj=where(mrativar[*,ired,ifinal] GT 0)
;splot,10^(loglam[*,ired,ifinal])[jj],(mratio[*,ired,ifinal])[jj],ps=3,yr=[0,15]
;foo=bspline_valu(loglam[*,ired,ifinal],sset_r)
;ii=sort(loglam[*,ired,ifinal])
;soplot,10^(loglam[*,ired,ifinal])[ii],foo[ii],color='red'

      ;----------
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
; SHOULD ALSO REJECT BASED UPON CHI^2 W.R.T. KURUCZ MODELS/ABSORP. LINES!!!???
   splog, 'Rejected ', nphoto-nfinal, ' of ', nphoto, ' std stars'

   ;----------
   ; Construct the final (B-splined) flux-calibration vectors

   for ifile=0, nfile-1 do begin
      ; Is this a blue CCD?
      ii = where(ifile EQ iblue, ct)
      if (ct EQ 1) then begin
         thisloglam = loglam[*,ifile,ifinal]
         thisset = sset_b
         thisflatarr = flatarr_b[*,ii,ifinal]
         thispres = pres_b[*,ii,ifinal]
      endif

      ; Is this a red CCD?
      ii = where(ifile EQ ired, ct)
      if (ct EQ 1) then begin
         thisloglam = loglam[*,ifile,ifinal]
         thisset = sset_r
         thisflatarr = flatarr_r[*,ii,ifinal]
         thispres = pres_r[*,ii,ifinal]
      endif
      thismratio = mratio[*,ifile,ifinal]
      thismrativar = mrativar[*,ifile,ifinal]

      ; Evaluate the B-spline for the stars at their measured wavelengths
      ; in this exposure, then modulated by the mean FLATARR
      ; for the stars in this exposure

      logmin = min(thisloglam[where(mrativar GT 0)], max=logmax)
      tmploglam = wavevector(logmin, logmax)
      flatarr_mean = 0 * tmploglam
      for i=0L, nfinal-1 do $
       flatarr_mean = flatarr_mean $
        + poly(tmploglam, thispres[*,0,i]) / nfinal
      tmpflux = bspline_valu(tmploglam, thisset) * flatarr_mean
      sset_tmp = bspline_iterfit(tmploglam, tmpflux, $
       nord=thisset.nord, bkpt=thisset.fullbkpt, maxiter=0)
;set_plot,'x'
;foo=bspline_valu(thisloglam, sset_tmp)
;splot,10^thisloglam,tmpflux,ps=3,yr=[0,20]
;soplot,10^(thisloglam)[ii],foo[ii]
if (min(tmpflux) LT 0 OR max(tmpflux) GT 1000) then stop ; ???

      ;----------
      ; Make plots

      xrange = [10.^logmin, 10.^logmax]
      djs_plot, 10.^tmploglam, tmpflux/flatarr_mean, $
       xrange=xrange, /xsytle, /nodata, $
       xtitle='Wavelength [Ang]', ytitle='Counts/(10^{-17}erg/cm^2/s/Ang'
      for j=0, nfinal-1 do $
       djs_oplot, 10.^thisloglam[*,0,j], thismratio[*,0,j], psym=3
      djs_oplot, 10.^tmploglam, tmpflux/flatarr_mean, color='red'

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
