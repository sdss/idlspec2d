
;   objname    - These must be spFrame file names all from either spectro-1
;                or spectro-2, but not both!
;   adderr     - Additional error to add to the formal errors, as a
;                fraction of the flux; default to 0.03 (3 per cent).
;------------------------------------------------------------------------------
; Read the Kurucz models at the specified log-lambda and resolution
; ISELECT - If set, then only return these model numbers.
; Note that this function allocates a ridiculous amount of memory
; for caching the oversampled spectra at many dispersions in GRIDFLUX.

function spflux_read_kurucz, loglam, dispimg, iselect=iselect1, $
 kindx_return=kindx_return

   common com_spflux_kurucz, kfile, kflux, kindx, kloglam, nmodel, $
    gridsig, gridlam, gridflux

   ;----------
   ; Read the high-resolution Kurucz models

   if (NOT keyword_set(kfile)) then begin
      ; Read the file with the high-resolution Kurucz models as generated
      ; by the procedure KURUCZ_FITSFILE.
      splog, 'Reading Kurucz models'
      kfile = filepath('kurucz_stds_raw_v5.fits', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
      krawflux = mrdfits(kfile, 0, hdr, /silent)  ; flux
      kindx = mrdfits(kfile, 1, /silent)
      dims = size(krawflux, /dimens)
      nrawpix = dims[0]
      nmodel = dims[1]
      waves = dindgen(nrawpix) * sxpar(hdr,'CD1_1') + sxpar(hdr,'CRVAL1')

      ; The models will be sub-sampled relative to the SDSS pix scale of 1d-4:
      subsamp = 5
      kdlog = 1.d-4 / subsamp

      ; Change units to erg/cm^2/s/Ang by dividing by the wavelength
      for imodel=0L, nmodel-1 do $
       krawflux[*,imodel] = krawflux[*,imodel] / waves

      ; These models are sampled linearly in **air** wavelength.
      ; Re-map them to linear in (vacuum) log-wavelength.
      airtovac, waves ; Remap wavelengths from air to vacuum
      rawloglam = alog10(waves)
      splog, 'Remapping Kurucz models to log-wavelengths'
      minlog1 = min(rawloglam, max=maxlog1)
      kloglam = wavevector(minlog1, maxlog1, binsz=kdlog)
      npix = n_elements(kloglam)
      kflux = fltarr(npix,nmodel)
      for imodel=0L, nmodel-1 do $
       kflux[*,imodel] = rebin_spectrum(krawflux[*,imodel], rawloglam, kloglam)

      ; Convolve the Kurucz models with a boxcar response
      ; representing the size of the SDSS pixels.
      splog, 'Convolving Kurucz models with SDSS pixel size'
      if (subsamp GT 1) then begin
         if ((subsamp MOD 2) EQ 0) then $
          kern = ([0.5, fltarr(subsamp-1) + 1.0, 0.5]) / subsamp $
         else $
          kern = (fltarr(subsamp) + 1.0) / subsamp
         for imodel=0L, nmodel-1 do begin
            kflux[*,imodel] = convol(kflux[*,imodel], kern, /center)
         endfor
      endif

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
      gridsig = 0.70 + findgen(21) * 0.05 ; span range [0.7,1.7]
      nkpix = 7 * subsamp + 1
      nres = n_elements(gridsig)
      gridlam = kloglam ; Keep the same wavelength mapping
      gridflux = fltarr(npix, nres, nmodel)
      for ires=0L, nres-1 do begin
         print, ires, nres
         kern = exp(-0.5 * (findgen(nkpix*2+1) - nkpix)^2 $
          / (gridsig[ires]*subsamp)^2)
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
; 0 = not near lines, 1 = near lines
; HWIDTH = half width in log-wavelength for masking stellar lines

function spflux_masklines, loglam, hwidth=hwidth, stellar=stellar, $
 telluric=telluric

   if (NOT keyword_set(hwidth)) then $
    hwidth = 5.7e-4 ; Default is to mask +/- 5.7 pix = 400 km/sec

   mask = bytarr(size(loglam,/dimens))

   if (keyword_set(stellar)) then begin
      starwave = [ $
       3830.0 , $ ; ? (H-7 is at 3835 Ang)
       3889.0 , $ ; H-6
       3933.7 , $ ; Ca_k
       3968.5 , $ ; Ca_H (and H-5 at 3970. Ang)
       4101.7 , $ ; H-delta
       4300.  , $ ; G-band
       4305.  , $ ; G-band
       4310.  , $ ; more G-band
       4340.5 , $ ; H-gamma
       4861.3 , $ ; H-beta
       5893.0 , $ ; Mg
       6562.8 , $ ; H-alpha
       8500.8 , $
       8544.6 , $
       8665.0 , $
       8753.3 , $
       8866.1 , $
       9017.5 , $
       9232.0 ]
      airtovac, starwave

      for i=0L, n_elements(starwave)-1 do begin
         mask = mask OR (loglam GT alog10(starwave[i])-hwidth $
          AND loglam LT alog10(starwave[i])+hwidth)
      endfor
   endif

   if (keyword_set(telluric)) then begin
      tellwave1 = [6850., 7150., 7560., 8105., 8930.]
      tellwave2 = [6960., 7350., 7720., 8240., 9030.]
      for i=0L, n_elements(tellwave1)-1 do begin
         mask = mask OR (loglam GT alog10(tellwave1[i]) $
          AND loglam LT alog10(tellwave2[i]))
      endfor
   endif

   return, mask
end
;------------------------------------------------------------------------------
; Divide the spectrum by a median-filtered spectrum.
; The median-filtered version is computed ignoring stellar absorp. features.

function spflux_medianfilt, loglam, objflux, objivar, width=width, $
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

      ; For the median-filter, ignore points near stellar absorp. features,
      ; but keep points near telluric bands.
      qgood = 1 - spflux_masklines(loglam[*,ispec], /stellar)

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
   npix = dims[0]
   if (ndim EQ 1) then nspec = 1 $
    else nspec = dims[1]

   ;----------
   ; Median-filter the object fluxes

   medflux = spflux_medianfilt(loglam, objflux, objivar, $
    width=filtsz, /reflect, newivar=medivar)
   sqivar = sqrt(medivar)

   ;----------
   ; Mask out the telluric bands

   sqivar = sqivar * (1 - spflux_masklines(loglam, /telluric))

   ;----------
   ; Load the Kurucz models into memory

   junk = spflux_read_kurucz(kindx=kindx)
   nmodel = n_elements(kindx)

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

   ;----------
   ; Return the best-fit model

   minchi2 = min(chivec, ibest)
   dof = total(sqivar NE 0)
   splog, 'Best-fit total chi2/DOF = ', minchi2/(dof>1)
   bestflux = modflux[*,*,ibest]

   ;----------
   ; Compute the chi^2 just around the stellar absorp. lines
   ; for the best-fit model star

   linesqivar = sqivar * spflux_masklines(loglam, hwidth=12e-4, /stellar)
   linechi2 = 0.
   for ispec=0L, nspec-1 do begin
      thismodel = spflux_medianfilt(loglam, modflux[*,ispec,ibest], $
       width=filtsz, /reflect)
      linechi2 = linechi2 + computechi2(medflux[*,ispec], $
       linesqivar[*,ispec], thismodel)
   endfor
   linedof = total(linesqivar NE 0)
   splog, 'Best-fit line chi2/DOF = ', linechi2/(linedof>1)

   kindx1 = create_struct(kindx[ibest], $
    'IMODEL', ibest, $
    'Z', zpeak, $
    'CHI2', minchi2, $
    'DOF', dof, $
    'LINECHI2', linechi2, $
    'LINEDOF', linedof)

   ;----------
   ; Plot the filtered object spectrum, overplotting the best-fit Kurucz model

   ; Select the observation to plot that has the highest S/N,
   ; and one that goes blueward of 4000 Ang.
   snvec = total(objflux * sqrt(objivar), 1) $
    * (10.^loglam[0,*] LT 4000 OR 10.^loglam[npix-1,*] LT 4000)
   junk = max(snvec, iplot) ; Best blue exposure

   snvec = total(objflux * sqrt(objivar), 1) $
    * (10.^loglam[0,*] GT 8600 OR 10.^loglam[npix-1,*] GT 8600)
   junk = max(snvec, jplot) ; Best red exposure

   csize = 0.85
   djs_plot, [3840., 4120.], [0.0, 1.4], /xstyle, /ystyle, /nodata, $
    xtitle='Wavelength [Ang]', ytitle='Normalized Flux'
   if (iplot[0] NE -1) then begin
      djs_oplot, 10^loglam[*,iplot], medflux[*,iplot]
      djs_oplot, 10^loglam[*,iplot], medmodel[*,iplot], color='red'
   endif
   xyouts, 3860, 1.25, kindx1.model, charsize=csize
   djs_xyouts, 4000, 0.3, charsize=csize, $
    string(minchi2/(dof>1), format='("Total \chi^2/DOF=",f5.2)')
   djs_xyouts, 4000, 0.2, charsize=csize, $
    string(linechi2/(linedof>1), format='("Lines \chi^2/DOF=",f5.2)')
   djs_xyouts, 3860, 0.1, string(kindx1.feh, kindx1.teff, kindx1.g, $
    zpeak*cspeed, $
    format='("Fe/H=", f4.1, "  T_{eff}=", f6.0, "  g=", f3.1, "  cz=",f5.0)'), $
    charsize=csize

   djs_plot, [8440., 9160.], [0.0, 1.4], /xstyle, /ystyle, /nodata, $
    xtitle='Wavelength [Ang]', ytitle='Normalized Flux'
   if (jplot[0] NE -1) then begin
      djs_oplot, 10^loglam[*,jplot], medflux[*,jplot]
      djs_oplot, 10^loglam[*,jplot], medmodel[*,jplot], color='red'
   endif

   return, bestflux
end
;------------------------------------------------------------------------------
function spflux_goodfiber, pixmask
   qgood = ((pixmask AND pixelmask_bits('NOPLUG')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('BADTRACE')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('BADFLAT')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('BADARC')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('MANYBADCOLUMNS')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('NEARWHOPPER')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('MANYREJECTED')) EQ 0)
   return, qgood
end

;------------------------------------------------------------------------------
function spflux_bspline, loglam, mratio, mrativar, outmask=outmask, $
 everyn=everyn, airmass=airmass

   isort = sort(loglam)
   nord = 4 ; ???

   ; Choose the break points using the EVERYN option, but masking
   ; out more pixels near stellar features just when selecting them.
   mask1 = 1 - spflux_masklines(loglam, hwidth=12.e-4, /stellar)
   ii = where(mrativar[isort] GT 0 AND mask1[isort] EQ 1)
   bkpt = 0
   fullbkpt = bspline_bkpts(loglam[isort[ii]], everyn=everyn, $
    bkpt=bkpt, nord=nord)

   outmask1 = 0
   if (keyword_set(airmass)) then begin
      x2 = airmass[isort]
   endif
   sset = bspline_iterfit(loglam[isort], mratio[isort], $
    invvar=mrativar[isort], lower=3, upper=3, fullbkpt=fullbkpt, $
    maxrej=ceil(0.05*n_elements(indx)), outmask=outmask1, nord=nord, $
    x2=x2, npoly=2*keyword_set(airmass))
   if (arg_present(outmask)) then begin
      outmask = bytarr(size(loglam,/dimens))
      outmask[isort] = outmask1
   endif

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

   npoly = 3 ; ???
   flatarr = fltarr(npix, nobj)
   pres = fltarr(npoly, nobj)
   for iobj=0L, nobj-1 do begin
      ii = where(newivar[*,iobj] GT 0, ct)
      if (ct GT npoly+1) then begin ; At least NPOLY+1 pixels for a fit...
         thisloglam = newloglam[ii]
         thisratio = newratio[ii,iobj] / meanratio[ii]
         thisivar = newivar[ii,iobj] * meanratio[ii]^2

; The following fit should probably have some kind of rejection !!!???
         ; The following is a weighted fit...
         pres1 = poly_fit(thisloglam-3.5, thisratio, npoly-1, $
          measure_errors=1./sqrt(thisivar))

         ; The following would be an unweighted fit...
;         pres1 = poly_fit(thisloglam-3.5d0, thisratio, npoly-1)

         flatarr[*,iobj] = poly(loglam[*,iobj]-3.5d0, pres1)
         pres[*,iobj] = reform(pres1, npoly)
       endif else begin
         flatarr[*,iobj] = 1
         pres[*,iobj] = 0
         pres[0,iobj] = 1
       endelse
   endfor

   pres = reform(pres, [npoly, dims[1:ndim-1]])
   return, reform(flatarr, dims)
end

;------------------------------------------------------------------------------
pro spflux_plotcalib, mratiologlam, mratioflux, mrativar, $
 fitloglam, fitflux, fitflux2, logrange=logrange, plottitle=plottitle

   xrange = 10.^logrange
   ii = where(fitloglam GE logrange[0] AND fitloglam LE logrange[1])
   yrange = [min(fitflux[ii])>0.04*median(fitflux[ii])/1.1, $
    max(fitflux[ii])*1.1]
   nfinal = (size(mratioflux, /dimens))[2]

   djs_plot, xrange, yrange, /xstyle, /ystyle, /nodata, /ylog, $
    xtitle='Wavelength [Ang]', ytitle='Counts/(10^{-17}erg/cm^2/s/Ang', $
    title=plottitle
   for k=0, nfinal-1 do begin
      jj = where(mratiologlam[*,0,k] GE logrange[0] $
       AND mratiologlam[*,0,k] LE logrange[1] $
       AND mrativar[*,0,k] GT 0, ct)
      if (ct GT 1) then $
       djs_oplot, 10.^mratiologlam[jj,0,k], mratioflux[jj,0,k], psym=3
   endfor
   djs_oplot, 10.^fitloglam[ii], fitflux[ii], color='green'
   if (total(fitflux2) GT 0) then $
    djs_oplot, 10.^fitloglam[ii], fitflux2[ii], color='red'

   return
end

;------------------------------------------------------------------------------
pro spflux_v5, objname, adderr=adderr, combinedir=combinedir

   if (n_elements(adderr) EQ 0) then adderr = 0.03
   nfile = n_elements(objname)

   ;----------
   ; Get the list of spectrograph ID and camera names

   camname = strarr(nfile)
   expnum = lonarr(nfile)
   spectroid = lonarr(nfile)
   for ifile=0, nfile-1 do begin
      spframe_read, objname[ifile], hdr=hdr
      camname[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
      spectroid[ifile] = strmid(camname[ifile],1,1)
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
      splog, 'WARNING: No SPECTROPHOTO or REDDEN stars for flux calibration'
      return
   endif

   ;----------
   ; Replace the magnitudes for the F stars with the PSF fluxes
   ; from the calibObj files !!!???

   plateid = sxpar(hdr, 'PLATEID')
   thismjd = sxpar(hdr, 'MJD')
   tsobj = plug2tsobj(plateid, plugmap=plugmap)
   plugmap[iphoto].mag = 22.5 - 2.5 * alog10(tsobj[iphoto].psfflux)

   ;----------
   ; Read the raw F-star spectra

   npix = 2048
   loglam = fltarr(npix, nfile, nphoto)
   objflux = fltarr(npix, nfile, nphoto)
   objivar = fltarr(npix, nfile, nphoto)
   dispimg = fltarr(npix, nfile, nphoto)
   airmass = fltarr(npix, nfile, 320)
   for ifile=0L, nfile-1 do begin
      spframe_read, objname[ifile], iphoto, wset=wset1, loglam=loglam1, $
       objflux=objflux1, objivar=objivar1, dispimg=dispimg1, $
       mask=mask1, hdr=hdr1, adderr=adderr

      ; Compute the airmass for every pixel of every object
      ; (every pixel is the same, of course)
      tai = sxpar(hdr1, 'TAI')
      for j=0, 319 do $
       airmass[*,ifile,j] = tai2airmass(plugmap[j].ra, plugmap[j].dec, tai=tai)

      ; Make a map of the size of each pixel in delta-(log10-Angstroms).
      ; Re-normalize the flux to ADU/(dloglam).
      ; Re-normalize the dispersion from /(raw pixel) to /(new pixel).
      correct_dlam, objflux1, objivar1, wset1, dlam=dloglam
      correct_dlam, dispimg1, 0, wset1, dlam=dloglam, /inverse

      ; Mask pixels on bad fibers
      objivar1 = objivar1 * spflux_goodfiber(mask1)

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

   !p.multi = [0,2,3]
   modflux = 0 * objflux
   for ip=0L, nphoto-1 do begin
      ; Find the best-fit model -- evaluated for each exposure [NPIX,NEXP]
      thismodel = spflux_bestmodel(loglam[*,*,ip], objflux[*,*,ip], $
       objivar[*,*,ip], dispimg[*,*,ip], kindx=thisindx)

      ; Also evaluate this model over a big wavelength range [3006,10960] Ang.
      tmploglam = 3.4780d0 + lindgen(5620) * 1.d-4
      tmpdispimg = 0 * tmploglam + 1.0 ; arbitrarily select this resolution
      tmpflux = spflux_read_kurucz(tmploglam, tmpdispimg, $
       iselect=thisindx.imodel)

      ; The returned models are redshifted, but not fluxed or
      ; reddening-corrected.  Do that now...
      extcurve = ext_odonnell(10.^loglam[*,*,ip], 3.1)
      thismodel = thismodel * extcurve
      extcurve = ext_odonnell(10.^tmploglam, 3.1)
      tmpflux = tmpflux * 10.^(-extcurve * 3.1 * ebv[ip] / 2.5)

      ; Now integrate the apparent magnitude for this spectrum,
      ; The units of FTHRU are such that m = -2.5*alog10(FTHRU) + (48.6-2.5*17)
      wavevec = 10.d0^tmploglam
      flambda2fnu = wavevec^2 / 2.99792e18
      fthru = filter_thru(tmpflux * flambda2fnu, waveimg=wavevec, /toair)
      thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)

;      scalefac = 10.^((thisindx.mag[2] - plugmap[iphoto[ip]].mag[2])/2.5) ; ???
      scalefac = 10.^((thismag[2] - plugmap[iphoto[ip]].mag[2])/2.5)
      thismodel = thismodel * scalefac

      modflux[*,*,ip] = thismodel
      if (ip EQ 0) then kindx = replicate( create_struct( $
       'PLATE', 0L, $
       'MJD', 0L, $
       'FIBERID', 0L, $
       'QGOOD', 0, $
       thisindx, $
       'MODELFLUX', fltarr(npix)), $
       nphoto)
      copy_struct_inx, thisindx, kindx, index_to=ip
      kindx[ip].plate = plateid
      kindx[ip].mjd = thismjd
      kindx[ip].fiberid = iphoto[ip] + 1 + 320 * (spectroid[0] - 1)
   endfor
   !p.multi = 0

   ;----------
   ; Keep track of which F stars are good

   qfinal = bytarr(nphoto) + 1B

   ; Start with a rejection of any stars with a bad chi^2/DOF either
   ; in the full spectrum or in just the absorp. line regions.
   ; Do not reject more than half the stars.
   chi2limit = 2.0 ; ???
   chi2list = (kindx.chi2 / (kindx.dof>1)) $
    > (kindx.linechi2 / (kindx.linedof>1))
   chi2list = chi2list + 100 * (kindx.linedof LT 10) ; Bad if < 10 pixels
   while (max(chi2list) GT chi2limit AND total(qfinal) GT nphoto/2.) do begin
      junk = max(chi2list, iworst)
      splog, 'Rejecting std star in fiber = ', iphoto[iworst] $
       + 320 * strmatch(camname[0],'*2') + 1
      chi2list[iworst] = 0
      qfinal[iworst] = 0B
   endwhile

   ;----------
   ; Loop over each exposure, and compute the PCA fit to MRATIO
   ; using outlier-rejection.
   ; Iterate, rejecting entire stars if they are terrible fits.

   iblue = where(strmatch(camname,'b*'), nblue)
   ired = where(strmatch(camname,'r*'), nred)

   qdone = 0L
   while (qdone EQ 0) do begin
      ifinal = where(qfinal,nfinal) ; This is the list of the good stars

      ;----------
      ; The MRATIO vectors are the "raw" flux-calib vectors for each expos+CCD

      mratio = objflux / modflux
      mrativar = objivar * modflux^2

      ; Ignore regions around the stellar features
      mrativar = mrativar * (1 - spflux_masklines(loglam, /stellar))

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
      ; Do the B-spline fits for the blue CCDs.

      everyn = nblue * nfinal * 5
      sset_b = spflux_bspline(loglam[*,iblue,ifinal], $
       mratio[*,iblue,ifinal], mrativar[*,iblue,ifinal], $
       everyn=everyn, outmask=mask_b)

      ;----------
      ; Do the B-spline fits for the red CCDs.
      ; Fit a 2-dimension B-spline using the airmass as the 2nd dimension,
      ; but only if the airmass spans at least 0.10 and there are at
      ; least 3 good stars.

      everyn = nred * nfinal * 1.5
      if (max(airmass) - min(airmass) GT 0.10 AND nfinal GE 3) then begin ; ???
         ; Get an airmass value for every *pixel* being fit
         thisair = airmass[*,ired,iphoto[ifinal]]
      endif else begin
         thisair = 0
      endelse
      sset_r = spflux_bspline(loglam[*,ired,ifinal], $
       mratio[*,ired,ifinal], mrativar[*,ired,ifinal], $
       everyn=everyn, outmask=mask_r, airmass=thisair)

      ;----------
      ; Find which star has the most pixels rejected, and reject
      ; that star if it's bad enough

      fracgood = fltarr(nfinal)
      for k=0L, nfinal-1 do $
       fracgood[k] = 0.5 * (mean(mask_b[*,*,k]) + mean(mask_r[*,*,k]))
      minfrac = min(fracgood, iworst)
      if (minfrac LT 0.80) then begin
         if (nfinal LE nphoto/2.) then begin
            splog, 'WARNING: Already rejected ', nphoto-nfinal, ' of ', $
             nphoto, ' std stars'
            qdone = 1B
         endif else begin
            splog, 'Rejecting std star in fiber = ', iphoto[ifinal[iworst]] $
             + 320 * strmatch(camname[0],'*2') + 1
            qfinal[ifinal[iworst]] = 0B
         endelse
      endif else begin
         qdone = 1B ; No other stars to reject
      endelse
   endwhile

   kindx[ifinal].qgood = 1
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
      if (tag_exist(thisset,'NPOLY')) then $
       x2 = airmass[ifile,iphoto[ifinal]] $
      else $
       x2 = 0

      ; Evaluate the B-spline for the stars at their measured wavelengths
      ; in this exposure, then modulated by the mean FLATARR
      ; for the stars in this exposure.
      ; We re-fit the B-spline to exactly recover what we had before,
      ; just modulated by the lower-order polynomial FLATARR.

      logmin = min(thisloglam[where(mrativar GT 0)], max=logmax)
      tmploglam = wavevector(logmin, logmax)
      flatarr_mean = 0 * tmploglam
      for i=0L, nfinal-1 do $
       flatarr_mean = flatarr_mean $
        + poly(tmploglam-3.5d0, thispres[*,0,i]) / nfinal
      if (keyword_set(x2)) then begin
         x2_min = min(airmass[*,ifile,iphoto[ifinal]], max=x2_max)
         splog, 'Exposure ', objname[ifile], $
          ' spans airmass range ', x2_min, x2_max
         tmpflux1 = bspline_valu(tmploglam, thisset, x2=x2_min+0*tmploglam) $
          * flatarr_mean
         tmpflux2 = bspline_valu(tmploglam, thisset, x2=x2_max+0*tmploglam) $
          * flatarr_mean
         calibset = bspline_iterfit([tmploglam,tmploglam], $
          [tmpflux1,tmpflux2], oldset=thisset, maxiter=0, $
          x2=[0*tmploglam+x2_min,0*tmploglam+x2_max])
      endif else begin
         tmpflux1 = bspline_valu(tmploglam, thisset) * flatarr_mean
         calibset = bspline_iterfit(tmploglam, tmpflux1, oldset=thisset, $
          maxiter=0)
         tmpflux2 = 0
      endelse

      ;----------
      ; Make plots of the spectro-photometry data for this exposure only,
      ; overplotting the global fit to all exposures in red.

      ; The following info is just used for the plot title
      platestr = string(sxpar(hdr1,'PLATEID'), format='(i4.4)')
      mjdstr = string(sxpar(hdr1,'MJD'), format='(i5.5)')
      plottitle = 'PLATE=' + platestr + ' MJD=' + mjdstr $
       + ' Spectro-Photo Calibration for ' + camname[ifile]

      !p.multi = [0,1,2]
      logrange = logmax - logmin
      spflux_plotcalib, $
       thisloglam, thismratio, thismrativar, $
       tmploglam, tmpflux1/flatarr_mean, tmpflux2/flatarr_mean, $
       logrange=(logmin+[0,1]*logrange/2.), plottitle=plottitle
      spflux_plotcalib, $
       thisloglam, thismratio, thismrativar, $
       tmploglam, tmpflux1/flatarr_mean, tmpflux2/flatarr_mean, $
       logrange=(logmin+[1,2]*logrange/2.)
      !p.multi = 0

      ;----------
      ; Create header cards describing the fit range

      hdr = ['']
      sxaddpar, hdr, 'WAVEMIN', 10.^logmin
      sxaddpar, hdr, 'WAVEMAX', 10.^logmax

      ;----------
      ; Generate the pixel map of the flux-calibration for this exposure+CCD

      spframe_read, objname[ifile], loglam=loglam1
      if (tag_exist(calibset,'NPOLY')) then x2 = airmass[*,ifile,*] $
       else x2 = 0
      calibimg = bspline_valu(loglam1, calibset, x2=x2)

      ; Set to zero any pixels outside the known flux-calibration region
      ibad = where(loglam1 LT logmin OR loglam1 GT logmax, nbad)
      if (nbad GT 0) then calibimg[ibad] = 0

      ;----------
      ; Write the output file

      ; Put the Kurucz models for this exposure in the output structure
      kindx.modelflux = reform(modflux[*,ifile,*], npix, nphoto)

      calibfile = djs_filepath(string(camname[ifile], expnum[ifile], $
       format='("spFluxcalib-", a2, "-", i8.8, ".fits")'), $
       root_dir=combinedir)
      mwrfits, calibimg, calibfile, hdr, /create
      mwrfits, calibset, calibfile
      mwrfits, kindx, calibfile
   endfor

   return
end
;------------------------------------------------------------------------------
