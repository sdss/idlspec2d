;+
; NAME:
;   spflux_v5
;
; PURPOSE:
;   Compute flux-calibration vectors for each CCD+exposure from std stars
;
; CALLING SEQUENCE:
;   spflux_v5, objname, [ adderr=, combinedir=, minfracthresh= ]
;
; INPUTS:
;   objname    - File names (including path) for spFrame files, all from
;                either spectro-1 or spectro-2, but not both!
;
; OPTIONAL INPUTS:
;   adderr     - Additional error to add to the formal errors, as a
;                fraction of the flux; default to 0.03 (3 per cent)
;   combinedir - Directory for output files
;   minfracthresh - the threshold for the minimum number of good pixels
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   One output file is written for every input file OBJNAME,
;   prefixed with spFluxcalib.  It containts the following:
;     HDU #0: Calibration image
;     HDU #1: Trace-set used to construct the former
;     HDU #2: Structure with info on each standard star
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   airtovac
;   bspline_bkpts()
;   bspline_iterfit()
;   bspline_valu()
;   combine1fiber
;   computechi2()
;   correct_dlam
;   djs_filepath()
;   djs_maskinterp()
;   djs_median()
;   djs_oplot
;   djs_plot
;   djs_xyouts
;   euler
;   ext_odonnell
;   fileandpath()
;   filter_thru()
;   find_nminima()
;   get_tai
;   mrdfits()
;   mwrfits
;   pixelmask_bits()
;   rebin_spectrum()
;   spframe_read
;   skymask()
;   spflux_read_kurucz()
;   splog
;   tai2airmass()
;   wavevector()
;
; INTERNAL SUPPORT ROUTINES:
;   spflux_masklines()
;   spflux_medianfilt()
;   spflux_bestmodel()
;   spflux_goodfiber()
;   spflux_bspline()
;   spflux_mratio_flatten()
;   spflux_plotcalib
;   sxaddpar
;   sxpar()
;
; REVISION HISTORY:
;   05-Feb-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
; Create a mask of 1's and 0's, where wavelengths that should not be used
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
      czero2 = where(thisback eq 0., count2)
      if count2 gt 0 then thisback[czero2] = 1.
      medflux[*,ispec] = objflux[*,ispec] / thisback
      if (arg_present(objivar)) then $
      newivar[*,ispec] = objivar[*,ispec] * thisback^2
   endfor

   return, medflux
end
;------------------------------------------------------------------------------
function spflux_bestmodel, loglam, objflux, objivar, dispimg, kindx=kindx1, $
 plottitle=plottitle

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

   mlines = spflux_masklines(loglam, hwidth=12e-4, /stellar)
   linesqivar = sqivar * mlines
   linechi2 = 0.
   for ispec=0L, nspec-1 do begin
      thismodel = spflux_medianfilt(loglam, modflux[*,ispec,ibest], $
       width=filtsz, /reflect)
      linechi2 = linechi2 + computechi2(medflux[*,ispec], $
       linesqivar[*,ispec], thismodel)
   endfor
   linedof = total(linesqivar NE 0)
   splog, 'Best-fit line chi2/DOF = ', linechi2/(linedof>1)

   ;----------
   ; Compute the median S/N for all the spectra of this object,
   ; and for those data just near the absorp. lines

   sn_median = median(objflux * sqrt(objivar))
   indx = where(mlines, ct)
   if (ct GT 1) then $
    linesn_median = median(objflux[indx] * sqrt(objivar[indx])) $
   else $
    linesn_median = 0.
   splog, 'Full median S/N = ', sn_median
   splog, 'Line median S/N = ', linesn_median

   kindx1 = create_struct(kindx[ibest], $
    'IMODEL', ibest, $
    'Z', zpeak, $
    'SN_MEDIAN', sn_median, $
    'CHI2', minchi2, $
    'DOF', dof, $
    'LINESN_MEDIAN', linesn_median, $
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
    xtitle='Wavelength [Ang]', ytitle='Normalized Flux', $
    title=plottitle
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
   nord = 3

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
    x2=x2, npoly=2*keyword_set(airmass), requiren=(everyn-1)>1)
   if (max(sset.coeff) EQ 0) then $
    message, 'B-spline fit failed!!'
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
   ; (Avoid divide-by- zeros.)

   if (ndim EQ 1) then begin
      meanratio = (newratio * newivar) / (newivar + (newivar EQ 0))
   endif else begin
      denom = total(newivar, 2)
      meanratio = total(newratio * newivar, 2) / (denom + (denom EQ 0))
   endelse

   qbadpix = meanratio LE 0
   ibadpix = where(qbadpix, nbadpix)
   if (nbadpix GT 0) then newivar[ibadpix,*] = 0

   ;--------
   ; Actually take this "mean" and turn it into something more like
   ; a median, to protect us against standard stars that have bad
   ; magnitudes from the imaging.

; Comment-out ???
;   igoodpix = where(qbadpix EQ 0)
;   if (ndim EQ 1) then medratio = newratio $
;    else medratio = djs_median(newratio, 2)
;   rescale = median( medratio[igoodpix] / meanratio[igoodpix] )
;   if (rescale LE 0) then begin
;      splog, 'Warning: RESCALE = ', rescale
;   endif else begin
;      meanratio = rescale * meanratio
;      splog, 'Rescale factor median/mean = ', rescale
;   endelse

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

         ; This fit requires no rejection, because this function falls
         ; within an iteration loop that rejects points.

         ; The following is a weighted fit...
;         pres1 = poly_fit(thisloglam-3.5d0, thisratio, npoly-1, $
;          measure_errors=1./sqrt(thisivar))

         ; The following would be an unweighted fit...
;         pres1 = poly_fit(thisloglam-3.5d0, thisratio, npoly-1)

         ; The following is an unweighted fit but with outlier-rejection...
         poly_iter, thisloglam-3.5d0, thisratio, npoly-1, 3., coeff=pres1

         flatarr[*,iobj] = poly(loglam[*,iobj]-3.5d0, pres1)
         pres[*,iobj] = reform(pres1, npoly)
       endif else begin
         flatarr[*,iobj] = 1
         pres[*,iobj] = 0
         pres[0,iobj] = 1
       endelse
   endfor

   if (ndim GT 1) then $
    pres = reform(pres, [npoly, dims[1:ndim-1]])
   return, reform(flatarr, dims)
end

;------------------------------------------------------------------------------
pro spflux_plotcalib, mratiologlam, mratioflux, mrativar, $
 fitloglam, fitflux, fitflux2, logrange=logrange, plottitle=plottitle

   xrange = 10.^logrange
   ii = where(fitloglam GE logrange[0] AND fitloglam LE logrange[1])
   yrange = [0.9 * min(fitflux[ii]), 1.1 * max(fitflux[ii])]
   if (size(mratioflux, /n_dimen) EQ 1) then nfinal = 1 $
    else nfinal = (size(mratioflux, /dimens))[2]

   djs_plot, xrange, yrange, /xstyle, /ystyle, /nodata, $
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
pro spflux_v5, objname, adderr=adderr, combinedir=combinedir, $
 minfracthresh=minfracthresh
   if (not keyword_set(minfracthresh)) then minfracthresh=0.80
   if (n_elements(adderr) EQ 0) then adderr = 0.03
   nfile = n_elements(objname)

   ;----------
   ; Get the list of spectrograph ID and camera names

   plateid = lonarr(nfile)
   mjd = lonarr(nfile)
   camname = strarr(nfile)
   expnum = lonarr(nfile)
   spectroid = lonarr(nfile)
   npixarr = lonarr(nfile)
   for ifile=0, nfile-1 do begin
      spframe_read, objname[ifile], hdr=hdr
      plateid[ifile] = strtrim(sxpar(hdr, 'PLATEID'),2)
      mjd[ifile] = strtrim(sxpar(hdr, 'MJD'),2)
      camname[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
      spectroid[ifile] = strmid(camname[ifile],1,1)
      expnum[ifile] = sxpar(hdr, 'EXPOSURE')
      npixarr[ifile] = sxpar(hdr, 'NAXIS1')
   endfor
   maxmjd = max(mjd)

   ;----------
   ; Figure out which objects are F stars.
   ; Assume that the plug map is the same for all exposures.

   spframe_read, objname[0], plugmap=plugmap, hdr=hdr
   objtype = strtrim(plugmap.objtype,2)
   iphoto = where((objtype EQ 'SPECTROPHOTO_STD' OR objtype EQ 'REDDEN_STD') $
    AND plugmap.offsetid EQ 1, nphoto)
   if (nphoto EQ 0) then begin
      splog, 'WARNING: No SPECTROPHOTO or REDDEN stars for flux calibration'
      return
   endif

   ;----------
   ; Read the raw F-star spectra

   npix = max(npixarr)
   nfiber = n_elements(plugmap) ; Number of fibers in one spectrograph
   loglam = fltarr(npix, nfile, nphoto)
   objflux = fltarr(npix, nfile, nphoto)
   objivar = fltarr(npix, nfile, nphoto)
   dispimg = fltarr(npix, nfile, nphoto)
   airmass = fltarr(npix, nfile, nfiber)
   for ifile=0L, nfile-1 do begin
      spframe_read, objname[ifile], iphoto, wset=wset1, loglam=loglam1, $
       objflux=objflux1, objivar=objivar1, dispimg=dispimg1, $
       mask=mask1, hdr=hdr1, adderr=adderr

      ; Compute the airmass for every pixel of every object
      ; (every pixel is the same, of course)
      get_tai, hdr1, tai_beg, tai, tai_end
      for j=0, nfiber-1 do $
       airmass[0:npixarr[ifile]-1,ifile,j] = tai2airmass(plugmap[j].ra, plugmap[j].dec, tai=tai)

      ; Make a map of the size of each pixel in delta-(log10-Angstroms).
      ; Re-normalize the flux to ADU/(dloglam).
      ; Re-normalize the dispersion from /(raw pixel) to /(new pixel).
      correct_dlam, objflux1, objivar1, wset1, dlam=dloglam
      correct_dlam, dispimg1, 0, wset1, dlam=dloglam, /inverse

      ; Mask pixels on bad fibers
      objivar1 = objivar1 * spflux_goodfiber(mask1)

      loglam[0:npixarr[ifile]-1,ifile,*] = loglam1
      ;it wont do having a tail of zeros in the wavelength so add some dummy values
      if (npix GT npixarr[ifile]) then begin
        dllam=loglam1[npixarr[ifile]-1,*]-loglam1[npixarr[ifile]-2,*]
        for j=0, nphoto-1 do $
           loglam[npixarr[ifile]:*,ifile,j] = loglam1[npixarr[ifile]-1,j]+dllam[0,j]*(1+findgen(npix-npixarr[ifile]))
      endif
      objflux[0:npixarr[ifile]-1,ifile,*] = objflux1
      ;hopefully the inverse variance of 0 of non-filled objects will indicate the uselessness
      ; of the extra
      objivar[0:npixarr[ifile]-1,ifile,*] = skymask(objivar1, mask1, mask1)
      dispimg[0:npixarr[ifile]-1,ifile,*] = dispimg1
   endfor

   ;----------
   ; Keep track of which F stars are good

   qfinal = bytarr(nphoto) + 1B

   ;----------
   ; For each star, find the best-fit model.

   !p.multi = [0,2,3]
   modflux = 0 * objflux
   for ip=0L, nphoto-1 do begin
      thisfiber = iphoto[ip] + 1 + nfiber * (spectroid[0] - 1)
      splog, prelog='Fiber '+string(thisfiber,format='(I4)')

      plottitle = 'PLATE=' + string(plateid[0], format='(i4.4)') $
       + ' MJD=' + string(maxmjd, format='(i5.5)') $
       + ' Spectro-Photo Star' $
       + ' Fiber ' + strtrim(thisfiber,2)

      ; Find the best-fit model -- evaluated for each exposure [NPIX,NEXP]
      thismodel = spflux_bestmodel(loglam[*,*,ip], objflux[*,*,ip], $
       objivar[*,*,ip], dispimg[*,*,ip], kindx=thisindx, plottitle=plottitle)

      ; Also evaluate this model over a big wavelength range [3006,10960] Ang.
      tmploglam = 3.4780d0 + lindgen(5620) * 1.d-4
      tmpdispimg = 0 * tmploglam + 1.0 ; arbitrarily select this resolution
      tmpflux = spflux_read_kurucz(tmploglam, tmpdispimg, $
       iselect=thisindx.imodel)

      ; The returned models are redshifted, but not fluxed or
      ; reddened.  Do that now...  we compare data vs. model reddened.
      extcurve1 = ext_odonnell(10.^loglam[*,*,ip], 3.1)
      thismodel = thismodel $
       * 10.^(-extcurve1 * 3.1 * plugmap[iphoto[ip]].sfd_ebv / 2.5)
      extcurve2 = ext_odonnell(10.^tmploglam, 3.1)
      tmpflux = tmpflux $
       * 10.^(-extcurve2 * 3.1 * plugmap[iphoto[ip]].sfd_ebv / 2.5)

      ; Now integrate the apparent magnitude for this spectrum,
      ; The units of FTHRU are such that m = -2.5*alog10(FTHRU) + (48.6-2.5*17)
      ; Note that these computed magnitudes, THISMAG, should be equivalent
      ; to THISINDX.MAG in the case of no reddening.
      wavevec = 10.d0^tmploglam
      flambda2fnu = wavevec^2 / 2.99792e18
      fthru = filter_thru(tmpflux * flambda2fnu, waveimg=wavevec, /toair)
      thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)

      ; Compute SCALEFAC = (plugmap flux) / (uncalibrated flux)
      if (tag_exist(plugmap, 'CALIBFLUX')) then begin
         scalefac = plugmap[iphoto[ip]].calibflux[2] $
          / 10.^((22.5-thismag[2])/2.5)
         ; Reject this star if we don't know its flux.
         if (plugmap[iphoto[ip]].calibflux[2] LE 0) then begin
            splog, 'Warning: Rejecting std star in fiber = ', $
             iphoto[ip] + 1 + nfiber * (spectroid[0] - 1), $
             ' with unknown calibObj flux'
            qfinal[ip] = 0
         endif
      endif else begin
         splog, 'WARNING: No CALIBFLUX for zero-pointing the fluxes'
         scalefac = 10.^((thismag[2] - plugmap[iphoto[ip]].mag[2])/2.5)
      endelse
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
      kindx[ip].plate = plateid[0]
      kindx[ip].mjd = maxmjd
      kindx[ip].fiberid = thisfiber
      splog, prelog=''
   endfor
   !p.multi = 0

   ;----------
   ; Reject any stars where more than 20% of the pixels marked are bad
   ; in any observation.

   fracgood = fltarr(nphoto)
   for ip=0L, nphoto-1 do begin
      for i=0L, nfile-1 do begin
         markasbad = (qfinal[ip]) AND (mean(objivar[*,i,ip] GT 0) LT minfracthresh)
         if (markasbad) then begin
            splog, 'Warning: Rejecting std star in fiber = ', $
             iphoto[ip] + 1 + nfiber * (spectroid[0] - 1), $
             ' with too many IVAR=0 pixels'
            qfinal[ip] = 0B
         endif
      endfor
   endfor
   ifinal = where(qfinal,nfinal) ; This is the list of the good stars
   if (nfinal EQ 0) then begin
      splog, 'ABORT: No good fluxing stars!'
      return
   endif

   ;----------
   ; Reject any stars with a bad chi^2/DOF either
   ; in the full spectrum or in just the absorp. line regions.
   ; Do not reject more than half the stars.

   chi2limit = 2.0 ; ???
   chi2list = (kindx.chi2 / (kindx.dof>1)) $
    > (kindx.linechi2 / (kindx.linedof>1))
   chi2list = chi2list + 100 * (kindx.linedof LT 10) ; Bad if < 10 pixels
   while (max(chi2list) GT chi2limit AND total(qfinal) GT nphoto/2.) do begin
      chi2max = max(chi2list, iworst)
      splog, 'Rejecting std star in fiber = ', $
       iphoto[iworst] + 1 + nfiber * (spectroid[0] - 1), $
       ' with chi2=', chi2max
      chi2list[iworst] = 0
      qfinal[iworst] = 0B
   endwhile

   ;----------
   ; Reject any stars with a very low S/N in the absorp. line regions.
   ; Do not reject more than half the stars.

   snlimit = 2.0 ; ???
   ; Do not reject any stars that are already rejected above...
   snlist = kindx.linesn_median + (snlimit+1) * (qfinal EQ 0)
   while (min(snlist) LT snlimit AND total(qfinal) GT nphoto/2.) do begin
      snmin = min(snlist, iworst)
      splog, 'Rejecting std star in fiber = ', $
       iphoto[iworst] + 1 + nfiber * (spectroid[0] - 1), $
       ' with median line S/N=', snmin
      snlist[iworst] = snlimit + 1
      qfinal[iworst] = 0B
   endwhile

   if (total(qfinal) EQ 0) then begin
      splog, 'ABORT: No good spectro-photo stars!'
      return
   endif

   ;----------
   ; Loop over each exposure, and compute the PCA fit to MRATIO
   ; using outlier-rejection.
   ; Iterate, rejecting entire stars if they are terrible fits.

   iblue = where(strmatch(camname,'b*'), nblue)
   ired = where(strmatch(camname,'r*'), nred)

   qdone = 0L
   iiter = 0L
   while (qdone EQ 0) do begin
      iiter = iiter + 1
      splog, 'Iteration #', iiter

      ifinal = where(qfinal,nfinal) ; This is the list of the good stars

      ;----------
      ; The MRATIO vectors are the "raw" flux-calib vectors for each expos+CCD

      mmask = modflux GT 0
      mratio = objflux / (modflux*mmask + (1-mmask))
      mrativar = objivar * modflux^2
      flatarr = 0 * mratio

      ; Ignore regions around the stellar features
      mrativar = mrativar * (1 - spflux_masklines(loglam, /stellar))

      ;----------
      ; For each camera (blue or red), divide-out a low-order polynomial from
      ; MRATIO each star to get them all to the same mean flux levels.
      ; This takes out gross throughput differences between exposures.
      ; Also, it will remove the ~5% large-scale spectrophotometry errors
      ; between invidual stars, both from spectrograph throughput variations
      ; and from slight mis-typing of the stars.

      if (nblue GT 0) then $
       flatarr[*,iblue,ifinal] = spflux_mratio_flatten(loglam[*,iblue,ifinal], $
        mratio[*,iblue,ifinal], mrativar[*,iblue,ifinal], pres=pres_b)

      if (nred GT 0) then $
       flatarr[*,ired,ifinal] = spflux_mratio_flatten(loglam[*,ired,ifinal], $
        mratio[*,ired,ifinal], mrativar[*,ired,ifinal], pres=pres_r)

      mratio[*,*,ifinal] = mratio[*,*,ifinal] / flatarr[*,*,ifinal]
      mrativar[*,*,ifinal] = mrativar[*,*,ifinal] * flatarr[*,*,ifinal]^2

      ;----------
      ; Do the B-spline fits for the blue CCDs.

      if (nblue GT 0) then begin
         everyn = nblue * nfinal * 10
         sset_b = spflux_bspline(loglam[*,iblue,ifinal], $
          mratio[*,iblue,ifinal], mrativar[*,iblue,ifinal], $
          everyn=everyn, outmask=mask_b)
      endif

      ;----------
      ; Do the B-spline fits for the red CCDs.
      ; Fit a 2-dimension B-spline using the airmass as the 2nd dimension,
      ; but only if the airmass spans at least 0.10 and there are at
      ; least 3 good stars.

      if (nred GT 0) then begin
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
      endif

      ;----------
      ; Find which star has the most pixels rejected, and reject
      ; that star if it's bad enough

      fracgood = fltarr(nfinal)
      for k=0L, nfinal-1 do begin
         if (nblue GT 0) then fracgood[k] += 0.5 * mean(mask_b[*,*,k])
         if (nred GT 0) then fracgood[k] += 0.5 * mean(mask_r[*,*,k])
      endfor
      minfrac = min(fracgood, iworst)
      if (minfrac LT minfracthresh) then begin
         if (nfinal LE nphoto/2.) then begin
            splog, 'WARNING: Already rejected ', nphoto-nfinal, ' of ', $
             nphoto, ' std stars'
            qdone = 1B
         endif else begin
            splog, 'Rejecting std star in fiber = ', $
             iphoto[ifinal[iworst]] + 1 + nfiber * (spectroid[0] - 1), $
             ' with fracgood=', minfrac
            qfinal[ifinal[iworst]] = 0B
         endelse
      endif else begin
         qdone = 1B ; No other stars to reject
      endelse
   endwhile

   kindx[ifinal].qgood = 1
   splog, 'Rejected ', nphoto-nfinal, ' of ', nphoto, ' std stars'

   ;----------
   ; Plot fluxing vectors and their polynomial offsets for individual stars
   ; in individual exposures.

   mratfit = 0 * mratio
   if (nblue GT 0) then $
    mratfit[*,iblue,ifinal] = bspline_valu(loglam[*,iblue,ifinal], sset_b)
   if (nred GT 0) then $
    if (tag_exist(sset_r,'NPOLY')) then $
     mratfit[*,ired,ifinal] = bspline_valu(loglam[*,ired,ifinal], sset_r, $
      x2=airmass[*,ired,iphoto[ifinal]]) $
    else $
     mratfit[*,ired,ifinal] = bspline_valu(loglam[*,ired,ifinal], sset_r)

   !p.multi = [0,1,2]
   explist = expnum[uniq(expnum, sort(expnum))]
   colorvec = ['default','red','green','blue','cyan','magenta','grey']
   xrange = 10^minmax(loglam[*,*,ifinal])
   ii = where(mrativar[*,*,ifinal] GT 0, ct)
   if (ct GT 1) then yrange = minmax((mratfit[*,*,ifinal])[ii]) $
    else yrange = minmax(mratfit[*,*,ifinal])
   plottitle = 'PLATE=' + string(plateid[0], format='(i4.4)') $
    + ' MJD=' + string(maxmjd, format='(i5.5)')
   for iexp=0, n_elements(explist)-1 do begin
      djs_plot, [0], [0], xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
       /nodata, xtitle='Wavelength [Ang]', ytitle='Flux-calib', $
       title=plottitle+' Exp #' + string(explist[iexp], format='(i8.8)')
      kk = where(expnum EQ explist[iexp], kct) ; blue+red files for this exp
      for j=0, nfinal-1 do begin
         thiscolor = colorvec[j MOD n_elements(colorvec)]
         for k=0, kct-1 do begin
            djs_oplot, 10^loglam[*,kk[k],ifinal[j]], $
;             mratio[*,kk[k],ifinal[j]] * flatarr[*,kk[k],ifinal[j]], $
             mratio[*,kk[k],ifinal[j]], $
             psym=3, color=thiscolor
;            djs_oplot, 10^loglam[*,kk[k],ifinal[j]], $
;;             mratfit[*,kk[k],ifinal[j]] * flatarr[*,kk[k],ifinal[j]], $
;             mratfit[*,kk[k],ifinal[j]], $
;             color=thiscolor
         endfor
         djs_xyouts, 0.9*xrange[0]+0.1*xrange[1], $
          yrange[0] + (j+1)*(yrange[1]-yrange[0])/(nfinal+1), $
          'Fiber '+strtrim(iphoto[ifinal[j]]+(spectroid[0]-1)*nfiber,2), $
          color=thiscolor
      endfor
   endfor
   !p.multi = 0

   ;----------
   ; Construct the final (B-splined) flux-calibration vectors

   for ifile=0, nfile-1 do begin
      splog, prelog=fileandpath(objname[ifile])
      ; Is this a blue CCD?
      ii = where(ifile EQ iblue, ct)
      if (ct EQ 1) then begin
         thisloglam = loglam[*,ifile,ifinal]
         thisset = sset_b
         thisflatarr = flatarr[*,iblue[ii],ifinal]
         thispres = pres_b[*,ii,ifinal]
      endif

      ; Is this a red CCD?
      ii = where(ifile EQ ired, ct)
      if (ct EQ 1) then begin
         thisloglam = loglam[*,ifile,ifinal]
         thisset = sset_r
         thisflatarr = flatarr[*,ired[ii],ifinal]
         thispres = pres_r[*,ii,ifinal]
      endif

      thismratio = mratio[*,ifile,ifinal]
      thismrativar = mrativar[*,ifile,ifinal]
      if (tag_exist(thisset,'NPOLY')) then $
       x2 = airmass[*,ifile,iphoto[ifinal]] $
      else $
       x2 = 0

      ; Evaluate the B-spline for the stars at their measured wavelengths
      ; in this exposure, then modulated by the mean FLATARR
      ; for the stars in this exposure.
      ; We re-fit the B-spline to exactly recover what we had before,
      ; just modulated by the lower-order polynomial FLATARR.

      logmin = min(thisloglam[where(thismrativar GT 0)], max=logmax)
      tmploglam = wavevector(logmin, logmax)
      flatarr_mean = 0 * tmploglam
      for i=0L, nfinal-1 do $
       flatarr_mean = flatarr_mean $
        + poly(tmploglam-3.5d0, thispres[*,0,i]) / nfinal

      ; Rather than re-generating the B-spline, I'll simply cheat and
      ; multiply the B-spline coefficients at their wavelengths.
      ; This is a bit of a hack, since there are NORD more values
      ; for FULLBKPT than there are for COEFF, so there's not exactly
      ; a 1-to-1 mapping between the two.
      tmpmult = interpol(flatarr_mean, tmploglam, thisset.fullbkpt)
      tmpmult = tmpmult[1:n_elements(thisset.fullbkpt)-thisset.nord]
      if (keyword_set(x2)) then begin
         for ipoly=0, thisset.npoly-1 do $
          thisset.coeff[ipoly,*] = thisset.coeff[ipoly,*] * tmpmult
         for ipoly=0, thisset.npoly-1 do $
          thisset.icoeff[ipoly,*] = thisset.icoeff[ipoly,*] * tmpmult
      endif else begin
         thisset.coeff = thisset.coeff * tmpmult
         thisset.icoeff = thisset.icoeff * tmpmult
      endelse

      if (keyword_set(x2)) then begin
         x2_min = min(x2, max=x2_max)
         splog, 'Exposure ', objname[ifile], $
          ' spans airmass range ', x2_min, x2_max
         tmpflux1 = bspline_valu(tmploglam, thisset, x2=x2_min+0*tmploglam)
         tmpflux2 = bspline_valu(tmploglam, thisset, x2=x2_max+0*tmploglam)
      endif else begin
         tmpflux1 = bspline_valu(tmploglam, thisset)
         tmpflux2 = 0
      endelse

      ;----------
      ; Make plots of the spectro-photometry data for this exposure only,
      ; overplotting the global fit to all exposures in red.

      ; The following info is just used for the plot title
      plottitle = 'PLATE=' + string(plateid[ifile], format='(i4.4)') $
       + ' MJD=' + string(mjd[ifile], format='(i5.5)') $
       + ' Spectro-Photo Calib for ' + camname[ifile] + '-' $
       + string(expnum[ifile], format='(i8.8)')

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
      if (tag_exist(thisset,'NPOLY')) then x2 = airmass[*,ifile,*] $
       else x2 = 0
      calibimg = float( bspline_valu(loglam1, thisset, x2=x2) )

      ; Set to zero any pixels outside the known flux-calibration region
      qbad = loglam1 LT logmin OR loglam1 GT logmax
      ibad = where(qbad, nbad)
      if (nbad GT 0) then calibimg[ibad] = 0

      minval = min(calibimg[where(qbad EQ 0)], max=maxval)
      if (maxval/minval GT 20. OR minval LT 0) then $
       splog, 'WARNING: Min/max fluxcalib = ', minval, maxval $
      else $
       splog, 'Min/max fluxcalib = ', minval, maxval

      ;----------
      ; Write the output file

      ; Put the Kurucz models for this exposure in the output structure
      kindx.modelflux = reform(modflux[*,ifile,*], npix, nphoto)

      calibfile = djs_filepath(string(camname[ifile], expnum[ifile], $
       format='("spFluxcalib-", a2, "-", i8.8, ".fits")'), $
       root_dir=combinedir)
      mwrfits, calibimg, calibfile, hdr, /create
      mwrfits, thisset, calibfile
      mwrfits, kindx, calibfile
      spawn, ['gzip','-f',calibfile], /noshell
      splog, prelog=''
   endfor

   return
end
;------------------------------------------------------------------------------
