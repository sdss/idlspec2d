;+
; NAME:
;   myfluxcalib
;
; PURPOSE:
;   Solve for flux-calibration vectors, one for each camera (blue and red)
;   on a single spectrograph.
;
; CALLING SEQUENCE:
;   myfluxcalib, filename, calibfile, colors=, adderr=
;
; INPUTS:
;   filename   - Name(s) of input file names, typically one blue and
;                one red file.  The pluggings need not be the same,
;                but there must be at least one good SPECTROPHOTO or
;                REDDEN standard in common.
;   calibfile  - Name(s) of output calibration files, one per FILENAME
;
; REQUIRED KEYWORDS:
;   colors     - The name(s) of each camera, 'b' for blue or 'r' for red;
;                typically, we pass both with COLORS=['b','r']
;
; OPTIONAL INPUTS:
;   adderr     - Fractional errors to add in quadrature to measured errors
;                at each point in each spectrum; a good value is 0.03 for 3%.
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   bspline_iterfit()
;   bspline_valu()
;   correct_dlam
;   divideflat
;   djs_diff_angle()
;   djs_filepath()
;   djs_maskinterp()
;   djs_median()
;   djs_oplot
;   djs_plot
;   fibermask_bits()
;   fileandpath()
;   mrdfits()
;   mwrfits
;   pca_solve()
;   readcol
;   splog
;   sxaddpar
;   sxpar()
;   traceset2xy
;
; INTERNAL SUPPORT ROUTINES
;   fluxfit()
;   qgoodfiber()
;   resortplugmap()
;
; REVISION HISTORY:
;   08-Sep-2000  Written by D. Schlegel & S. Burles
;-
;------------------------------------------------------------------------------
; Re-sort the 2nd plug-map to match that of the first file, so
; that PLUGMAP2[INDX] = PLUGMAP1
function resortplugmap, plugmap1, plugmap2

   nfiber = n_elements(plugmap1)
   indx = lonarr(nfiber)

   for ifiber=0, nfiber-1 do begin
      adist = djs_diff_angle(plugmap1.ra, plugmap1.dec, $
       plugmap2[ifiber].ra, plugmap2[ifiber].dec, units='degrees')
      indx[ifiber] = $
       (where(adist LT 2./3600. AND strtrim(plugmap1.objtype,2) NE 'NA'))[0]
   endfor

   return, indx
end

;------------------------------------------------------------------------------
; Actually to the fit of the observed spectrum to that of an F8 star.
function fluxfit, loglam, objflux, objivar, color=color, refmag=refmag

   wave = 10^loglam
   logmin = min(loglam)
   logmax = max(loglam)

   ;----------
   ; Read the spectrum of an F8 star

   f8file = filepath('f8vspline.dat', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, f8file, f8wave, f8flux

   ;----------
   ; Scale the flux
   ;       at v=0, flux at 5556A is 956 photons/cm^2/A/s
   ;               = 3.52e-9 ergs/cm^2/s/A
   ;       scale to r', with 10^(-r'/2.5)
   ;       and return in units to 1e-17 ergs/cm/s/A
   ;       so factor in exponent is 10^((21.37 - r')/2.5)
   ;
   ;       AB magnitude, the scaling this assumes
   ;       the AB magnitude of BD+17 is 9.4 at 5560
   ;
   ;       we then use f_nu = 10^-0.4*(AB+48.6)
   ;       and then f_lambda = f_nu * c / (lambda)^2
   ; 
   ;       c is 3.0e18 Ang/s
   ;       lambda is in Ang

   oldzeropoint = 21.37

   ; The original file f8v.dat gives 21.2985, with continuum only:
   ; Absorption in r band only gives a 0.003 difference!  We will
   ; ignore for the time being
   ; 35.4771 = alog10(1.0e17 * 3.0e18)

   fiducial_flux = filter_thru(f8flux * f8wave^2, waveimg=f8wave)
   zeropoint = -(alog10(fiducial_flux[2]) - 35.4771)*2.5 - 48.6

   if zeropoint LT 21.0 OR zeropoint GT 21.5 then begin
     splog, 'WARNING: has the flux calibration standard changed?'
     splog, 'Please address, for now, resetting to old zero point'
     zeropoint = oldzeropoint
   endif 

   scalefac = 10.^((zeropoint - refmag) / 2.5)
   f8flux = f8flux * scalefac

   ;----------
   ; Divide the input spectrum by that of the F8 star

   intrinspl = spl_init(alog10(f8wave), f8flux)
   f8spline = spl_interp(alog10(f8wave), f8flux, intrinspl, loglam)
   fitflux = objflux / f8spline

   ;----------
   ; Mask out around absorption lines

   absfile = filepath('f8v.abs', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, absfile, absmin, absmax

   abmask = bytarr(n_elements(objflux)) + 1b
   for i=0, n_elements(absmin)-1 do $
    abmask = abmask * (wave LT absmin[i] OR wave GT absmax[i])

   ;----------
   ; Select break points for spline

   if (color EQ 'b') then bkptfile = filepath('blue.bkpts', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   if (color EQ 'r') then bkptfile = filepath('red.bkpts', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, bkptfile, allbkpts
   allbkpts = alog10(allbkpts) ; Convert to log-10 Angstroms

   ibk = where(allbkpts GE logmin AND allbkpts LE logmax, nbk)
   if (nbk LT 4) then $
    message, 'Error selecting break points'
   allbkpts = [logmin, allbkpts[ibk[1:nbk-2]], logmax]

   ;----------
   ; Do the spline fit

   indx = where(abmask)
   if (keyword_set(objivar)) then invvar = objivar[indx] $
    else invvar = 0
   sset = bspline_iterfit(loglam[indx], fitflux[indx], $
    invvar=invvar, nord=4, bkpt=allbkpts, upper=3, lower=3, $
    maxrej=ceil(0.05*n_elements(objflux)))

   ;----------
   ; QA plot

   yplot = bspline_valu(loglam,sset)
   djs_plot, [min(wave)-100,max(wave)+100], [0,1.1*max(yplot)], /nodata, $
    /xstyle, /ystyle, $
    xtitle='\lambda [A]', ytitle='Counts / (10^{-17}erg/cm/s/Ang)', $
    title='Spectro-Photo PCA for ' + color + '-camera'
   ymid = total(0.5 * !y.crange)
   for i=0, n_elements(absmin)-1 do begin
      djs_oplot, [absmin[i],absmin[i]], !y.crange, color='blue'
      djs_oplot, [absmax[i],absmax[i]], !y.crange, color='blue'
      djs_oplot, [absmin[i],absmax[i]], [ymid,ymid], color='blue'
   endfor
   djs_oplot, 10^loglam, fitflux
   djs_oplot, wave, yplot, color='red'
   djs_oplot, 10^allbkpts, bspline_valu(allbkpts,sset), psym=4, color='red'

   return, sset
end

;------------------------------------------------------------------------------
function qgoodfiber, fibermask
   qgood = ((fibermask AND fibermask_bits('NOPLUG')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADTRACE')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADFLAT')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADARC')) EQ 0) $
       AND ((fibermask AND fibermask_bits('MANYBADCOLUMNS')) EQ 0) $
       AND ((fibermask AND fibermask_bits('NEARWHOPPER')) EQ 0) $
       AND ((fibermask AND fibermask_bits('MANYREJECTED')) EQ 0)
   return, qgood
end

;------------------------------------------------------------------------------
pro myfluxcalib, filename, calibfile, colors=colors, adderr=adderr

   dloglam = 1.0d-4 ; ???

   nfile = n_elements(filename)
   if (n_elements(calibfile) NE nfile) then $
    message, 'Dimensions of FILENAME and CALIBFILE do not agree'

   ;----------
   ; Assign scores to each object based upon its color relative
   ; to the color BD+17 4708.

   bd17mag = [10.56, 9.64, 9.35, 9.25, 9.23]
   bd17color = bd17mag[0:3] - bd17mag[1:4]

   ;----------
   ; Read the plug maps and masks and decide which star(s) to use

   for ifile=0, nfile-1 do begin
      thismask = mrdfits(filename[ifile],2)
      thisplug = mrdfits(filename[ifile],5)
      goodfiber = qgoodfiber(thismask[0,*])

      if (ifile EQ 0) then begin
         ; Save the first plug-map for comparison to the others
         plugmap = thisplug
         plugindx = lindgen(n_elements(thisplug))

         ; Save the first good-fiber mask
         goodmask = goodfiber

         ; Assign scores based upon colors and magnitudes
         colordiff = thisplug.mag[0:3] - thisplug.mag[1:4] 
         for i=0, 3 do $
          colordiff[i,*] = colordiff[i,*] - bd17color[i]
         colordiff = sqrt(total(colordiff^2,1))
         colorscore = thisplug.mag[2] + 40.0 * colordiff ; Lower score is better
      endif else begin
         ; Re-sort the plug-map to match that of the first file
         indx = resortplugmap(plugmap, thisplug)
         plugindx = [ [plugindx], [indx] ]

         ; Multiply this good-fiber mask with the others,
         ; setting to zero if the object does not exist in this plugging.
         goodmask = goodmask * goodfiber[indx>0] * (indx NE -1)
      endelse
   endfor

   ;----------
   ; Select all flux-calibration stars

   iphoto = where( ( strtrim(plugmap.objtype) EQ 'SPECTROPHOTO_STD' OR $
    strtrim(plugmap.objtype) EQ 'REDDEN_STD') AND goodmask, nphoto)

   if (nphoto EQ 0) then $
    message, 'No SPECTROPHOTO or REDDEN stars for flux calibration'

   ;----------

   nnew = lonarr(nfile)

   for ifile=0, nfile-1 do begin

      ;----------
      ; Read in the flux, errors, and wavelengths

      objflux = mrdfits(filename[ifile], 0, hdr)
      objivar = mrdfits(filename[ifile], 1)
      objmask = mrdfits(filename[ifile], 2)
      wset = mrdfits(filename[ifile], 3)
      traceset2xy, wset, xx, loglam

      ;----------
      ; Read header info

      cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
      expstr = string(sxpar(hdr,'EXPOSURE'), format='(i8.8)')

      ;----------
      ; Do not fit where the spectrum may be dominated by sky-sub residuals.

      objivar = skymask(objivar, objmask)
objmask = 0 ; Free memory

      ;----------
      ; Add an additional error term equal to ADDERR of the flux.

      if (keyword_set(adderr)) then begin
         gmask = objivar NE 0 ; =1 for good points
         objivar = 1.0 / ( 1.0/(objivar + (1-gmask)) $
          + (adderr * (objflux>0))^2 ) * gmask
      endif

      ;----------
      ; Make a map of the size of each pixel in delta-(log10-Angstroms),
      ; and re-normalize the flux to ADU/(dloglam)

      correct_dlam, objflux, objivar, wset, dlam=dloglam

      ;----------
      ; Apply flux-correction factor between spectro-photometric exposure
      ; and this exposure.

      spectroid = strmid(cameras,1,1)
      cc = fileandpath(calibfile[ifile], path=combinedir)
 
      corrset = mrdfits($
       djs_filepath('spFluxcorr-'+expstr+'-'+spectroid+'.fits', $
        root_dir=combinedir), 1)
      traceset2xy, corrset, loglam, corrimg

      divideflat, objflux, objivar, 1.0/corrimg, $
       minval=0.1/median(corrimg)

      ;----------
      ; Re-bin the spectro-photo stars to the same wavelength mapping

      newloglam = wavevector(min(loglam[*,indx]), max(loglam[*,indx]), $
       binsz=dloglam)
      nnew[ifile] = n_elements(newloglam)

      newflux = fltarr(nnew[ifile],nphoto)
      newivar = fltarr(nnew[ifile],nphoto)

      for iobj=0, nphoto-1 do begin
         combine1fiber, loglam[*,iphoto[iobj]], $
          objflux[*,iphoto[iobj]], objivar[*,iphoto[iobj]], $
          newloglam=newloglam, binsz=dloglam, newflux=flux1, newivar=ivar1
         newflux[*,iobj] = flux1
         newivar[*,iobj] = ivar1
      endfor

      newflux = djs_maskinterp(newflux, newivar EQ 0, iaxis=0, /const)

      ;----------
      ; Append the data from all files together

      if (ifile EQ 0) then begin
         allloglam = newloglam
         allflux = newflux
         allivar = newivar
      endif else begin
         allloglam = [allloglam, newloglam]
         allflux = [allflux, newflux]
         allivar = [allivar, newivar]
      endelse

   endfor

   ;----------
   ; Find the PCA solution for all files simultaneously
   ; Use a rejection scheme that rejects up to 1% of all spectral bins
   ; (~40 pixels) for each 10% of the spectrum (~800 pixels).

   npix = (size(allflux,/dimens))[0]

   ;-------------------------
   ; With so few spectra for a PCA solution, try to reject some of the
   ; most deviant points first.  Do this by re-normalizing all the spectra
   ; to the same mean, and replace the discrepant points with the mean
   ; of the good points (and mask them).

   ; Find the normalization factor for each spectrum, MULTFAC.
   multfac = fltarr(nphoto)
   for iobj=0, nphoto-1 do $
    multfac[iobj] = $
     median( djs_median(allflux[0:nnew[0]-1,iobj], $
      width=99, boundary='reflect')) + $
     median( djs_median(allflux[nnew[0]:npix-1,iobj], $
      width=99, boundary='reflect'))

   ; At each wavelength, find the median normalized flux value.
   medvec = fltarr(npix)
   for ipix=0, npix-1 do begin
      normflux = allflux[ipix,*] / multfac
;      normivar = allivar[ipix,*] * multfac^2
      igood = where(allivar[ipix,*] NE 0, ngood)
      if (ngood EQ 0) then medvec[ipix] = median([normflux]) $
       else medvec[ipix] = median([normflux[igood]])
   endfor

   ; Now for each spectrum, mask deviant points.
   ; Look for locally deviant points in the ratio of each object spectrum
   ; to the median spectrum.  This allows there to be a different color
   ; for each object, showing up as a slope in that ratio.
   ; Only reject points that are locally deviant by more than 50%.
   nreject = 0
   for iobj=0, nphoto-1 do begin
      ratvec = allflux[*,iobj] / multfac[iobj] - medvec
      ratvec[0:nnew[0]-1] = ratvec[0:nnew[0]-1] $
       - djs_median(ratvec[0:nnew[0]-1], $
       width=99, boundary='reflect')
      ratvec[nnew[0]:npix-1] = ratvec[nnew[0]:npix-1] $
       - djs_median(ratvec[nnew[0]:npix-1], $
       width=99, boundary='reflect')
      ibad = where(abs(ratvec) GT 0.50 * medvec, nbad) ; <-- 50% threshhold
      ; Mask the deviant points and replace with the median value
      if (nbad GT 0) then begin
         allflux[ibad,iobj] = medvec[ibad] * multfac[iobj]
         allivar[ibad,iobj] = 0
         nreject = nreject + nbad
      endif
   endfor
   splog, 'Rejected points from median spectrum = ', nreject

   ; Use a large number of iterations below to be sure that we converge
   ; in de-weighting very bad points.
   pcaflux = pca_solve(allflux, allivar, $
    niter=30, nkeep=1, usemask=usemask, eigenval=eigenval, acoeff=acoeff, $
    maxiter=5, upper=5, lower=5, maxrej=ceil(0.01*npix), $
    groupsize=ceil(npix/5.))

   maxmask = max(usemask)
   if (maxmask LE 3) then minuse = 1 $
    else if (maxmask EQ 4) then minuse = 2 $
    else minuse = 3

   ;----------
   ; Determine the reference (r-band) magnitude for the PCA spectrum.
   ; We get an estimate from each of the input spectro-photo spectra.
   ; Then choose the median of these estimates.

   refmag = plugmap[iphoto].mag[2] + 2.5 * alog10(acoeff)
   splog, 'Estimates of spectro-photo PCA r-mag = ', refmag
   splog, 'Dispersion of spectro-photo PCA r-mag = ', $
    (nphoto EQ 1 ? 0 : stdev(refmag))
   refmag = median([refmag], /even)

   ;----------
   ; Set up for QA plots

   !p.multi = [0,1,nfile]

   ;----------
   ; Do the actual fits

   i1 = 0
   for ifile=0, nfile-1 do begin

      i2 = i1 + nnew[ifile] - 1

      ;----------
      ; Find indices of good data points for this eigen-spectrum

      indx = lindgen(i2-i1+1) + i1
      indx = indx[ where(usemask[indx] GE minuse) ]

      ; Compare this mean (PCA) spectrum to an F star spectrum,
      ; w/out using any errors but allowing some rejection.
      calibset = fluxfit(allloglam[indx], pcaflux[indx], $
       color=colors[ifile], refmag=refmag)

      ;----------
      ; Create header cards describing the fit range

      hdr = ['']
      wavemin = 10.^min(allloglam[indx])
      wavemax = 10.^max(allloglam[indx])
      sxaddpar, hdr, 'WAVEMIN', wavemin
      sxaddpar, hdr, 'WAVEMAX', wavemax

      mwrfits, 0, calibfile[ifile], hdr, /create
      mwrfits, calibset, calibfile[ifile]

      i1 = i2 + 1
   endfor

   !p.multi = 0

   return
end
;------------------------------------------------------------------------------
