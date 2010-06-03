;+
; NAME:
;   solvefilter
;
; PURPOSE:
;   Solve for the 2.5-m imaging filter curves by using the spectra.
;
; CALLING SEQUENCE:
;   solvefilter, [ filttype=, filternum=, plate=, mjd=, $
;    starerr=, qsoerr=, wavemin=, wavemax=, magrej=, sncut=, maxiter=, $
;    fluxpath=, value=, fixed= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   filttype   - Type of functional form for the filter curve.  Options are:
;                  'sdss': Modified SDSS filter curve (default).  3 params.
;                  'tanh': Function with tanh() shape at edges. 5 params.
;   filternum  - Filter number, 1=g, 2=r, 3=i; default to 3.
;   plate      - Plate number(s); if not specified, then select all DR1 plates
;                with number > 431.
;   mjd        - MJD for each PLATE.
;   starerr    - Fractional error to add in quadrature to photometric errors
;                for stars; default to 0.05; if <=0, then do not use stars.
;   qsoerr     - Fractional error to add in quadrature to photometric errors
;                for QSOs; default to 0.15; if <=0, then do not use QSOs.
;   wavemin    - Minimum wavelength for spectra during computation;
;                default to 3800 Ang.
;   wavemax    - Maximum wavelength for spectra during computation;
;                default to 9300 Ang.
;   magrej     - Reject any objects where the raw photo vs. spectro magnitude
;                difference is more than MAGREJ from the median difference
;                for that plate.  This will reject wild outliers, which are
;                often objects where there is a bright blend but the PHOTO
;                flux is only for a fainter child.  Default value is 0.5 mag.
;                Or, QSOs that have varied.
;   sncut      - Minimum SN_MEDIAN (median S/N per pixel) for spectroscopic
;                objects used in sample; default to 2.0
;   maxiter    - Maximum number of iterations in call to MPFIT(); default to 200
;   fluxpath   - Path name for spPlate files used for reading the spectra.
;                The spZ files are still read from $BOSS_SPECTRO_REDUX/$PLATE
;                regardless of this keyword.
;   value      - Initial guess values for the fit parameters.  The default
;                values are chosen to closely match the Gunn Jun-2001 curves.
;   fixed      - A vector of elements set to 0 for each parameter to be fit,
;                and 1 for each parameter value to fix.  Default to fitting
;                all parameters.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   For FILTTYPE='tanh', the following function response(loglam) is fit:
;     m = (a[4] - 1) / (a[1] - a[0])
;     b = 1 - m * a[0]
;     Filter = tahn((loglam - a[0])*a[2]) + 1)
;            * tahn((a[1] - loglam)*a[3]) + 1)
;            * (m * loglam + b)
;   There are a total of five a[] parameters above.
;   The default initial-guess values for g-band are:
;     a = [alog10(3950), alog10(5325), 100, 140, 2.0]
;   The default initial-guess values for r-band are:
;     a = [alog10(5580), alog10(6750), 130, 160, 1.2]
;   The default initial-guess values for i-band are:
;     a = [alog10(6915), alog10(8210), 150, 220, 0.55]
;
;   For FILTTYPE='sdss', we fit a modified version of the Gunn Jun-2001
;   filter curves.  The wavelength scale is remapped with a shift and
;   rescaling, which has the effect of moving the filter edges.  The
;   filter shape is also multiplied by a function that is linear in
;   log-wavelength, which has the effect of changing the broad-band
;   slope of the filter.  Given a filter curve Gunnfilt(loglam), it is
;   re-mapped as follows:
;     slopeterm = (loglam - 3.5)^theta[2]
;     Filter(loglam) = Gunnfilt((loglam + a[1]) * a[0]) * slopeterm
;   The default initial-guess values for the three a[] parameters are always:
;     a = [0, 1.0, 0.01]
;   The effect of a positive a[2] makes the filter slope more upwards
;   with wavelength, and a negative a[2] makes it slope more downwards.
;
;   Iteratively solve for the SDSS 2.5-m filter curves, using one of
;   the several possible parameterizations as specified by FILTTYPE.
;   For each possible filter curve, we regress the spectroscopic magnitude
;   vs. the photometric magnitude.  In order to not be sensitive to
;   photometric calibration errors, each group of objects (same RUN, RERUN,
;   CAMCOL, PLATE, SPECTROGRAPHID) is allowed to have a floating zero-point
;   offset.  These offsets are plotted as MAGOFFSET in one of the final plots.
;
;   We only use spectroscopically-confirmed stars and QSOs, not any galaxies.
;   Anything targetted as a galaxy is rejected.  Any blended objects are
;   rejected.  Anything with the following bits set is rejected:
;   OBJECT2_SATUR_CENTER, OBJECT2_INTERP_CENTER, OBJECT2_PSF_FLUX_INTERP.
;   These last cuts remove stars with CRs in the core on the i-band images
;   that were targetted as QSOs -- this was actually due to a bug in PHOTO
;   that called stars CRs when the seeing was too good, then the incorrect
;   CR-removal made the i-band 2 mags fainter, and these things were targetted
;   as QSOs.
;
;   All calculations are done in vacuum wavelengths, but then converted
;   to air wavelengths at the end.
;
; EXAMPLES:
;   Solve for the i-band filter using Tremonti's re-reductions of the spectra:
;     IDL> solvefilter, filternum=3, fluxpath='/scr/wire50/cat/recalib/kurucz'
;
; BUGS:
;   I should average together more telluric spectra for better S/N.
;   I should use the extinction coeff for each imaging night.
;   Do Christy's flux-calibrations improve things?
;   Do any objects argue for light leaks?
;   Should we more heavily weight the QSOs?
;
; DATA FILES:
;   $IDLUTILS/data/filters/sdss_jun2001_$FILTER_atm.dat
;   $BOSS_SPECTRO_REDUX/0432/spFrame-r2-00007466.fits*  (for telluric-correction)
;   $BOSS_SPECTRO_REDUX/$PLATE/spPlate-$PLATE-$MJD.fits
;   $BOSS_SPECTRO_REDUX/$PLATE/spZbest-$PLATE-$MJD.fits
;   $BOSS_SPECTRO_REDUX/plates/tsObj*-$PLATE.fit
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_int2bin()
;   djs_maskinterp()
;   djs_oplot
;   djs_plot
;   djs_xyouts
;   headfits()
;   idlspec2d_version()
;   mpfit()
;   mrdfits()
;   readcol
;   readspec
;   sdss_run2mu()  (in photoop product)
;   skymask()
;   splog
;   tai2airmass()
;   traceset2xy
;   wavevector()
;
; INTERNAL SUPPORT ROUTINES:
;   solvefiltshape()
;   solvefiltfn()
;
; REVISION HISTORY:
;   05-Nov-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
forward_function mpfit, solvefiltfn

;------------------------------------------------------------------------------
; Construct the filter curve corresponding to this set of parameters
; Return the Gunn filter curve if THETA is not set.

function solvefiltshape, theta, loglam

   common com_solvefilt, groupnum, bigloglam, taugunn, $
    bigflux, photoflux, photoinvsig, tauextinct, ntot, nbigpix, $
    airmass, gunnfilt, filternum, filttype, groupratio, spectroflux

   if (NOT keyword_set(theta)) then $
    return, gunnfilt[*,filternum]

   case filttype of
   'tanh': begin
      mm = (theta[4] - 1.d0) / (theta[1] - theta[0])
      bb = 1.d0 - mm * theta[0]
      fcurve = (tanh((loglam - theta[0])*theta[2]) + 1.d0) $
             * (tanh((theta[1] - loglam)*theta[3]) + 1.d0) $
             * (mm * loglam + bb)
      end
   'sdss': begin
      slopeterm = (bigloglam - 3.5d0)^theta[2]
      linterp, bigloglam, gunnfilt[*,filternum] * slopeterm, $
       (bigloglam - theta[0]) * theta[1], fcurve
      end
   endcase

   if (total(finite(fcurve)) NE n_elements(fcurve)) then $
    message, 'NaN in filter shape'

   fcurve = fcurve > 0

   return, fcurve
end
;------------------------------------------------------------------------------
function solvefiltfn, theta

   common com_solvefilt, groupnum, bigloglam, taugunn, $
    bigflux, photoflux, photoinvsig, tauextinct, ntot, nbigpix, $
    airmass, gunnfilt, filternum, filttype, groupratio, spectroflux

   ngroup = max(groupnum) + 1
   groupratio = dblarr(ngroup)

   ; Construct the filter curve corresponding to this set of parameters
   fcurve = solvefiltshape(theta, bigloglam)
   sumfilt = total(fcurve)

   ; Change from f_lambda to f_nu
   flambda2fnu = 10^(2*bigloglam) / 2.99792d18 * 10^((48.6 - 2.5*17.)/2.5)

   ; Loop through each group of spectra, integrate over the filter curve,
   ; and minimize the spectro/photo flux normalization for that group.
   spectroflux = dblarr(ntot)
   leftall = dblarr(ntot)
   rightall = dblarr(ntot)
   for igroup=0L, ngroup-1 do begin
      indx = where(groupnum EQ igroup, nthis)

      ; Get the extinction curve for these objects.  If THETA is undefined,
      ; then use the Gunn extinction curve.
      meanair = mean(airmass[indx])
;      if (NOT keyword_set(theta)) then fextinct = exp(-meanair * taugunn) $
      if (NOT keyword_set(theta)) then fextinct = 1. $ ; ???
       else fextinct = exp(-meanair * tauextinct)
      fmult = fcurve * flambda2fnu * fextinct

      spectroflux[indx] = $
       total(bigflux[*,indx] * rebin(fmult,nbigpix,nthis), 1) $
       / (sumfilt + (sumfilt LE 0))
      leftval = spectroflux[indx] * photoinvsig[indx]
      rightval = photoflux[indx] * photoinvsig[indx]
      denom = total(leftval^2)
      if (denom GT 0) then groupratio[igroup] = total(leftval * rightval) / denom $
       else groupratio[igroup] = 1
      leftval = leftval * groupratio[igroup]
      leftall[indx] = leftval
      rightall[indx] = rightval
   endfor

   ; Return a vector of all the chi's.
   return, leftall - rightall
end
;------------------------------------------------------------------------------
pro solvefilter, filttype=filttype1, filternum=filternum1, $
 plate=plate, mjd=mjd, starerr=starerr, qsoerr=qsoerr, $
 wavemin=wavemin, wavemax=wavemax, $
 magrej=magrej, sncut=sncut, maxiter=maxiter, fluxpath=fluxpath, $
 value=value, fixed=fixed

   ;----------
   ; Set common block

   common com_solvefilt, groupnum, bigloglam, taugunn, $
    bigflux, photoflux, photoinvsig, tauextinct, ntot, nbigpix, $
    airmass, gunnfilt, filternum, filttype, groupratio, spectroflux

   if (keyword_set(filttype1)) then filttype = filttype1 $
    else filttype = 'sdss'
   if (n_elements(filternum1) NE 0) then filternum = filternum1 $
    else filternum = 3 ; Default to i-band
   if (n_elements(starerr) EQ 0) then starerr = 0.05
   if (n_elements(qsoerr) EQ 0) then qsoerr = 0.15
   if (NOT keyword_set(wavemin)) then wavemin = 3800.d0
   if (NOT keyword_set(wavemax)) then wavemax = 9300.d0
   if (NOT keyword_set(magrej)) then magrej = 0.5
   if (NOT keyword_set(sncut)) then sncut = 2.0
   if (NOT keyword_set(maxiter)) then maxiter = 200
   wcovcut = 0.32
   filtname = ['u','g','r','i','z']
   if (starerr LE 0) then usestars = 0 $
    else usestars = 1
   if (qsoerr LE 0) then useqsos = 0 $
    else useqsos = 1

   if (filternum LT 1 OR filternum GT 3) then $
    message, 'I only can cope with FILTERNUM=1,2, or 3'

   t0 = systime(1)

   bigloglam = wavevector(alog10(wavemin), alog10(wavemax))
   nbigpix = n_elements(bigloglam)

   ;----------
   ; Read Gunn's measurement of the filter -- including the telluric bands.

   gunnfilt = dblarr(nbigpix,n_elements(filtname))
   for ifilt=0, n_elements(filtname)-1 do begin
      filename = filepath('sdss_jun2001_'+filtname[ifilt]+'_atm.dat', $
       root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
      readcol, filename, fwave1, fthru1, fthru2, fthru3, fext1, /silent

      ; Convert wavelengths to vacuum.
      airtovac, fwave1

      fthru = 0.5 * (fthru1 + fthru2) / fext1 ; Average these two columns
      linterp, alog10(fwave1), fthru, bigloglam, gunnfilt1
      gunnfilt[*,ifilt] = gunnfilt1

      ; Assemble the data to get the atmospheric extinction curve
      if (ifilt EQ 0) then begin
         fwave = fwave1
         fext = fext1
      endif else begin
         fwave = [fwave, fwave1]
         fext = [fext, fext1]
      endelse
   endfor

   isort = uniq(fwave, sort(fwave))
   fwave = fwave[isort]
   fext = fext[isort]
   linterp, alog10(fwave), fext, bigloglam, gunnextinct
   taugunn = -alog(gunnextinct) / 1.3 ; Scale from 1.3 to 1.0 airmasses

   ;----------
   ; Set the structure to pass to MPFIT

   blankpar = {value:0.D, fixed:0, limited:[0b,0b], $
    limits:[0.D,0], mpmaxstep: 0.D}

   case filttype of
   'tanh': begin
      parinfo = replicate(blankpar, 5)
      if (filternum EQ 1) then $
       parinfo.value = [alog10(3950), alog10(5325), 100, 140, 2.0]
      if (filternum EQ 2) then $
       parinfo.value = [alog10(5580), alog10(6750), 130, 160, 1.2]
      if (filternum EQ 3) then $
       parinfo.value = [alog10(6915), alog10(8210), 150, 220, 0.55]
      parinfo.limited = [[0,1], [1,0], [1,0], [1,0], [1,0]]
      logmid = 0.5 * (parinfo[0].value + parinfo[1].value)
      parinfo.limits = [[0,logmid-5.d-4], [logmid+5.d-4,0], $
       [0.0,0], [0.0,0], [0,0]]
      parinfo.mpmaxstep = [5.d-4, 5.d-4, 10., 10., 0.05]
      end
   'sdss': begin
      parinfo = replicate(blankpar, 3)
      parinfo.value = [0.d-5, 1.000, 0.01]
      parinfo.mpmaxstep = [5.d-4, 0.01, 0.02]
      end
   else: message, 'Unknown FILTTYPE'
   endcase

   if (keyword_set(value)) then parinfo.value = value
   if (keyword_set(fixed)) then parinfo.fixed = fixed

   ;----------
   ; Construct the atmospheric extinction curve

   ; Start with a simple expression for the extinction
;   tausimple = 10^(12.4 - 3.6 * bigloglam)

   ; Start with Gunn's extinction curve, but interpolating
   ; over the telluric bands which will be replaced with
   ; high-resolution spectra of those features.
   bigwave = 10^bigloglam
   vactoair, bigwave
   tmask = (bigwave GT 6800 AND bigwave LT 7000) $
        OR (bigwave GT 7100 AND bigwave LT 7400) $
        OR (bigwave GT 7550 AND bigwave LT 7750) $
        OR (bigwave GT 8050 AND bigwave LT 8350)
   tausimple = djs_maskinterp(taugunn, tmask, /const)

   framefile = filepath('spFrame-r2-00007466.fits*', $
    root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir='0432')
   fcalibfile = filepath('spFluxcalib-r2-00007466.fits*', $
    root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir='0432')
   framefile = (findfile(framefile))[0]
   fcalibfile = (findfile(fcalibfile))[0]
   hdr = headfits(framefile)
   thisair = sxpar(hdr, 'AIRMASS')
   wset = mrdfits(framefile, 3)
   traceset2xy, wset, xx, tloglam
   telluric = mrdfits(fcalibfile, 0)
   linterp, tloglam[*,0], telluric[*,0], bigloglam, tellcorr
   tautelluric = -alog(tellcorr)/thisair ; Scale back to one airmass
   tauextinct = (tausimple + tautelluric) > 0

   ;----------
   ; Find the list of good plates if not provided.
   ; Default to using all DR1 plates with PLATE>431.

   if (keyword_set(plate)) then begin
      splog, 'Using user-supplied list of plates'
      plist = replicate(create_struct('plate', 0L, 'mjd', 0L), $
       n_elements(plate))
      plist.plate = plate
      if (keyword_set(mjd)) then begin
         plist.mjd = mjd
      endif else begin
         ; Determine the MJD for each plate number
         readspec, plate, replicate(1,n_elements(plate)), mjd=mjd1
         plist.mjd = mjd1
      endelse
   endif else begin
      splog, 'Find list of good plates'
      platelist, plist=plist
      iuse = where(strmatch(plist.public,'*DR1*') AND plist.plate GT 431 $
       AND strmatch(plist.status1d,'Done*'))
      plist = plist[iuse]
   endelse

   ;----------
   ; Make certain that the data exists for all these plates,
   ; especially if PATH is set to read the flux vectors from
   ; some non-standard directory.

   splog, 'Test existence of flux data for all plates'
   nplate = n_elements(plist)
   qkeep = bytarr(nplate)
   for iplate=0, nplate-1 do begin
      readspec, plist[iplate].plate, mjd=plist[iplate].mjd, path=fluxpath, $
       objhdr=objhdr
      if (keyword_set(objhdr)) then qkeep[iplate] = 1B
   endfor
   ikeep = where(qkeep, nplate)
   if (nplate EQ 0) then $
    message, 'No plate data found in specified path'
   plist = plist[ikeep]
   splog, 'Number of usable plates = ', nplate

   ;----------
   ; Loop over each plate, and accumulate groups of objects

   nplate = n_elements(plist)
   gcounter = 0L
   for iplate=0, nplate-1 do begin
      splog, 'Reading redshift data for plate #', plist[iplate].plate, $
       ' (', iplate+1, ' of ', nplate, ')'
      readspec, plist[iplate].plate, mjd=plist[iplate].mjd, $
       zans=zans, tsobj=tsobj, plug=plug, /silent

      ; Identify objects targetted as galaxies (we don't want these)
      qgalaxy = (plug.primtarget AND 2L^6+2L^7+2L^8+2L^5+2L^26) NE 0

      ; Identify objects that are not blends
      bflag = djs_int2bin(ulong(tsobj.objc_flags), ndigit=32)
      bflag2 = djs_int2bin(ulong(tsobj.objc_flags2), ndigit=32)
      qblend = transpose(bflag[3,*] EQ 1 AND bflag[6,*] EQ 0)
      qbright = transpose(bflag[1,*])
      qchild = transpose(bflag[4,*])
      qsingle = (qblend EQ 0) AND (qbright EQ 0) AND (qchild EQ 0)

      ; Find which objects have interpolated pixels near the center
      ; (we'll throw them out)
      qinterp = transpose(bflag2[11,*] OR bflag2[12,*] OR bflag2[15,*])

      ; Compute the airmass for each object
;      junk = sdss_run2mu(tsobj.run, tsobj.field, tai=tai)
;      airmass1 = tai2airmass(zans.plug_ra, zans.plug_dec, tai=tai)
;      if (min(airmass1) LT 0.99 OR max(airmass1) GT 3.5) then $
;       message, 'Invalid AIRMASS'
      airmass1 = tsobj.airmass[2]

      ; Group objects with the same run+rerun+camcol+plate+spectrographid
      idstring = string(tsobj.run) + string(tsobj.rerun) $
       + string(tsobj.camcol) + string(zans.plate) + string(plug.spectrographid)
      idlist = idstring[ uniq(idstring, sort(idstring)) ]
      ngroup = n_elements(idlist)

      ; Reject wild mag outliers, which are often objects where there
      ; is a bright blend but the PHOTO flux is only for a fainter child
      magdiff = - 2.5 * alog10(zans.spectroflux[filternum]) $
       - tsobj.psfflux[filternum]
      meddiff = median(magdiff)

      qstar = strmatch(zans.class,'STAR*')
      qqso = strmatch(zans.class,'QSO*')

      for igroup=0, ngroup-1 do begin
         indx = where(idstring EQ idlist[igroup])

         igood = where(zans[indx].zwarning EQ 0 $
          AND qsingle[indx] EQ 1 $
          AND qgalaxy[indx] EQ 0 $
          AND qinterp[indx] EQ 0 $
          AND (qstar[indx]*usestars OR qqso[indx]*useqsos) $
          AND magdiff[indx] GT meddiff-magrej $
          AND magdiff[indx] LT meddiff+magrej $
          AND zans[indx].wcoverage GT wcovcut $
          AND zans[indx].sn_median GT sncut, ngood)

         if (ngood GE 2) then begin
            if (gcounter EQ 0) then begin
               groupnum = replicate(gcounter, ngood)
               zall = zans[indx[igood]]
               tsall = tsobj[indx[igood]]
               airmass = airmass1[indx[igood]]
            endif else begin
               groupnum = [groupnum, replicate(gcounter, ngood)]
               zall = [zall, zans[indx[igood]]]
               tsall = [tsall, tsobj[indx[igood]]]
               airmass = [airmass, airmass1[indx[igood]]]
            endelse
            gcounter = gcounter + 1
         endif
      endfor
   endfor
   ntot = n_elements(zall)
   splog, 'Number of spectra = ', ntot
   splog, 'Number of groups = ', max(groupnum)+1

   ;----------
   ; Construct the big matrices

   bigflux = fltarr(nbigpix, ntot)

   ;----------
   ; Read in the actual spectra

   splog, 'Reading spectra'
   for iplate=0, nplate-1 do begin
      splog, 'Reading spectra for plate #', plist[iplate].plate, $
       ' (', iplate+1, ' of ', nplate, ')'
      indx = where(zall.plate EQ plist[iplate].plate $
       AND zall.mjd EQ plist[iplate].mjd, nthis)
      if (nthis GT 0) then begin
         ; Read in the spectra, and interpolate over bad points
         fiberid = zall[indx].fiberid
         readspec, plist[iplate].plate, mjd=plist[iplate].mjd, fiberid, $
          flux=objflux, loglam=loglam, invvar=objivar, $
          andmask=andmask, ormask=ormask, path=fluxpath, /align
         npix = n_elements(loglam)
         objivar = skymask(objivar, andmask, ormask)
         objflux = djs_maskinterp(objflux, objivar LE 0, iaxis=0, /const)

         if (bigloglam[0] LT loglam[0]) then begin
            i1 = (where(bigloglam GE loglam[0]))[0]
            j1 = 0L
         endif else begin
            i1 = 0L
            j1 = (where(loglam GE bigloglam[0]))[0]
         endelse
         ncopy = (nbigpix - i1) < (npix - j1)
         bigflux[i1:i1+ncopy-1,indx] = objflux[j1:j1+ncopy-1,*]
      endif
   endfor
objflux = 0
andmask = 0
ormask = 0
objivar = 0

   ;----------
   ; Decide upon the object counts and errors from PHOTO.

;   photoflux = 10.d0^(-tsall.psfcounts[filternum]/2.5)
;   photoflerr = tsall.psfcountserr[filternum] * abs(photoflux)
   photoflux = tsall.psfflux[filternum]
   photoflivar = tsall.psfflux_ivar[filternum]
   qgood = photoflivar GT 0

   ; Add an additional error term
   qstar = strmatch(zall.class,'STAR*')
   qqso = strmatch(zall.class,'QSO*')
;   photoinvsig = 1. / sqrt( photoflerr^2 $
;    + (qstar * starerr * abs(photoflux))^2 $
;    + (qqso * qsoerr * abs(photoflux))^2 )
   photoinvsig = qgood / sqrt( 1./(photoflivar+qgood-1) $
    + (qstar * starerr * abs(photoflux))^2 $
    + (qqso * qsoerr * abs(photoflux))^2 )

   ; Now convert these to AB flux, according to the numbers derived
   ; by Hogg on 13 Aug 2002.
   aboffsets = [-0.042, 0.036, 0.015, 0.013, -0.002]
   photoflux = photoflux * 10.d0^(-aboffsets[filternum]/2.5)

   ;----------
   ; Compute the chi^2 for the initial guess parameters

   origchi =  solvefiltfn(parinfo.value)
   dof = ntot - n_elements(parinfo) - ngroup
   origrchi2 = total(origchi^2) / dof
   splog, 'Start Chi2/DOF = ', origrchi2

   ;----------
   ; Do the actual fit to the filter curve

   t1 = systime(1)
   theta = mpfit('solvefiltfn', $
    parinfo=parinfo, perror=perror, maxiter=maxiter, $
    nfev=nfev, niter=niter, status=status)
   t2 = systime(1)
   chivec = solvefiltfn(theta)
   chi2pdof = total(chivec^2) / dof

   splog, 'Time for non-linear fitting = ', t2-t1, ' sec'
   splog, 'Number of iterations = ', niter
   splog, 'Number of function evaluations = ', nfev
   splog, 'Fit values = ', theta
   splog, 'Fit errors = ', perror
   splog
   splog, 'Median |chi| = ', median(abs(chivec))
   splog, 'Final Chi2/DOF = ', chi2pdof

   ;----------
   ; Identify the 10 most deviant points

   nworst = 10
   iworst = (reverse(sort(abs(chivec))))[0:nworst-1]

   ;----------
   ; Reconstruct the filters at 1.3 airmasses

   ; Gunn filter curve
   finitial = solvefiltshape(parinfo.value, bigloglam) * exp(-1.3 * tauextinct)
   finitial = finitial * mean(gunnfilt[*,filternum]) / mean(finitial)
   junk1 = solvefiltfn() ; Force evaluation of spectroflux
   initdiff1 = -2.5 * alog10(spectroflux / photoflux)
   initdiff2 = -2.5 * alog10(spectroflux * groupratio[groupnum] / photoflux)

   ; 1st guess filter curve
   fguess = solvefiltshape(parinfo.value, bigloglam) * exp(-1.3 * tauextinct)
   fguess = fguess * mean(gunnfilt[*,filternum]) / mean(fguess)
   junk1 = solvefiltfn(parinfo.value) ; Force evaluation of spectroflux
   guessdiff1 = -2.5 * alog10(spectroflux / photoflux)
   guessdiff2 = -2.5 * alog10(spectroflux * groupratio[groupnum] / photoflux)

   ; Best-fit filter curve
   fbest = solvefiltshape(theta, bigloglam) * exp(-1.3 * tauextinct)
   fbest = fbest * mean(gunnfilt[*,filternum]) / mean(fbest)
   fbest = fbest * (fbest GT 0) + 0.0 * (fbest LE 0) ; Get rid of values -0.00
   junk2 = solvefiltfn(theta) ; Force evaluation of spectroflux
   magdiff1 = -2.5 * alog10(spectroflux / photoflux)
   magdiff2 = -2.5 * alog10(spectroflux * groupratio[groupnum] / photoflux)

   ;----------
   ; Derive the reduced chi^2 for each plate

   rchi2plate = fltarr(nplate)
   for iplate=0, nplate-1 do begin
      indx = where(zall.plate EQ plist[iplate].plate $
       AND zall.mjd EQ plist[iplate].mjd, nthis)
      if (nthis GT 0) then begin
         ngroup1 = n_elements(uniq(groupnum[indx]))
         thisdof = nthis - ngroup1
         rchi2plate[iplate] = total(chivec[indx]^2) / thisdof
      endif
   endfor

   ;----------
   ; Compute the magnitude shift that we needed to apply to each
   ; group of spectroscopic mags to agree with the photo mags

   magoffset = -2.5 * alog10(groupratio)

   ;----------
   ; Make plots

   datestring = strlowcase(string((strsplit(systime(),/extract))[[2,1,4]], $
    format='(i2.2,a,a)'))
   plottitle = 'Best-Fit ' + filtname[filternum]+'-band Filter ' + datestring
   plotfile = 'sdss_djs_' + datestring + '_' + filtname[filternum] + '.ps'

   csize = 1.1
   dfpsplot, plotfile, /square, /color

   xrange = [ bigwave[(where(fbest GT 0.01))[0]] - 300, $
    bigwave[(reverse(where(fbest GT 0.01)))[0]] + 300 ]
   plot, bigwave, gunnfilt[*,filternum]>fguess>fbest, /nodata, $
    xtitle='Air Wavelength [Ang]', ytitle='Filter Response at 1.3 Airmass', $
    charsize=csize, xrange=xrange, /xstyle, title=plottitle
   djs_oplot, bigwave, gunnfilt[*,filternum], color='cyan'
   djs_oplot, bigwave, fguess, color='red'
   djs_oplot, bigwave, fbest, color='green'
   xplot = total(!x.crange * [0.95,0.05])
   yplot = !y.crange[1]
   djs_xyouts, xplot, 0.32*yplot, 'Mamoru/Gunn Jun-2001', $
    charsize=csize, color='cyan'
   djs_xyouts, xplot, 0.26*yplot, 'Initial guess for fit' $
     + ' \chi^2_r=' + string(origrchi2,format='(f6.3)'), $
    charsize=csize, color='red'
   djs_xyouts, xplot, 0.20*yplot, 'Schlegel Best-Fit '+datestring $
    + ' \chi^2_r=' + string(chi2pdof,format='(f6.3)'), $
    charsize=csize, color='green'
   djs_xyouts, xplot, 0.14*yplot, 'DOF =' + string(dof) $
    + '  \Delta \chi^2 =' + string((origrchi2-chi2pdof)*dof), $
    charsize=csize

   !p.multi = [0,1,2]
   iplot = where(rchi2plate GT 0)
   djs_plot, [plist[iplot].plate], [rchi2plate[iplot]], psym=4, $
    xtitle='Plate Number', ytitle='\chi^2 / DOF', $
    charsize=csize, title=plottitle

   ; These are the values that we would *subtract* from the spectro mags
   iuniq = uniq(groupnum)
   plot, [zall[iuniq].plate], [magoffset], psym=4, $
    xtitle='Plate Number', ytitle='(SPECTRO - PHOTO) Mag Offset per group', $
    charsize=csize, title=plottitle
   !p.multi = 0

   !p.multi = [0,1,2]
   plot, zall.plate + zall.fiberid/1000., guessdiff1, psym=3, $
    xtitle='Plate Number', ytitle='(SPECTRO - PHOTO) w/out mag offsets', $
    charsize=csize, title='Initial Guess Filter'
   oplot, !x.crange, [0,0]
   plot, zall.plate + zall.fiberid/1000., magdiff1, psym=3, $
    yrange=!y.crange, /ystyle, $ ; Use same Y plotting limits as above
    xtitle='Plate Number', ytitle='(SPECTRO - PHOTO) w/out mag offsets', $
    charsize=csize, title=plottitle
   oplot, !x.crange, [0,0]

   plot, zall.plate + zall.fiberid/1000., guessdiff2, psym=3, $
    xtitle='Plate Number', ytitle='(SPECTRO - PHOTO) w/ mag offsets', $
    charsize=csize, title='Initial Guess Filter'
   oplot, !x.crange, [0,0]
   plot, zall.plate + zall.fiberid/1000., magdiff2, psym=3, $
    yrange=!y.crange, /ystyle, $ ; Use same Y plotting limits as above
    xtitle='Plate Number', ytitle='(SPECTRO - PHOTO) w/ mag offsets', $
    charsize=csize, title=plottitle
   oplot, !x.crange, [0,0]
   ; Label the NWORST worst points, according to their chi-deviation
   djs_xyouts, total(!x.crange*[0.95,0.05]), total(!y.crange*[0.08,0.92]), $
    'Worst outliers by \chi^2 in red', color='red'
   djs_oplot, zall[iworst].plate+zall[iworst].fiberid/1000., $
    magdiff[iworst], psym=4, color='red'
   for i=0, nworst-1 do $
    djs_xyouts, zall[iworst[i]].plate+zall[iworst[i]].fiberid/1000., $
     magdiff2[iworst[i]], string(zall[iworst[i]].plate, $
     zall[iworst[i]].mjd, zall[iworst[i]].fiberid, $
     format='(" ",i4,"/",i5,"-",i3," ")'), orient=90, $
     align=(magdiff2[iworst[i]] LT 0), color='red'

   istar = where(strmatch(zall.class,'STAR*'), nstar)
   iqso = where(strmatch(zall.class,'QSO*'), nqso)
;   photocolor = tsall.psfcounts[2] - tsall.psfcounts[3]
   photocolor = -2.5 * alog10(tsall.psfflux[2] / tsall.psfflux[3])

   if (nstar GT 1) then begin
      plot, photocolor[istar], guessdiff2[istar], psym=3, charsize=csize, $
       xtitle='(r-i) for stars', ytitle='(SPECTRO - PHOTO) w/ mag offsets', $
       title='Initial Guess Filter'
      oplot, !x.crange, [0,0]
      isort = sort(photocolor[istar])
      djs_oplot, photocolor[istar[isort]], $
       djs_median(guessdiff2[istar[isort]],width=51<nstar,boundary='reflect'), $
       color='red'
      djs_xyouts, total(!x.crange*[0.95,0.05]), total(!y.crange*[0.08,0.92]), $
       'RMS = ' + string(stddev(guessdiff2[istar])) + ' mag', $
       color='red', charsize=csize
      djs_xyouts, total(!x.crange*[0.95,0.05]), total(!y.crange*[0.14,0.86]), $
       'Running median of 51 pts', color='red', charsize=csize
      plot, photocolor[istar], magdiff2[istar], psym=3, charsize=csize, $
       xtitle='(r-i) for stars', ytitle='(SPECTRO - PHOTO) w/ mag offsets', $
       title=plottitle
      oplot, !x.crange, [0,0]
      djs_oplot, photocolor[istar[isort]], $
       djs_median(magdiff2[istar[isort]],width=51<nstar,boundary='reflect'), $
       color='green'
      djs_xyouts, total(!x.crange*[0.95,0.05]), total(!y.crange*[0.08,0.92]), $
       'RMS = ' + string(stddev(magdiff2[istar])) + ' mag', $
       color='green', charsize=csize
      djs_xyouts, total(!x.crange*[0.95,0.05]), total(!y.crange*[0.14,0.86]), $
       'Running median of 51 pts', color='green', charsize=csize
   endif

   if (nqso GT 1) then begin
      plot, zall[iqso].z, guessdiff2[iqso], psym=3, charsize=csize, $
       xtitle='z for QSOs', ytitle='(SPECTRO - PHOTO) w/ mag offsets', $
       title='Initial Guess Filter'
      oplot, !x.crange, [0,0]
      isort = sort(zall[iqso].z)
      djs_oplot, zall[iqso[isort]].z, $
       djs_median(guessdiff2[iqso[isort]],width=51<nqso,boundary='reflect'), $
       color='red'
      djs_xyouts, total(!x.crange*[0.95,0.05]), total(!y.crange*[0.08,0.92]), $
       'RMS = ' + string(stddev(guessdiff2[iqso])) + ' mag', $
       color='red', charsize=csize
      djs_xyouts, total(!x.crange*[0.95,0.05]), total(!y.crange*[0.14,0.86]), $
       'Running median of 51 pts', color='red', charsize=csize
      plot, zall[iqso].z, magdiff2[iqso], psym=3, charsize=csize, $
       xtitle='z for QSOs', ytitle='(SPECTRO - PHOTO) w/ mag offsets', $
       title=plottitle
      oplot, !x.crange, [0,0]
      djs_oplot, zall[iqso[isort]].z, $
       djs_median(magdiff2[iqso[isort]],width=51<nqso,boundary='reflect'), $
       color='green'
      djs_xyouts, total(!x.crange*[0.95,0.05]), total(!y.crange*[0.08,0.92]), $
       'RMS = ' + string(stddev(magdiff2[iqso])) + ' mag', $
        color='green', charsize=csize
      djs_xyouts, total(!x.crange*[0.95,0.05]), total(!y.crange*[0.14,0.86]), $
       'Running median of 51 pts', color='green', charsize=csize
      !p.multi = 0
   endif

   dfpsclose

   ;----------
   ; Write the filter curve to a file

   outfile = 'sdss_djs_' + datestring + '_' + filtname[filternum] + '.dat'
   openw, olun, outfile, /get_lun
   printf, olun, '# Filename = ' + outfile
   printf, olun, '# Generated by SOLVEFILTER in idlspec2d ' + idlspec2d_version()
   printf, olun, '# Generated on ' + systime()
   printf, olun, '#'
   printf, olun, '# FILTER = ' + filtname[filternum]
   printf, olun, '# FILTTYPE = ' + filttype
   printf, olun, '# STARERR = ' + string(starerr)
   printf, olun, '# QSOERR = ' + string(qsoerr)
   printf, olun, '# SNCUT = ' + string(sncut)
   printf, olun, '# WAVEMIN = ' + string(wavemin)
   printf, olun, '# WAVEMAX = ' + string(wavemax)
   printf, olun, '# MAXITER = ' + string(maxiter)
   printf, olun, '# BOSS_SPECTRO_REDUX = ' + getenv('BOSS_SPECTRO_REDUX')
   if (keyword_set(fluxpath)) then printf, olun, '# FLUXPATH = ' + fluxpath
   printf, olun, '#'
   printf, olun, '# Number of plates = ' + string(nplate)
   printf, olun, '# Number of spectra = ' + string(ntot)
   printf, olun, '# Number of groups = ' + string(max(groupnum)+1)
   printf, olun, '# Start values = ' + string(parinfo.value, format='(99e14.6)')
   printf, olun, '# Fit values = ' + string(theta, format='(99e14.6)')
   printf, olun, '# Fit errors = ' + string(perror, format='(99e14.6)')
   printf, olun, '# Number of iterations = ' + string(niter)
   printf, olun, '#'
   printf, olun, '# Median |chi| = ' + string(median(abs(chivec)))
   printf, olun, '# Start Chi2/DOF = ', origrchi2
   printf, olun, '# Final Chi2/DOF = ' + string(chi2pdof)
   printf, olun, '#'
   printf, olun, '# lambda  response response resnoa   xatm(1.3)'
   for i=0, nbigpix-1 do $
    printf, olun, bigwave[i], fbest[i], fbest[i], $
     fbest[i] * exp(1.3 * tauextinct[i]), exp(-1.3 * tauextinct[i]), $
     format='(f8.2,4f9.5)'
   close, olun
   free_lun, olun

   return
end
;------------------------------------------------------------------------------
