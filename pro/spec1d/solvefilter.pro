
; COMMENTS:
;   We only use unblended stars and QSOs.
;   All calculations are done in vacuum wavelengths, but then converted
;     to air wavelengths at the end.
;
; EXAMPLES:
;   Solve for the i-band filter using Tremonti's re-reductions of the spectra:
;     IDL> solvefilter, filtnum=3, fluxpath='/scr/wire50/cat/recalib/kurucz'

; What are the outlier objects?
; Plot the residuals as a fn of plate number
; Do Christy's flux-calibrations improve things?
; Which particular objects argue for a blue light leak (using chi^2)?
; Should we more heavily weight the QSOs?
; I should average together more telluric spectra for better S/N.
; Use the extinction coeff for each imaging night.

forward_function mpfit, mpfitfun, solvefiltfn

;------------------------------------------------------------------------------
; Construct the filter curve corresponding to this set of parameters
function solvefiltshape, theta, loglam

   common com_solvefilt, groupnum, bigloglam, $
    bigflux, photoflux, photoinvsig, tauextinct, ntot, nbigpix, $
    airmass, gunnfilt, filternum, filttype, groupratio, spectroflux

   case filttype of
   'tanh': begin
      mm = (theta[4] - 1.d0) / (theta[1] - theta[0])
      bb = 1.d0 - mm * theta[0]
      fcurve = (tanh((loglam - theta[0])*theta[2]) + 1.d0) $
             * (tanh((theta[1] - loglam)*theta[3]) + 1.d0) $
             * (mm * loglam + bb)
      end
   'sdss': begin
      slopeterm = max(bigloglam) + theta[2] * bigloglam
      linterp, bigloglam, gunnfilt[*,filternum] * slopeterm, $
       (bigloglam - theta[0]) * theta[1], fcurve
      end
   endcase

   if (total(finite(fcurve)) NE n_elements(fcurve)) then $
    message, 'NaN in filter shape'

   return, fcurve
end
;------------------------------------------------------------------------------
function solvefiltfn, theta

   common com_solvefilt, groupnum, bigloglam, $
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

      ; Get the extinction curve for these objects
      meanair = mean(airmass[indx])
      fextinct = exp(-meanair * tauextinct)
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
 plate=plate, mjd=mjd, adderr=adderr, wavemin=wavemin, wavemax=wavemax, $
 sncut=sncut, maxiter=maxiter, fluxpath=fluxpath

   ;----------
   ; Set common block

   common com_solvefilt, groupnum, bigloglam, $
    bigflux, photoflux, photoinvsig, tauextinct, ntot, nbigpix, $
    airmass, gunnfilt, filternum, filttype, groupratio, spectroflux

   if (keyword_set(filttype1)) then filttype = filttype1 $
    else filttype = 'sdss'
   if (n_elements(filternum1) NE 0) then filternum = filternum1 $
    else filternum = 3 ; Default to i-band
   if (NOT keyword_set(adderr)) then adderr = 0.05 ; Fractional error to add
   if (NOT keyword_set(wavemin)) then wavemin = 3800.d0
   if (NOT keyword_set(wavemax)) then wavemax = 9300.d0
   if (NOT keyword_set(sncut)) then sncut = 2.0
   if (NOT keyword_set(maxiter)) then maxiter = 200
   filtname = ['u','g','r','i','z']

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
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
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

   case filttype of
   'tanh': begin
      parinfo = replicate({value:0.D, fixed:0, limited:[0b,0b], $
       limits:[0.D,0], mpmaxstep: 0.D}, 5)
      if (filternum EQ 1) then $
       parinfo.value = [alog10(3950), alog10(5325), 2.0, 100, 140]
      if (filternum EQ 2) then $
       parinfo.value = [alog10(5580), alog10(6750), 1.2, 130, 160]
      if (filternum EQ 3) then $
       parinfo.value = [alog10(6915), alog10(8210), 0.55, 150, 220]
      parinfo.limited = [[0,1], [1,0], [1,0], [1,0], [1,0]]
      logmid = 0.5 * (parinfo[0].value + parinfo[1].value)
      parinfo.limits = [[0,logmid-5.d-4], [logmid+5.d-4,0], $
       [0.0,0], [0.0,0], [0,0]]
      parinfo.mpmaxstep = [5.d-4, 5.d-4, 10., 10., 0.05]
      end
   'sdss': begin
      parinfo = replicate({value:0.D, fixed:0, limited:[0b,0b], $
       limits:[0.D,0], mpmaxstep: 0.D}, 3)
      parinfo.value = [-10.d-5, 1.001, 0.01]
      parinfo.mpmaxstep = [5.d-4, 0.01, 0.02]
      end
   else: message, 'Unknown FILTTYPE'
   endcase

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
    root_dir=getenv('SPECTRO_DATA'), subdir='0432')
   framefile = (findfile(framefile))[0]
   hdr = headfits(framefile)
   thisair = sxpar(hdr, 'AIRMASS')
   wset = mrdfits(framefile, 3)
   traceset2xy, wset, xx, tloglam
   telluric = mrdfits(framefile, 8)
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
      iuse = where(strmatch(plist.public,'*DR1*') AND plist.plate GT 431)
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
      qblend = transpose(bflag[3,*] EQ 1 AND bflag[6,*] EQ 0)
      qbright = transpose(bflag[1,*])
      qchild = transpose(bflag[4,*])
      qsingle = (qblend EQ 0) AND (qbright EQ 0) AND (qchild EQ 0)

      ; Compute the airmass for each object
      junk = sdss_run2mu(tsobj.run, tsobj.field, tai=tai)
      airmass1 = tai2airmass(zans.plug_ra, zans.plug_dec, tai=tai)
      if (min(airmass1) LT 0.99 OR max(airmass1) GT 3) then $
       message, 'Invalid AIRMASS'

      ; Group objects with the same run+rerun+camcol+plate+spectrographid
      idstring = string(tsobj.run) + string(tsobj.rerun) $
       + string(tsobj.camcol) + string(zans.plate) + string(plug.spectrographid)
      idlist = idstring[ uniq(idstring, sort(idstring)) ]
      ngroup = n_elements(idlist)

      for igroup=0, ngroup-1 do begin
         indx = where(idstring EQ idlist[igroup])

         ; Reject wild mag outliers, which are often objects where there
         ; is a  bright blend but the PHOTO flux is only for a fainter child
         magdiff = tsobj[indx].psfcounts[filternum] $
          + 2.5 * alog10(zans[indx].counts_spectro[filternum])

         igood = where(zans[indx].zwarning EQ 0 $
          AND qsingle[indx] EQ 1 $
          AND qgalaxy[indx] EQ 0 $
          AND ( strmatch(zans[indx].class,'STAR*') $
           OR strmatch(zans[indx].class,'QSO*') ) $
          AND magdiff GT median([magdiff])-0.5 $
          AND magdiff LT median([magdiff])+0.5 $
          AND zans[indx].wcoverage GT 0.30 $
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

   photoflux = 10.d0^(-tsall.psfcounts[filternum]/2.5)
   photoflerr = tsall.psfcountserr[filternum] * abs(photoflux)
   photoinvsig = 1. / sqrt( photoflerr^2 + (adderr*abs(photoflux))^2 )

   ; Now convert these to AB flux, according to the numbers derived
   ; by Hogg on 13 Aug 2002.
   aboffsets = [-0.042, 0.036, 0.015, 0.013, -0.002]
   photoflux = photoflux * 10.d0^(-aboffsets[filternum]/2.5)

   ;----------
   ; Do the actual fit to the filter curve

   t1 = systime(1)
   theta = mpfit('solvefiltfn', $
    parinfo=parinfo, perror=perror, maxiter=maxiter, $
    nfev=nfev, niter=niter, status=status)
   chivec = solvefiltfn(theta)
   dof = n_elements(chivec) - n_elements(theta) - ngroup
   chi2pdof = total(chivec^2) / dof

   splog, 'Time for non-linear fitting = ', systime(1)-t1, ' sec'
   splog, 'Number of iterations = ', niter
   splog, 'Number of function evaluations = ', nfev
   splog, 'Fit values = ', theta
   splog, 'Fit errors = ', perror
   splog
   splog, 'Median |chi| = ', median(abs(chivec))
   splog, 'Chi2/DOF = ', chi2pdof

   ;----------
   ; Identify the 10 most deviant points

   nworst = 10
   iworst = (reverse(sort(abs(chivec))))[0:nworst-1]

   ;----------
   ; Reconstruct the filters at 1.3 airmasses

   fguess = solvefiltshape(parinfo.value, bigloglam) * exp(-1.3 * tauextinct)
   fguess = fguess * mean(gunnfilt[*,filternum]) / mean(fguess)

   fbest = solvefiltshape(theta, bigloglam) * exp(-1.3 * tauextinct)
   fbest = fbest * mean(gunnfilt[*,filternum]) / mean(fbest)
   fbest = fbest * (fbest GT 0) + 0.0 * (fbest LE 0) ; Get rid of values -0.00

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

   csize = 2
   dfpsplot, plotfile, /square, /color

   djs_plot, plist.plate, rchi2plate, psym=4, $
    xtitle='Plate Number', ytitle='\chi^2 / DOF', $
    charsize=csize, title=plottitle

   ; These are the values that we would *subtract* from the spectro mags
   iuniq = uniq(groupnum)
   plot, zall[iuniq].plate, magoffset, psym=4, $
    xtitle='Plate Number', ytitle='(SPECTRO - PHOTO) Mag Offset per group', $
    charsize=csize, title=plottitle

   magdiff = -2.5 * alog10(spectroflux / photoflux)
   plot, zall.plate + zall.fiberid/1000., magdiff, psym=3, $
    xtitle='Plate Number', ytitle='(SPECTRO - PHOTO) w/out mag offsets', $
    charsize=csize, title=plottitle
   magdiff = -2.5 * alog10(spectroflux * groupratio[groupnum] / photoflux)
   plot, zall.plate + zall.fiberid/1000., magdiff, psym=3, $
    xtitle='Plate Number', ytitle='(SPECTRO - PHOTO) w/ mag offsets', $
    charsize=csize, title=plottitle
   ; Label the NWORST worst points, according to their chi-deviation
   djs_oplot, zall[iworst].plate+zall[iworst].fiberid/1000., $
    magdiff[iworst], psym=4, color='red'
   for i=0, nworst-1 do $
    djs_xyouts, zall[iworst[i]].plate+zall[iworst[i]].fiberid/1000., $
     magdiff[iworst[i]], string(zall[iworst[i]].plate, $
     zall[iworst[i]].mjd, zall[iworst[i]].fiberid, $
     format='(i4,"/",i5,"-",i3," ")'), orient=90, align=1.0, color='red'

   xrange = [ bigwave[(where(fbest GT 0.02))[0]] - 300, $
    bigwave[(reverse(where(fbest GT 0.02)))[0]] + 300 ]
   plot, bigwave, gunnfilt[*,filternum], $
    xtitle='Air Wavelength [Ang]', ytitle='Filter Response at 1.3 Airmass', $
    charsize=csize, xrange=xrange, /xstyle, title=plottitle
;   djs_oplot, bigwave, fguess, color='red'
   djs_oplot, bigwave, fbest, color='green'
   xplot = total(!x.crange * [0.95,0.05])
   yplot = !y.crange[1]
   djs_xyouts, xplot, 0.92*yplot, 'Gunn Jun-2001 Curve', charsize=0.75*csize
   djs_xyouts, xplot, 0.86*yplot, 'Schlegel Best-Fit '+datestring, $
   charsize=0.75*csize, color='green'

   dfpsclose

   ;----------
   ; Write the filter curve to a file

   outfile = 'sdss_djs_' + datestring + '_' + filtname[filternum] + '.dat'
   openw, olun, outfile, /get_lun
   printf, olun, '# Filename = ' + outfile
   printf, olun, '# Generated on ' + systime()
   printf, olun, '#'
   printf, olun, '# Number of plates = ' + string(nplate)
   printf, olun, '# Number of spectra = ' + string(ntot)
   printf, olun, '# Number of groups = ' + string(max(groupnum)+1)
   printf, olun, '# Fit values = ' + string(theta, format='(99f)')
   printf, olun, '# Fit errors = ' + string(perror, format='(99f)')
   printf, olun, '# Median |chi| = ' + string(median(abs(chivec)))
   printf, olun, '# Chi2/DOF = ' + string(chi2pdof)
   if (keyword_set(fluxpath)) then print, olun, '# FLUXPATH = ' + fluxpath
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
