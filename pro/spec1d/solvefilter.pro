
; What are the outlier objects?
; Do Christy's flux-calibrations improve things?
; Which particular objects argue for a blue light leak (using chi^2)?

forward_function mpfit, mpfitfun, solvefiltfn

;------------------------------------------------------------------------------
; Construct the filter curve corresponding to this set of parameters
function solvefiltshape, theta, loglam

   mm = (theta[4] - 1.d0) / (theta[1] - theta[0])
   bb = 1.d0 - mm * theta[0]
   fcurve = (tanh((loglam - theta[0])*theta[2]) + 1.d0) $
          * (tanh((theta[1] - loglam)*theta[3]) + 1.d0) $
          * (mm * loglam + bb)

   if (total(finite(fcurve)) NE n_elements(fcurve)) then $
    message, 'NaN in filter shape'

   return, fcurve
end
;------------------------------------------------------------------------------
function solvefiltfn, theta

   common com_solvefilt, groupnum, bigloglam, $
    bigflux, photoflux, photoinvsig, tauextinct, ntot, nbigpix, $
    airmass

   ; Construct the filter curve corresponding to this set of parameters
   fcurve = solvefiltshape(theta, bigloglam)

   ; Change from f_lambda to f_nu
   flambda2fnu = 10^(2*bigloglam) / 2.99792d18

   ; Loop through each group of spectra, integrate over the filter curve,
   ; and minimize the spectro/photo flux normalization for that group.
   leftall = dblarr(ntot)
   rightall = dblarr(ntot)
   for igroup=0L, max(groupnum) do begin
      indx = where(groupnum EQ igroup, nthis)

      ; Get the extinction curve for these objects
      meanair = mean(airmass[indx])
      fextinct = exp(-meanair * tauextinct)
      fmult = fcurve * flambda2fnu * fextinct

      leftval = total(bigflux[*,indx] * rebin(fmult,nbigpix,nthis), 1) $
       * photoinvsig[indx]
      rightval = photoflux[indx] * photoinvsig[indx]
      denom = total(leftval^2)
      if (denom GT 0) then ratio = total(leftval * rightval) / denom $
       else ratio = 1
      leftval = leftval * ratio
      leftall[indx] = leftval
      rightall[indx] = rightval
   endfor

   ; Return a vector of all the chi's.
   return, leftall - rightall
end
;------------------------------------------------------------------------------
pro solvefilter

filternum = 3 ; ???
;path = '/scr/wire50/cat/recalib/kurucz' ; ???

   if (NOT keyword_set(adderr)) then adderr = 0.06 ; Fractional error to add
   if (NOT keyword_set(wavemin)) then wavemin = 3800.d0
   if (NOT keyword_set(wavemax)) then wavemax = 9300.d0
   filtname = ['u','g','r','i','z']

   ;----------
   ; Set the parameters used for the initial guesses

   wave1 = [0, 0, 0, 6902., 0]
   wave2 = [0, 0, 0, 8194., 0]

   t0 = systime(1)

   ;----------
   ; Set common block

   common com_solvefilt, groupnum, bigloglam, $
    bigflux, photoflux, photoinvsig, tauextinct, ntot, nbigpix, $
    airmass

   bigloglam = wavevector(alog10(wavemin), alog10(wavemax))

   ;----------
   ; Append Gunn's measurement of the filter -- including the telluric bands.

   filename = filepath('sdss_jun2001_'+filtname[filternum]+'_atm.dat', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, filename, fwave, fthru1, fthru2, fthru3, fext, /silent
   fthru = 0.5 * (fthru1 + fthru2) ; Average these two columns
   linterp, alog10(fwave), fthru, bigloglam, gunnfilt

   ;----------
   ; Construct the atmospheric extinction curve

   ; Start with a simple expression for the extinction
   tausimple = 10^(12.3 - 3.6 * bigloglam)

   framefile = filepath('spFrame-r2-00007466.fits.gz', $
    root_dir=getenv('SPECTRO_DATA'), subdir='0432')
   hdr = headfits(framefile)
   thisair = sxpar(hdr, 'AIRMASS')
   wset = mrdfits(framefile, 3)
   traceset2xy, wset, xx, tloglam
   telluric = mrdfits(framefile, 8)
   linterp, tloglam[*,0], telluric[*,0], bigloglam, tellcorr
   tautelluric = -alog(tellcorr/thisair) ; Scale back to one airmass

   tauextinct = (tausimple + tautelluric) > 0

;   linterp, alog10(fwave), fext, bigloglam, gunnextinct
;   taugunn = -alog(gunnextinct) / 1.2 ; at 1.2 airmasses

   ;----------
   ; Find the list of good plates

   splog, 'Find list of good plates'
   platelist, plist=plist
   i = where(strmatch(plist.platequality,'good*') $
    AND strmatch(plist.smearuse,'T*') $
    AND plist.plate GT 431 AND plist.mjd LT 52220)
;i = where(strmatch(plist.platequality,'good*') $
; AND (plist.plate EQ 323))
i = where(strmatch(plist.public,'*DR1*') AND plist.plate GT 431)
;i = where(plist.plate EQ 323 OR plist.plate EQ 324)
;i = i[0:2] ; ???
   plist = plist[i]

   ;----------
   ; Make certain that the data exists for all these plates,
   ; especially if PATH is set to read the flux vectors from
   ; some non-standard directory.

   splog, 'Test existence of flux data for all plates'
   nplate = n_elements(plist)
   qkeep = bytarr(nplate)
   for iplate=0, nplate-1 do begin
      readspec, plist[iplate].plate, mjd=plist[iplate].mjd, path=path, $
       objhdr=objhdr
      if (keyword_set(objhdr)) then qkeep[iplate] = 1B
   endfor
   ikeep = where(qkeep, nplate)
   plist = plist[ikeep]
   splog, 'Number of usable plates = ', nplate

   ;----------
   ; Loop over each plate, and accumulate groups of objects

   nplate = n_elements(plist)
   gcounter = 0L
   for iplate=0, nplate-1 do begin
      splog, 'Reading redshift data for plate #', plist[iplate].plate, $
       ' (', iplate, ' of ', nplate, ')'
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
          AND zans[indx].sn_median GT 3, ngood)

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

   nbigpix = n_elements(bigloglam)
   bigflux = fltarr(nbigpix, ntot)

   ;----------
   ; Read in the actual spectra

   splog, 'Reading spectra'
   for iplate=0, nplate-1 do begin
      splog, 'Reading spectra for plate #', plist[iplate].plate, $
       ' (', iplate, ' of ', nplate, ')'
      indx = where(zall.plate EQ plist[iplate].plate $
       AND zall.mjd EQ plist[iplate].mjd, nthis)
      if (nthis GT 0) then begin
         ; Read in the spectra, and interpolate over bad points
         fiberid = zall[indx].fiberid
         readspec, plist[iplate].plate, mjd=plist[iplate].mjd, fiberid, $
          flux=objflux, loglam=loglam, invvar=objivar, $
          andmask=andmask, ormask=ormask, path=path, /align
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

   ;----------
   ; Decide upon the object counts and errors from PHOTO.

   photoflux = 10.d0^(-tsall.psfcounts[filternum]/2.5)
   photoflerr = tsall.psfcountserr[filternum] * abs(photoflux)
   photoinvsig = 1. / sqrt( photoflerr^2 + (adderr*abs(photoflux))^2 )

   ;----------
   ; Set the structure to pass to MPFIT

   parinfo = replicate({value:0.D, fixed:0, limited:[0b,0b], $
    limits:[0.D,0], mpmaxstep: 0.D}, 5)
maxiter = 100
   loglam1 = alog10(wave1[filternum]+0.d0)
   loglam2 = alog10(wave2[filternum]+0.d0)
   logmid = 0.5 * (loglam1 + loglam2)
   parinfo.value = [loglam1, loglam2, 200., 200., 0.55]
   parinfo.limited = [[0,1], [1,0], [1,0], [1,0], [1,0]]
   parinfo.limits = [[0,logmid-5.d-4], [logmid+5.d-4,0], $
    [0.0,0], [0.0,0], [0,0]]
   parinfo.mpmaxstep = [5.d-4, 5.d-4, 10., 10., 0.05]

   t1 = systime(1)
   theta = mpfit('solvefiltfn', $
    parinfo=parinfo, perror=perror, maxiter=maxiter, $
    nfev=nfev, niter=niter, status=status)
   splog, 'Time for non-linear fitting = ', systime(1)-t1, ' sec'
   splog, 'Number of iterations = ', niter
stop
   chivec = solvefiltfn(theta)
print, perror
print, total(chivec^2) / (n_elements(chivec) - n_elements(theta))
   fbest = solvefiltshape(theta, bigloglam) * exp(-1.2 * tauextinct)
   fguess = solvefiltshape(parinfo.value, bigloglam) * exp(-1.2 * tauextinct)
splot, 10^bigloglam, gunnfilt/max(gunnfilt)
soplot, 10^bigloglam, fbest/max(fbest), color='green'
soplot, 10^bigloglam, fguess/max(fguess), color='red'

stop
   return
end
;------------------------------------------------------------------------------
