;+
; NAME:
;   star_dvelocity
;
; PURPOSE:
;   Solve for velocity shifts between multiple exposures.
;
; CALLING SEQUENCE:
;   res = star_dvelocity(plate, mjd=, [ vmin=, vmax=, fiberid=, path= ] )
;
; INPUTS:
;   plate      - Plate number
;   mjd        - MJD for this plate
;
; OPTIONAL INPUTS:
;   vmin       - Minimum velocity to consider; default to -700.
;   vmax       - Maxmimum velocity to consider; default to 700.
;   fiberid    - If set, then only fit for the objects specified by
;                these 1-indexed fiber IDs.  Set to 0 to not
;                fit any objects (but return empty 640-element output
;                structures); default to fitting all.
;   path       - Override all path information with this directory name
;                (for reading the spCFrame files)
;
; OUTPUTS:
;   res        - Output structure with result for each object [NOBJECT].
;                The structure contains the following
;                best-fit Elodie spectrum:
;                  VELOCITY_TIME[10] - Timestamp for each exposure
;                  VELOCITY[10]      - Velocity for each exposure [km/s]
;                  VELOCITY_ERR[10]  - Error in velocity
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Velocities are computed relative to the first exposure.  Each exposure
;   is compared to every other exposure, with no assumption about the kind
;   of object.  The object need not be a star.
;
; EXAMPLES:
;   Fit the velocity shifts for the white dwarf star in SDSS publication #535:
;     IDL> res=star_dvelocity(518,mjd=52282,fiberid=285)
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;   dvelocity_struct()
;   dvelocity_fn()
;
; REVISION HISTORY:
;   09-Jul-2005  Written by D. Schlegel, LBL
;------------------------------------------------------------------------------
function dvelocity_struct

   nmax = 20

   result = create_struct( $
    name = 'VSHIFT', $
    'velocity_tai', dblarr(nmax), $
    'velocity'    , fltarr(nmax), $
    'velocity_err', fltarr(nmax), $
    'velocity_covar', fltarr(nmax,nmax) $
   )

   return, result
end
;------------------------------------------------------------------------------
function dvelocity_fn, value, loglam=loglam, objflux=objflux, objivar=objivar, $
 ilist=ilist, newloglam=newloglam

   nexp = max(ilist) + 1L
   if (nexp LE 1) then return, 0.

   nnewpix = n_elements(newloglam)
   newflux = fltarr(nnewpix, nexp)
   newivar = fltarr(nnewpix, nexp)

   ; Resample all of the spectra to the same sampling, but shifting
   ; by the (test) velocity of each exposure.
   for iexp=0L, nexp-1L do begin
      ii = where(ilist EQ iexp)
      combine1fiber, loglam[*,ii]-value[iexp], objflux[*,ii], objivar[*,ii], $
       newloglam=newloglam, newflux=newflux1, newivar=newivar1
      newflux[*,iexp] = newflux1
      newivar[*,iexp] = newivar1
   endfor

   ; Compute the chi vector between all possible pairs of exposures
   npair = nexp * (nexp - 1) / 2
   chivec = fltarr(npair * nnewpix)
   ipair = 0L
   for i=0, nexp-2 do begin
      for j=i+1, nexp-1 do begin
         qgood = newivar[*,i] GT 0 AND newivar[*,j] GT 0
         thisivar = qgood / (1./(newivar[*,i] + (qgood EQ 0)) $
          + 1./(newivar[*,j] + (qgood EQ 0)))
         chivec[ipair*nnewpix : (ipair+1)*nnewpix-1] = $
          (newflux[*,i] - newflux[*,j]) * sqrt(thisivar)
         ipair = ipair + 1L
      endfor
   endfor

;   chivec = fltarr((nexp-1) * nnewpix)
;   for i=1, nexp-1 do $
;    chivec[(i-1)*nnewpix : i*nnewpix-1] = (newflux[*,0] - newflux[*,i]) $
;     * sqrt(newivar[*,i])

   return, chivec
end

;------------------------------------------------------------------------------
function star_dvelocity, plate, mjd=mjd, fiberid=fiberid1, path=path

   if (n_params() LT 1) then begin
      print, 'Syntax - res = star_dvelocity(plate, mjd=, $'
      print, ' [ fiberid=, path= ] )'
      return, 0
   endif
   if (NOT keyword_set(vmin)) then vmin = -700.
   if (NOT keyword_set(vmax)) then vmax = 700.

   stime0 = systime(1)

   ;----------
   ; Create the output structure, and return empty data if FITINDX=[-1].

   ntot = 640L
   res_best = replicate(dvelocity_struct(), ntot)

   if (n_elements(fiberid1) EQ 1 AND keyword_set(fiberid1) EQ 0) $
    then return, res_best ; If FIBERID was set to zero
   if (keyword_set(fiberid1)) then fiberid = fiberid1 $
    else fiberid = lindgen(ntot)

   ;----------
   ; Compute the velocities for each object

   nobj = n_elements(fiberid)
   for iobj=0L, nobj-1L do begin
      splog, 'Object number ', iobj, ' of ', nobj, ', Fiber #', fiberid[iobj]

      readonespec, plate, mjd=mjd, fiberid[iobj], $
       framehdr=framehdr, loglam=loglam, flux=objflux, invvar=objivar, $
       mask=mask, path=path
      nfile = n_elements(framehdr)
      objivar = objivar * ((mask AND pixelmask_bits('COMBINEREJ')) EQ 0)
;objivar = skymask(objivar, mask, mask) ; ???

      ; Make the wavelenghts ascending for each spectrum
      for ifile=0L, nfile-1L do begin
         if (loglam[0,ifile] GT loglam[1,ifile]) then begin
            loglam[*,ifile] = reverse(loglam[*,ifile])
            objflux[*,ifile] = reverse(objflux[*,ifile])
            objivar[*,ifile] = reverse(objivar[*,ifile])
            mask[*,ifile] = reverse(mask[*,ifile])
         endif
      endfor

      ;----------
      ; Get the list of exposure numbers

      expnum = lonarr(nfile)
      for ifile=0L, nfile-1L do $
       expnum[ifile] = sxpar(*framehdr[ifile],'EXPOSURE')
      iuniq = uniq(expnum, sort(expnum))
      nexp = n_elements(iuniq)
      for i=0L, nexp-1L do $
       res_best[fiberid[iobj]-1].velocity_tai[i] = $
        sxpar(*framehdr[iuniq[i]],'TAI')

      ; Define ILIST as 0...NEXP-1 for each input file (spectrum)
      explist = expnum[iuniq]
      ilist = lonarr(nfile)
      for ifile=0L, nfile-1L do ilist[ifile] = where(expnum[ifile] EQ explist)

      minlam = min(loglam, max=maxlam)
      newloglam = wavevector(minlam, maxlam)

      if (nexp GT 1) then begin

         ;----------
         ; Compute the relative velocity between each possible pair of exposures

subsamp = 10 ; ???
width = 1.e-4 ; ???
maxshift = 20.d-4 ; ???
minchi2diff = 5. ; ???
         binsz = 1.0d-4/subsamp
         bigloglam = wavevector(minlam - maxshift, maxlam + maxshift, $
          binsz=binsz)
         nbigpix = n_elements(bigloglam)
         ipix1 = long(maxshift/binsz) $
          + lindgen((maxlam-minlam)/1.0d-4) * subsamp
         nlag = 2 * long(maxshift/binsz) + 1
         lags = lindgen(nlag) - (nlag-1)/2

         ; Resample all of the spectra on the subsampled wavelength grid
         bigflux = fltarr(nbigpix, nexp)
         bigivar = fltarr(nbigpix, nexp)
         for iexp=0L, nexp-1L do begin
            ii = where(ilist EQ iexp)
            combine1fiber, loglam[*,ii], objflux[*,ii], objivar[*,ii], $
             newloglam=bigloglam, newflux=bigflux1, newivar=bigivar1
            bigflux[*,iexp] = bigflux1
            bigivar[*,iexp] = bigivar1
         endfor

         npair = nexp * (nexp - 1) / 2
         mmatrix = dblarr(npair+1,nexp)
         mmatrix[0,0] = 1.d0 ; Force velocity of first exposure to zero
         bvec = dblarr(npair+1)
         sqivar = dblarr(npair+1)
         ipair = 1L
         for i1=0, nexp-2 do begin
            for i2=i1+1, nexp-1 do begin
               chivec = dblarr(nlag)
               for ilag=0L, nlag-1L do begin
                  qgood = bigivar[ipix1,i1] GT 0 AND $
                   bigivar[ipix1+lags[ilag],i2] GT 0
                  thisivar = qgood / (1./(bigivar[ipix1,i1] + (qgood EQ 0)) $
                   + 1./(bigivar[ipix1+lags[ilag],i2] + (qgood EQ 0)))
                  chivec[ilag] = total( $
                   (bigflux[ipix1,i1] - bigflux[ipix1+lags[ilag],i2])^2 $
                   * thisivar )
               endfor

               ; Find the best-fit shift in log-wavelength
               xpeak1 = find_nminima(chivec, lags*1.d-4/subsamp, xerr=xerr1, $
                width=width, nfind=2, ypeak=ypeak1) ;, /debug,/doplot)
               if (n_elements(ypeak1) GT 1) then $
                if (ypeak1[1] GT 0 AND ypeak1[1] LT ypeak1[0]+minchi2diff $
                 and abs(xpeak1[0]-xpeak1[1]) LT xerr1[0]) then $
                 xerr1[0] = abs(xpeak1[0]-xpeak1[1])
               mmatrix[ipair,i1] = 1.
               mmatrix[ipair,i2] = -1.
               bvec[ipair] = xpeak1[0]
               if (xerr1[0] GT 0) then sqivar[ipair] = 1. / xerr1[0]
               ipair = ipair + 1
            endfor
         endfor
         sqivar[0] = total(sqivar) ; Heavily weight the datum that the velocity
                                   ; of the first exposure is zero
         junk = computechi2(bvec, sqivar, mmatrix, acoeff=fitvelocity1, var=var)

         ;----------
         ; Do the fit

         parinfo = {value: 0.D, fixed: 0, limited: [0b,0b], $
          limits: [0.d0,0.d0], tied: '', step: 2.e-6, mpmaxstep: 1.e-4 }
         parinfo = replicate(parinfo, nexp)
         parinfo.value = fitvelocity1
; ???
         parinfo[0].fixed = 1
;         for i=1L, nexp-1 do parinfo[0].tied = parinfo[0].tied $
;          + ' - P(' + strtrim(string(i),2) + ')'

         ftol = 1d-20
         gtol = 1d-20
         xtol = 1d-20
         maxiter = 50L
         functargs = { loglam: loglam, objflux: objflux, objivar: objivar, $
          ilist: ilist, newloglam: newloglam }
         fitval = mpfit('dvelocity_fn', parinfo=parinfo, functargs=functargs, $
          perror=perror, maxiter=maxiter, ftol=ftol, gtol=gtol, xtol=xtol, $
          niter=niter, status=status, covar=covar, /quiet)

         pixscale = alog(10.) * 2.99792458e5 ; (km/sec)/(10e-4 log-wave)
         res_best[fiberid[iobj]-1].velocity[0:nexp-1] = fitval * pixscale
         res_best[fiberid[iobj]-1].velocity_err[0:nexp-1] = perror * pixscale
         res_best[fiberid[iobj]-1].velocity_covar[0:nexp-1,0:nexp-1] = $
          covar * pixscale^2
      endif
   endfor

   splog, 'Total time for STAR_DVELOCITY = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'

   ; Free memory
   for ifile=0L, nfile-1L do $
    ptr_free, framehdr[ifile]

   return, res_best
end
;------------------------------------------------------------------------------
