;+
; NAME:
;   nebular_pcasky
;
; PURPOSE:
;   Solve for the PCA components of the night sky emission
;
; CALLING SEQUENCE:
;   nebular_pcasky, [ waverange=, wavefit=, outfile=, niter=, nkeep= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   waverange  - Wavelength range for output; default to [3750.,9250.] Ang
;   wavefit    - Fitting wavelength range; default to [4000.,7700.] Ang
;   outfile    - Output file; default to 'pcasky.fits'
;   niter      - Parameter for PCA_SOLVE(); default to 10
;   nkeep      - Parameter for PCA_SOLVE(); default to 4
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Make use of the instrumental response at each pixel ???
;
;   Since we do not adequately describe the auroral emission, the
;   three [O I] lines at 5577, 6300, 6363 are simply masked from the
;   fit and interpolated over, such that they don't appear at all
;   in the output PCA spectra.
;   For and example of strong auroal activity, see the spectrum 286/51999-422.
;
; PROCEDURES CALLED:
;   combine1fiber
;   mwrfits
;   pca_solve()
;   platelist()
;   readonespec
;   readspec
;   splog
;   sxpar()
;   wavevector()
;
; REVISION HISTORY:
;   05-Dec-2006  Written by A. West & D. Schlegel, Berkeley
;-
;------------------------------------------------------------------------------
pro nebular_pcasky, waverange=wrange1, wavefit=wfit1, $
 outfile=outfile1, niter=niter1, nkeep=nkeep1

   wrange = n_elements(wrange1) EQ 2 ? wrange1 : [3750.,9250.]
   wfit = n_elements(wfit1) EQ 2 ? wfit1 : [4000.,7700.]
   outfile = keyword_set(outfile1) ? outfile1 : 'pcasky.fits'
   niter = n_elements(niter1) GT 0 ? niter1 : 10
   nkeep = n_elements(nkeep1) GT 0 ? nkeep1 : 4

   newloglam = wavevector(alog10(wrange[0]), alog10(wrange[1]))

   ;----------
   ; Select plates with little nebular emission for computing sky PCA
   ; This region selected based upon Finkbeiner's H-alpha maps.

   platelist, plist=plist
   euler, plist.ra, plist.dec, ll, bb, 1
   ikeep = where(strmatch(plist.status1d,'Done*') $
    AND strmatch(plist.platequality,'Bad*') EQ 0 $
    AND ll GT 270 AND ll LT 290 AND bb GT 50, nplate)
   plist = plist[ikeep]
   splog, 'Number of plates = ', nplate

   ;----------
   ; Determine which are sky fibers on these plates

   splog, 'Reading ZANS structures'
   readspec, plist.plate, mjd=plist.mjd, zans=zans, /silent
   qsky = strmatch(zans.objtype, 'SKY*')
   isky = where(qsky, nsky)
   zans = zans[isky]

   ;----------
   ; Count the total number of sky spectra, one per exposure.
   ; (We will be combining the blue+red spectra below.)

   ntot = 0L
   for iplate=0L, nplate-1L do begin
      print, format='("Counting plate ",i4," of ",i4,a,$)', iplate, nplate, $
       string(13b)
      ii = where(zans.plate EQ plist[iplate].plate $
       AND zans.mjd EQ plist[iplate].mjd AND qsky, nn)
      if (nn GT 0) then begin
         readspec, plist[iplate].plate, mjd=plist[iplate].mjd, $
          objhdr=objhdr, /silent
         expid = long(strmid(sxpar(objhdr, 'EXPID*'),3,8))
         nexp = n_elements(uniq(expid, sort(expid)))
         ntot += nexp * nn
      endif
   endfor
   print
   splog, 'Total number of sky spectra = ', ntot

   ;----------
   ; Read in all the individual exposures for these sky spectra

   i = 0L
   newflux = fltarr(n_elements(newloglam),ntot)
   newivar = fltarr(n_elements(newloglam),ntot)
   for iplate=0L, nplate-1L do begin
      print, format='("Reading plate ",i4," of ",i4,a,$)', iplate, nplate, $
       string(13b)
      ii = where(zans.plate EQ plist[iplate].plate $
       AND zans.mjd EQ plist[iplate].mjd AND qsky, nn)
      for j=0L, nn-1L do begin
         readonespec, zans[ii[j]].plate, mjd=zans[ii[j]].mjd, $
          zans[ii[j]].fiberid, $
          flux=flux1, invvar=invvar1, loglam=loglam1, $
          sky=sky1, framehdr=framehdr, /silent

         ; Revert the heliocentric correction, and put us back in the
         ; rest frame of planet Earth.
         for iexp=0, n_elements(framehdr)-1 do begin
            heliov = sxpar(*framehdr[iexp], 'HELIO_RV')
            loglam1[*,iexp] += alog10(1.d0 + heliov/2.99792458d5)
         endfor

         ; Combine the blue+red spectra for each exposure
         expid = lonarr(n_elements(framehdr))
         for k=0, n_elements(framehdr)-1 do $
          expid[k] = sxpar(*framehdr[k], 'EXPOSURE')
         explist = expid[uniq(expid,sort(expid))]
         for iexp=0, n_elements(explist)-1 do begin
            ithis = where(expid EQ explist[iexp])
            combine1fiber, loglam1[*,ithis], flux1[*,ithis]+sky1[*,ithis], $
             invvar1[*,ithis], binsz=1d-4, newloglam=newloglam, $
             newflux=newflux1, newivar=newivar1
            newflux[*,i] = newflux1
            newivar[*,i] = newivar1
            i++
         endfor
      endfor
   endfor
   print

   ngoodpt = total(newivar NE 0,2)

   ;----------
   ; Mask out wavelengths outside the fitting domain
   ; and mask out the three [O_I] lines

   fitmask = (newloglam GE alog10(wfit[0])) * (newloglam LE alog10(wfit[1]))

   lambda = [5577.339, 6300.304, 6363.776]
   linesig = [4.5d-4, 4.5d-4, 3.5e-4]
   airtovac, lambda
   linemask = bytarr(n_elements(newloglam)) + 1B
   for i=0, n_elements(lambda)-1 do $
    linemask *= abs(newloglam - alog10(lambda[i])) GT linesig[i]

   for i=0L, ntot-1L do $
    newivar[*,i] *= fitmask * linemask

   ;----------
   ; Do the PCA solution

   res = pca_solve(newflux, newivar, niter=niter, nkeep=nkeep, usemask=usemask)

   ; Discard any wavelength with fewer than 20 good measurements
   ; (before my addditional masking)
   for ikeep=0, nkeep-1 do res[*,ikeep] *= (ngoodpt GE 10)

   ; Interpolate over the masked emission lines
   for ikeep=0, nkeep-1 do $
    res[*,ikeep] = djs_maskinterp(res[*,ikeep], linemask EQ 0, /const)

   ;----------
   ; Write a FITS file

   mwrfits, float(res), outfile, /create
   mwrfits, newloglam, outfile

   return
end
;------------------------------------------------------------------------------
