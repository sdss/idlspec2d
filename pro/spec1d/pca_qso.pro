;------------------------------------------------------------------------------
pro pca_qso

   minwave = 3300
   maxwave = 8800
   snmax = 100
   niter = 10
   nkeep = 4

   get_juldate, jd
   mjdstr = string(long(jd-2400000L), format='(i5)')
   outfile = 'spEigenQSO-' + mjdstr + '.fits'

   ;----------
   ; Read the input spectra

   eigenfile = filepath('eigeninput_qso.dat', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, eigenfile, skip=2, plate, mjd, fiber, zfit, format='(L,L,L,D)'
;readcol, eigenfile, skip=2, plate, mjd, fiber, zfit, format='(L,L,L,D)', $
; numline=20 ; ???

   readspec, plate, fiber, mjd=mjd, flux=objflux, invvar=objivar, $
    andmask=andmask, plugmap=plugmap, loglam=objloglam

   nobj = (size(objflux, /dimens))[1]
   objdloglam = objloglam[1] - objloglam[0]

   if (keyword_set(snmax)) then begin
      ifix = where(objflux^2 * objivar GT snmax^2)
      if (ifix[0] NE -1) then objivar[ifix] = (snmax/objflux[ifix])^2
   endif

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-subtraction
   ; residuals.  Grow that mask by 2 pixels in each direction.

   skymask = (andmask AND pixelmask_bits('BRIGHTSKY')) NE 0
   for iobj=0, nobj-1 do $
    skymask[*,iobj] = smooth(float(skymask[*,iobj]),5) GT 0
   ibad = where(skymask)
andmask = 0 ; Free memory
   if (ibad[0] NE -1) then objivar[ibad] = 0

   ;----------

   pcaflux = pca_solve(objflux, objivar, objloglam, zfit, $
    niter=niter, nkeep=nkeep, newloglam=newloglam, eigenval=eigenval)

   indx = where(10^newloglam GE minwave AND 10^newloglam LE maxwave)
   sxaddpar, hdr, 'OBJECT', 'QSO'
   sxaddpar, hdr, 'COEFF0', newloglam[indx[0]]
   sxaddpar, hdr, 'COEFF1', objdloglam
   for i=0, n_elements(eigenval)-1 do $
    sxaddpar, hdr, 'EIGEN'+strtrim(string(i),1), eigenval[i]
   mwrfits, float(pcaflux[indx,*]), outfile, hdr, /create

   return
end
;------------------------------------------------------------------------------
