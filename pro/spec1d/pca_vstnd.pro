;------------------------------------------------------------------------------
pro pca_vstnd

;   wavemin = 3300
;   wavemax = 8800
   snmax = 100
   niter = 10
   nkeep = 2

   get_juldate, jd
   mjdstr = string(long(jd-2400000L), format='(i5)')
   outfile = 'spEigenVstnd-' + mjdstr + '.fits'

   ;----------
   ; Read the input spectra

   eigenfile = filepath('eigeninput_vstnd.dat', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   djs_readcol, eigenfile, skip=2, plate, mjd, fiber, zfit, format='(L,L,L,D)'

   readspec, plate, fiber, mjd=mjd, flux=objflux, invvar=objivar, $
    andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=objloglam

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-sub residuals.

   objivar = skymask(objivar, andmask, ormask)
andmask = 0 ; Free memory
ormask = 0 ; Free memory

   nobj = (size(objflux, /dimens))[1]
   objdloglam = objloglam[1] - objloglam[0]

   if (keyword_set(snmax)) then begin
      ifix = where(objflux^2 * objivar GT snmax^2)
      if (ifix[0] NE -1) then objivar[ifix] = (snmax/objflux[ifix])^2
   endif

   ;----------

   pcaflux = pca_solve(objflux, objivar, objloglam, zfit, $
    wavemin=wavemin, wavemax=wavemax, $
    niter=niter, nkeep=nkeep, newloglam=newloglam, eigenval=eigenval)
   pcaflux = float(pcaflux)

   sxaddpar, hdr, 'OBJECT', 'VELOCITY_STANDARD'
   sxaddpar, hdr, 'COEFF0', newloglam[0]
   sxaddpar, hdr, 'COEFF1', objdloglam
   for i=0, n_elements(eigenval)-1 do $
    sxaddpar, hdr, 'EIGEN'+strtrim(string(i),1), eigenval[i]

   mwrfits, pcaflux, outfile, hdr, /create

   return
end
;------------------------------------------------------------------------------
