;------------------------------------------------------------------------------
pro pca_lrgtest, platenums, nkeep=nkeep, $
 zmin=zmin, zmax=zmax, outfile=outfile, loz=loz, midz=midz, hiz=hiz

   wavemin = 3600.
   wavemax = 6200.
   snmax = 100
   niter = 10
   if (NOT keyword_set(nkeep)) then nkeep = 2
   minuse = 5

   if (NOT keyword_set(platenums)) then begin
      platelist, plist=plist
      platenums = plist[ where(plist.qsurvey) ].plate
   endif

   if (keyword_set(loz)) then begin
      outfile = 'spLRG-loz.fits'
      zmin = 0.15
      zmax = 0.25
   endif else if (keyword_set(midz)) then begin
      outfile = 'spLRG-midz.fits'
      zmin = 0.30
      zmax = 0.40
   endif else if (keyword_set(hiz)) then begin
      outfile = 'spLRG-hiz.fits'
      zmin = 0.40
      zmax = 0.50
   endif

   if (NOT keyword_set(zmin)) then zmin = 0.0
   if (NOT keyword_set(zmax)) then zmax = 0.5
   if (NOT keyword_set(outfile)) then $
    outfile = string(long(zmin*100), long(zmax*100), $
     format='("spLRG_", i3.3, "_", i3.3, ".fits")')

   ;----------
   ; Read the input spectra

   readspec, platenums, zans=zans, plug=plug
   indx = where(((plug.primtarget AND 2LL^5) NE 0 $
           OR (plug.primtarget AND 2LL^26) NE 0) $
    AND strtrim(zans.class) EQ 'GALAXY' $
    AND zans.z GE zmin AND zans.z LE zmax)
   zans = zans[indx]
   splog, 'Number of objects = ', n_elements(zans)

   readspec, zans.plate, zans.fiberid, mjd=zans.mjd, $
    flux=objflux, invvar=objivar, $
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
   ; Set the new wavelength mapping here...

   newloglam = wavevector(alog10(wavemin), alog10(wavemax), binsz=objdloglam)

   ;----------
   ; Do PCA solution

   pcaflux = pca_solve(objflux, objivar, objloglam, zans.z, $
;    wavemin=wavemin, wavemax=wavemax, $
    niter=niter, nkeep=nkeep, newloglam=newloglam, eigenval=eigenval, $
    usemask=usemask)
   pcaflux = float(pcaflux)

   ;----------
   ; Fill in bad data with a running median of good data

   qgood = usemask GE minuse
   igood = where(qgood, ngood)
   ibad = where(qgood EQ 0, nbad)
   medflux = 0 * pcaflux
   if (nbad GT 0) then begin
      for i=0, nkeep-1 do begin
         medflux[igood,i] = $
          djs_median(pcaflux[igood,i], width=51, boundary='nearest')
         medflux[*,i] = djs_maskinterp(medflux[*,i], qgood EQ 0, /const)
      endfor
      pcaflux[ibad,*] = medflux[ibad,*]
   endif

   ;----------
   ; Write output file

   sxaddpar, hdr, 'OBJECT', 'LRG z=' + string(zmin) + ' ' + string(zmax)
   sxaddpar, hdr, 'COEFF0', newloglam[0]
   sxaddpar, hdr, 'COEFF1', objdloglam
   sxaddpar, hdr, 'NGALAXY', n_elements(zans), ' Number of galaxies in fit'
   for i=0, n_elements(eigenval)-1 do $
    sxaddpar, hdr, 'EIGEN'+strtrim(string(i),1), eigenval[i]
   for i=0, n_elements(platenums)-1 do $
    sxaddhist, 'Plate '+string(platenums[i]), hdr

   mwrfits, pcaflux, outfile, hdr, /create

   return
end
;------------------------------------------------------------------------------
