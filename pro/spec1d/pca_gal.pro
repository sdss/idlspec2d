;------------------------------------------------------------------------------
pro pca_gal

   wavemin = 1850. ; Good to z=1
   wavemax = 9300.
   snmax = 100
   niter = 10
   nkeep = 4
   minuse = 10

   get_juldate, jd
   mjdstr = string(long(jd-2400000L), format='(i5)')
   outfile = 'spEigenGal-' + mjdstr + '.fits'
   plotfile = 'spEigenGal-' + mjdstr + '.ps'

   dfpsplot, plotfile, /color, /landscape
   colorvec = ['default', 'red', 'green', 'blue']

   ;----------
   ; Read the input spectra

   eigenfile = filepath('eigeninput_gal.dat', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   djs_readcol, eigenfile, skip=2, plate, mjd, fiber, zfit, format='(L,L,L,D)'

   readspec, plate, fiber, mjd=mjd, flux=objflux, invvar=objivar, $
    andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=objloglam

   ;----------
   ; Insist that all of the requested spectra exist

   imissing = where(plugmap.fiberid EQ 0, nmissing)
   if (nmissing GT 0) then begin
      for i=0, nmissing-1 do $
       print, 'Missing plate=', slist[imissing[i]].plate, $
        ' mjd=', slist[imissing[i]].mjd, $
        ' fiber=', slist[imissing[i]].fiberid
      message, string(nmissing) + ' missing object(s)'
   endif

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

   pcaflux = pca_solve(objflux, objivar, objloglam, zfit, $
;    wavemin=wavemin, wavemax=wavemax, $
    niter=niter, nkeep=nkeep, newloglam=newloglam, eigenval=eigenval, $
    acoeff=acoeff, usemask=usemask)
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
   ; Make plots

   djs_plot, 10^newloglam, pcaflux[*,0], $
    xrange=minmax(10^newloglam), yrange=minmax(pcaflux), /xstyle, $
    color=colorvec[0], $
    xtitle='Wavelength [Ang]', ytitle='Flux [arbitrary units]', $
    title='Galaxies: Eigenspectra', /xlog
   for i=1, nkeep-1 do $
    djs_oplot, 10^newloglam, pcaflux[*,i], $
     color=colorvec[i MOD n_elements(colorvec)]

   aratio10 = acoeff[1,*] / acoeff[0,*]
   aratio20 = acoeff[2,*] / acoeff[0,*]
   djs_plot, aratio10, aratio20, /nodata, $
    xtitle='Eigenvalue Ratio (a_1/a_0)', $
    ytitle='Eigenvalue Ratio (a_2/a_0)', $
    title='Galaxies: Eigenvalue Ratios'
   for j=0, n_elements(aratio10)-1 do $
    djs_xyouts, aratio10[j], aratio20[j], align=0.5, $
     string(plate[j], fiber[j], format='(i4,"-",i3)'), $
     color=colorvec[j MOD n_elements(colorvec)]

   ;----------
   ; Write output file

   sxaddpar, hdr, 'OBJECT', 'GALAXY'
   sxaddpar, hdr, 'COEFF0', newloglam[0]
   sxaddpar, hdr, 'COEFF1', objdloglam
   for i=0, n_elements(eigenval)-1 do $
    sxaddpar, hdr, 'EIGEN'+strtrim(string(i),1), eigenval[i]

   mwrfits, pcaflux, outfile, hdr, /create

   dfpsclose

   return
end
;------------------------------------------------------------------------------
