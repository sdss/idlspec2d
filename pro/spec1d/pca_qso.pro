;------------------------------------------------------------------------------
pro pca_qso

   wavemin = 460
   wavemax = 9300
   snmax = 100
   niter = 10
   nkeep = 4
   minuse = 3

   get_juldate, jd
   mjdstr = string(long(jd-2400000L), format='(i5)')
   outfile = 'spEigenQSO-' + mjdstr + '.fits'
   plotfile = 'spEigenQSO-' + mjdstr + '.ps'

   dfpsplot, plotfile, /color
   colorvec = ['default', 'red', 'green', 'blue']

   ;----------
   ; Read the input spectra

   eigenfile = filepath('eigeninput_qso.dat', $
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
   ; Set the new wavelength mapping here...

   newloglam = wavevector(alog10(wavemin), alog10(wavemax), binsz=objdloglam)

   ;----------
   ; Do PCA solution -- all components at once.

   ; The following would solve for all the eigen-vectors at once.
   ; This can result in an unphysical 1st eigencomponent, probably
   ; because each spectrum only covers a small range of rest wavelength.

;   pcaflux = pca_solve(objflux, objivar, objloglam, zfit, $
;    wavemin=wavemin, wavemax=wavemax, $
;    niter=niter, nkeep=nkeep, newloglam=newloglam, eigenval=eigenval)

   ;----------
   ; Do PCA solution -- but one component at a time.

   for ikeep=0, nkeep-1 do begin
      splog, 'Solving for eigencomponent #', ikeep+1, ' of ', nkeep

      pcaflux1 = pca_solve(objflux, objivar, objloglam, zfit, $
       wavemin=wavemin, wavemax=wavemax, $
       niter=niter, nkeep=1, newloglam=newloglam, $
       newflux=newflux, newivar=newivar, eigenval=eigenval1, $
       acoeff=acoeff1, usemask=usemask)
      if (ikeep EQ 0) then saveflux = newflux

      if (ikeep EQ 0) then begin
         pcaflux = pcaflux1
         eigenval = eigenval1
;         acoeff = acoeff1
      endif else begin
         pcaflux = [[pcaflux], [pcaflux1]]
; These eigenvalues are not normalized properly!!!???
         eigenval = [eigenval, eigenval1]
;         acoeff = [acoeff, acoeff1]
      endelse

      ; Re-solve for the coefficients using all PCA components so far
      acoeff = fltarr(ikeep+1,nobj)
      for iobj=0, nobj-1 do begin
         junk = computechi2(saveflux[*,iobj], sqrt(newivar[*,iobj]), $
          pcaflux, acoeff=theta)
         acoeff[*,iobj] = theta
      endfor

      objloglam = 0 ; Prevent re-binning of spectra on subsequent calls
                    ; to PCA_SOLVE().
;      objflux = newflux - acoeff1 ## pcaflux1
      objflux = saveflux - acoeff ## pcaflux
      objivar = newivar

   endfor
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
    title='QSOs: Eigenspectra', /xlog
   for i=1, nkeep-1 do $
    djs_oplot, 10^newloglam, pcaflux[*,i], $
     color=colorvec[i MOD n_elements(colorvec)]

   aratio10 = acoeff[1,*] / acoeff[0,*]
   aratio20 = acoeff[2,*] / acoeff[0,*]
   djs_plot, aratio10, aratio20, /nodata, $
    xtitle='Eigenvalue Ratio (a_1/a_0)', $
    ytitle='Eigenvalue Ratio (a_2/a_0)', $
    title='QSOs: Eigenvalue Ratios'
   for j=0, n_elements(aratio10)-1 do $
    djs_xyouts, aratio10[j], aratio20[j], align=0.5, $
     string(plate[j], fiber[j], format='(i4,"-",i3)'), $
     color=colorvec[j MOD n_elements(colorvec)]

   ;----------
   ; Write output file

   sxaddpar, hdr, 'OBJECT', 'QSO'
   sxaddpar, hdr, 'COEFF0', newloglam[0]
   sxaddpar, hdr, 'COEFF1', objdloglam
   for i=0, n_elements(eigenval)-1 do $
    sxaddpar, hdr, 'EIGEN'+strtrim(string(i),1), eigenval[i]

   mwrfits, pcaflux, outfile, hdr, /create

   dfpsclose

   return
end
;------------------------------------------------------------------------------
