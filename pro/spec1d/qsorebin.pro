;------------------------------------------------------------------------------
pro qsorebin

   objdloglam = 1.d-4
   wavemin = 525.
   wavemax = 9300.
   minuse = 40

   ;----------
   ; Read the template files
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

   tfile = filepath('qso.template', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   djs_readcol, tfile, twave, tflux, terr, lq, uq, usemask, $
    format='(D,F,F,F,F,L)'
   tloglam0 = alog10(twave[0])
   tdloglam = alog10(twave[1] / twave[0])
   igood = where(terr GT 0)
   tivar = 0 * terr
   tivar[igood] = 1. / (terr[igood])^2

   ;----------
   ; Set the new wavelength mapping here...

   newloglam = wavevector(alog10(wavemin), alog10(wavemax), binsz=objdloglam)

   ;----------
   ; Re-bin to the same pixel bin size.

   newflux = interpol(tflux, alog10(twave), newloglam)
   newivar = interpol(tivar, alog10(twave), newloglam)
   newuse = interpol(float(usemask), alog10(twave), newloglam)

   ;----------
   ; Fill in bad data with a running median of good data

   qgood = newuse GE minuse AND newivar GT 0
   igood = where(qgood, ngood)
   ibad = where(qgood EQ 0, nbad)
   medflux = 0 * newflux
   if (nbad GT 0) then begin
      medflux[igood] = djs_median(newflux[igood], width=51, boundary='nearest')
      medflux = djs_maskinterp(medflux, qgood EQ 0, /const)
      newflux[ibad] = medflux[ibad]
   endif

   ;----------
   ; Write output file

   sxaddpar, hdr, 'COEFF0', newloglam[0]
   sxaddpar, hdr, 'COEFF1', objdloglam
   mwrfits, float(newflux), 'spEigenQSO.fits', hdr, /create

   return
end
;------------------------------------------------------------------------------

