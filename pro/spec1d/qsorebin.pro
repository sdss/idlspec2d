pro qsorebin

   objdloglam = 1.d-4

   ;----------
   ; Read the template files
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

   tfile = filepath('qso.template', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   djs_readcol, tfile, twave, tflux, terr, lq, uq, npt, $
    format='(D,F,F,F,F,L)'
   tloglam0 = alog10(twave[0])
   tdloglam = alog10(twave[1] / twave[0])
   igood = where(terr GT 0)
   tivar = 0 * terr
   tivar[igood] = 1. / (terr[igood])^2

   ; Test if the binning is the same as for the objects.
   ; If not, then re-bin to the same pixel bin size.

   if (abs((tdloglam-objdloglam)/objdloglam) GT 1.e-5) then begin
      ii = where(tivar GT 0 AND npt GT 10 AND twave GT 800)
      tloglam0 = alog10(twave[ii[0]])
      nnewbin = alog10(max(twave[ii]) / min(twave[ii])) / objdloglam
      newloglam = tloglam0 + lindgen(nnewbin) * objdloglam
      newflux = interpol(tflux[ii], alog10(twave[ii]), newloglam)
      newivar = interpol(tivar, alog10(twave), newloglam)
   endif

   sxaddpar, hdr, 'COEFF0', newloglam[0]
   sxaddpar, hdr, 'COEFF1', objdloglam
   mwrfits, newflux, 'spEigenQSO.fits', hdr, /create

stop
end
