;+
; NAME:
;   apofluxcalib
;
; PURPOSE:
;   Generate the flux-calibration vectors for use by APOPLOT.
;
; CALLING SEQUENCE:
;   apofluxcalib
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The output files spFluxcalib-$CAMERA.fits should be moved to
;   the directory $IDLSPEC2D_DIR/examples for use by APOPLOT.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   bspline_iterfit()
;   bspline_valu()
;   headfits()
;   mrdfits()
;   mwrfits
;   traceset2xy
;
; REVISION HISTORY:
;   04-Dec-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro apofluxcalib

   root_dir = '/u/dss/spectro/test/0406'
   camname = ['b1','b2','r1','r2']
   flatfile = 'spFlat-' + camname + '-00006788.fits'
   arcfile = 'spArc-' + camname + '-00006789.fits'
   objfile = 'spFrame-' + camname + '-00006790.fits'
   fcalibfile = 'spFluxcalib-0406-51817-' + camname + '.fits'
   outfile = 'spFluxcalib-' + camname + '.fits'

   for icam=0, n_elements(camname)-1 do begin
      ; Read the wavelength solution and superflat
      wset = mrdfits(filepath(objfile[icam], root_dir=root_dir), 3)
      superfit = mrdfits(filepath(objfile[icam], root_dir=root_dir), 6)
      traceset2xy, wset, xpos, loglam

      ; Read the flux-calibration vector
      calibhdr = headfits(filepath(fcalibfile[icam], root_dir=root_dir))
      cwavemin = sxpar(calibhdr, 'WAVEMIN')
      cwavemax = sxpar(calibhdr, 'WAVEMAX')
      calibset = mrdfits(filepath(fcalibfile[icam], root_dir=root_dir), 1)
      calibfac = bspline_valu(loglam, calibset)

      ; Now produce the values by which we **divide** the raw spectra
      fluxvector = calibfac * superfit
      isort = sort(loglam)
      loglam = loglam[isort]
      fluxvector = fluxvector[isort]

      fset = bspline_iterfit(loglam, fluxvector, nord=4, bkspace=1.d-4)
      mwrfits, 0, outfile[icam], calibhdr, /create
      mwrfits, fset, outfile[icam]

   endfor

   return
end
;------------------------------------------------------------------------------
