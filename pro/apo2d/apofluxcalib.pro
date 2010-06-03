;+
; NAME:
;   apofluxcalib
;
; PURPOSE:
;   Generate the flux-calibration vectors for use by APOPLOT.
;
; CALLING SEQUENCE:
;   apofluxcalib, [ platenum, mjd= ]
;
; INPUTS:
;   platenum - Plate number for obtaining the fluxing vectors; default to 406
;   mjd      - Modified Julian Date for above plate; default to plate 406
;              on MJD 51817 if neither is specified
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
;   For reference (excerpted from PR #6766):
;     Here is the list of 10 main-survey plates with the
;     highest S/N per minute in all cameras (where I've looked
;     up to plate 650 in the Spectro-2D v5 beta reductions):
;              398       51789
;              402       51793
;              406       51817
;              406       51900
;              411       51817
;              416       51811
;              418       51817
;              431       51877
;              436       51883
;              439       51877
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
;   readspec
;   traceset2xy
;
; REVISION HISTORY:
;   04-Dec-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro apofluxcalib, platenum, mjd=mjd

   if (NOT keyword_set(platenum)) then begin
      platenum = 406
      if (NOT keyword_set(mjd)) then mjd = 51817
   endif

   ; Determine the best exposure number
   readspec, platenum, mjd=mjd, objhdr=hdr
   if (NOT keyword_set(hdr)) then $
    message, 'No header found for PLATE=', platenum, ' MJD=', mjd
   expnum = sxpar(hdr, 'BESTEXP')
   splog, 'Choosing exposure number ', expnum

   root_dir = getenv('BOSS_SPECTRO_REDUX')
   platestr = string(platenum, format='(i4.4)')
   expstr = string(expnum, format='(i8.8)')

   camname = ['b1','b2','r1','r2']
   objfile = 'spFrame-' + camname + '-' + expstr + '.fits'
   fcalibfile = 'spFluxcalib-' + camname + '-' + expstr + '.fits'
   outfile = 'spFluxcalib-' + camname + '.fits'

   for icam=0, n_elements(camname)-1 do begin
      thisfile = (findfile( $
       filepath(objfile[icam], root_dir=root_dir, subdir=platestr)+'*'))[0]
      thatfile = (findfile( $
       filepath(fcalibfile[icam], root_dir=root_dir, subdir=platestr)+'*'))[0]
      if (keyword_set(thisfile) * keyword_set(thatfile) EQ 0) then $
       message, 'Could not find all input files'

      ; Read the wavelength solution and superflat
      objhdr = headfits(thisfile)
      wset = mrdfits(thisfile, 3)
      traceset2xy, wset, xpos, loglam

      ; Read the flux-calibration vector
      calibhdr = headfits(thatfile)
      cwavemin = sxpar(calibhdr, 'WAVEMIN')
      cwavemax = sxpar(calibhdr, 'WAVEMAX')
      calibset = mrdfits(thatfile, 1)
      calibfac = bspline_valu(loglam, calibset)

      ; Now produce the values by which we **divide** the raw spectra
      fluxvector = calibfac
      isort = sort(loglam)
      loglam = loglam[isort]
      fluxvector = fluxvector[isort]

      fset = bspline_iterfit(loglam, fluxvector, nord=4, bkspace=1.d-4)
      sxaddpar, calibhdr, 'EXPTIME', sxpar(objhdr,'EXPTIME')
      mwrfits, 0, outfile[icam], calibhdr, /create
      mwrfits, fset, outfile[icam]
   endfor

   return
end
;------------------------------------------------------------------------------
