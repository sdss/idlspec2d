;+
; NAME:
;   fcalib_default
;
; PURPOSE:
;   Return a default flux-calibration vector
;
; CALLING SEQUENCE:
;   calibfac = fcalib_default(camname, loglam, [exptime] )
;
; INPUTS:
;   camname    - Camera name; 'b1', 'b2', 'r1' or 'r2'
;   loglam     - Wavelengths [log-10 Angstroms]
;
; OPTIONAL INPUTS:
;   exptime    - Rescale the fluxing vector to this exposure time [sec];
;                default to the value in the spFluxcalib file header
;                (probably 900 sec)
;
; OUTPUTS:
;   calibfac   - Flux-calibration values in units of
;                electrons/(e-17 erg/cm^2/Ang); set to 0 if outside
;                the knonwn wavelengths; the spectra in electrons (not ADU)
;                should be divided by this; note that this calibration vector
;                is for the specified exposure time, so there is not time unit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/spFluxcalib-$CAMERA.fits
;
; PROCEDURES CALLED:
;   bspline_valu()
;   headfits()
;   mrdfits()
;   sxpar()
;
; REVISION HISTORY:
;   17-Dec-2005  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function fcalib_default, camname, loglam, exptime

   if (n_params() LT 2) then $
    message, 'Must specify CAMNAME and LOGLAM'

   ; Read a canonical flux-calibration vector
   fcalibfile = filepath('spFluxcalib-'+camname[0]+'.fits', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
   fcalibfile = (findfile(fcalibfile+'*'))[0]
   if (NOT keyword_set(fcalibfile)) then $
    message, 'No default flux-calib file found for CAMERA=', camname
   calibhdr = headfits(fcalibfile)
   cwavemin = sxpar(calibhdr, 'WAVEMIN')
   cwavemax = sxpar(calibhdr, 'WAVEMAX')
   calibset = mrdfits(fcalibfile, 1, /silent)
   calibfac = bspline_valu(loglam, calibset) $
    * (loglam GE alog10(cwavemin) AND loglam LE alog10(cwavemax))
   if (keyword_set(exptime)) then $
    calibfac = calibfac * exptime / sxpar(calibhdr,'EXPTIME')

   return, calibfac
end
;------------------------------------------------------------------------------
