;+
; NAME:
;   undo_smear
;
; PURPOSE:
;   Create correction vectors for 1 plate which can be used to remove the
;   the effects of "smear" from the spectra.
;
; CALLING SEQUENCE:
;   undo_smear, plate, mjd, loglam, objtype
;
; INPUTS:
;    plate   -- plate ID number
;    mjd     -- MJD of the plate of interest
;    loglam  -- log10(wavelength in Angstroms) [npix]
;    objtype -- vector of object types from plugmap [nfiber]  
;
; OUTPUTS:
;   A vector of corrections is returned with dimensions [npix, nfiber].  
;   The flux array should be *multiplied* by this correction. 
;
; COMMENTS:
;   Spectro2d currently measures the ratio of each science exposure to the 
;   "smear" exposure and stores the results as a set of third order 
;   polynomials in files named 'spFluxcorr-$EXPOSURE-$SPECTROGRAPH.fits'
;   The ratios are used by spectro2d to correct every fiber such that it's
;   low order SED matches that of the smear.  This is desirable for point
;   sources because the may have lost light from the 3" science aperture which
;   is recovered in the 5"x8" effective aperture of the smear.  However, this
;   correction is NOT appropriate for extended sources such as galaxies since 
;   the larger smear aperture sees a different part of the object which could
;   have different colors than the region observed in the science aperture.
;
;   Fortunately it is easy to undo the smear correction.  We can tie the
;   SEDs of the spectra to the "calib" image instead of the smear.  The
;   "calib" image is simply defined as the highest S/N science exposure.
;   This piece of code figures out which exposure is the smear and which is 
;   the calib by reading the logfile spDiagcomb-$PLATE-$MJD.log. The ratio
;   of smear/calib is found in the spFluxcorr-$EXPOSURE-$SPECTROGRAPH.fits 
;   file for the calib exposure.  To correct for the effects of smear we
;   want to multiply the object flux by the low order fit to the ratio of 
;   calib and smear images (the inverse of the "spFluxcorr" vector).  
;
;   Which objects should have smear removed is the matter of some debate, 
;   since it depends on just how point like or extended the objects are.  
;   Thus this program does not apply the correction for smear, it simply 
;   returns it for each fiber in the plate.
;
;   Because the smear and calib images have different exposure times, their 
;   ratio is not close to unity, even for sources which have lost/gained 
;   very little light.  However, we want our smear correction to be mostly 
;   a change in color rather than a change in total flux since the spPlate
;   files are already flux calibrated.  We achieve this by normalizing the
;   the smear correction such that the mean for the spectro-photo standards is 
;   ~1 in the r-band.  The standard stars are identified as those objects 
;   having objtype = '*_STD'.
; 
; BUGS:
;
; PROCEDURES CALLED:
;    meanclip()
;    mrdfits()
;
; REVISION HISTORY:
; 28-Oct-2002  Written by C. Tremonti (JHU)
;-
;-------------------------------------------------------------------------------

function undo_smear, plate, mjd, loglam, objtype, spectro_data_dir = $
         spectro_data_dir

  platestr = string(plate,format='(i4.4)')
  mjdstr = string(mjd,format='(i5.5)')

  if not keyword_set(spectro_data_dir) then $
    spectro_data_dir = getenv('SPECTRO_DATA')
  if (NOT keyword_set(spectro_data_dir)) then $
      message, 'Environment variable SPECTRO_DATA must be set!'

  ;-----------------------------------------------------------------------------
  ; Grep the log file to determine the smear and calib images used
  ;-----------------------------------------------------------------------------

  logfile = filepath('spDiagcomb-' + platestr + '-' + mjdstr + '.log', $
            root_dir=spectro_data_dir, subdirectory=platestr)

  spawn, 'grep "Select smear image" ' + logfile, smeartxt
  spawn, 'grep "Select calib image" ' + logfile, calibtxt

  words = strsplit(smeartxt[0], '#', /extract)
  smearexp1 = long(words[1])
  words = strsplit(smeartxt[1], '#', /extract)
  smearexp2 = long(words[1])

  words = strsplit(calibtxt[0], '#', /extract)
  calibexp1 = long(words[1])
  words = strsplit(calibtxt[1], '#', /extract)
  calibexp2 = long(words[1])

  ;-----------------------------------------------------------------------------
  ; If the smear and calib images are the same then return all 1's
  ;-----------------------------------------------------------------------------

  npix = n_elements(loglam)
  smearcor = fltarr(npix, 640) + 1

  if (smearexp1 eq calibexp1) and (smearexp2 eq calibexp2) then return, smearcor

  ;-----------------------------------------------------------------------------
  ; Read in the "Fluxcorr" vectors for the calib images (these are a 
  ; polynomial fit to the ratio of smear/calib)
  ;-----------------------------------------------------------------------------
   
  calibfile1 = filepath('spFluxcorr-' + string(calibexp1, format = '(I8.8)') $
               + '-1.fits', root_dir=spectro_data_dir, subdirectory=platestr)
  calibcorset1 = mrdfits(calibfile1, 1)
  traceset2xy, calibcorset1, rebin(loglam, npix, 320), calibcorimg1

  calibfile2 = filepath('spFluxcorr-' + string(calibexp2, format = '(I8.8)') $
               + '-2.fits', root_dir=spectro_data_dir, subdirectory=platestr)
  calibcorset2 = mrdfits(calibfile2, 1)
  traceset2xy, calibcorset2, rebin(loglam, npix, 320), calibcorimg2

  calibcorimg = [[calibcorimg1], [calibcorimg2]] 

  ;-----------------------------------------------------------------------------
  ; Normalize the corrections such that the mean of the spectro photo stars
  ; is 1 in the r-band.  (This is necessary because of the different exposure
  ; times of the smear and calib images.)
  ;-----------------------------------------------------------------------------

  stds = where(strmatch(objtype, '*STD*'), nstds)
  rmed = fltarr(nstds)
  normwl = where(10.0^loglam gt 5800 and 10.0^loglam lt 6800)

  for istd = 0, nstds - 1 do $
     rmed[istd] = median(calibcorimg[normwl, stds[istd]])

  meanclip, rmed, rmean

  calibcorimg = calibcorimg / rmean

  smearcor = 1.0 / calibcorimg

  return, smearcor
end
