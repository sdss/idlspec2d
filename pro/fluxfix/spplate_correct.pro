;+
; NAME: 
;   spplate_correct
;
; PURPOSE:
;   Correct the flux calibration of the spPlate files for DR1
;
; CALLING SEQUENCE:
;   spplate_correct, plate, mjd, spectro_data_dir=, /remove_smear, $
;                   /zeropoint_cor
;
; INPUTS:
;    plate  -- plate ID number
;    mjd    -- MJD of the plate of interest
;
; OPTIONAL KEYWORDS
;    spectro_data_dir -- directory path to original spPlate files
;    remove_smear     -- set to remove the smear correction from all fibers
;    zeropoint_cor    -- correct the flux zeropoint as well as the color
;                          
; OUTPUTS:
;    All outputs are written to the current directory.  These include:
;    -- two FITS binary tables containing info about the standard stars called
;       spStd-$PLATE-$MJD-$SPECTROGRAPH.fits
;    -- a plot of the flux correction vectors called sphotocorr-$PLATE-$MJD.ps
;    -- a plot comparing magnitudes synthesized from the spectra with fiber
;       mags from photo (same as produced by idlspec2d) spSN2d-$PLATE-$MJD.ps
;    -- a modified version of the spPlate-$PLATE-$MJD.fits file
;    
; COMMENTS:
;    The spPlate files are modified to have improved spectrophotometry.  There
;    are two corrections: one which modifies the way in which the smear 
;    corrections are applied, and one which improves flux calibration. 
; 
;    The first correction which may be applied removes the effects of using
;    smear observations to calibrate the original data.  Because galaxies
;    are extended sources the 5x8 arcsecond smear aperture contains much more
;    light than the 3" science aperture.  It is not desirable to calibrate
;    the low order SED of galaxies to the light in the smear aperture 
;    because it could have a different color than that of the science aperture 
;    resulting in a mismatch between the low-order SED and the high order 
;    spectral features.  
;
;    Fortunately each fiber's SED can be easily corrected to match the "calib" 
;    image instead of the smear.  (The "calib" image is the highest S/N science 
;    image.)  This correction has the form of a 3rd order legendre polynomial
;    and is applied on a fiber to fiber basis.  Stars and other point sources 
;    are not greatly affected by the removal of the smear correction, since any 
;    light lost to atmospheric dispersion is effectively recovered by the flux 
;    calibration. 
;
;    The second correction which is applied is the flux correction itself.
;    A flux calibration "correction vector" is derived independently for each of 
;    the two spectrographs by ratioing the standard star spectra with stellar 
;    atmosphere models of the appropriate spectral type.  Each half-plate
;    is then multiplied by one correction vector (which is made from the 
;    average of many stars).  This correction is necessary, because the 
;    standard stars observed are not all identical in type to BD+17 4708, as 
;    was first assumed.
;  
;    By design, the flux correction changes the color of the spectra -- not
;    the absolute flux level.  The flux zeropoint will also be corrected if
;    the keyword "zeropoint_cor" is set AND the keyword "smear_remove" is
;    set.  When smear is removed, all of the spectra should have r-band 
;    synthetic magnitudes approximately equal to the fiber mag.  However there 
;    are sometimes slight offsets due to the manner in which the original 
;    Spectro2d flux calibration was done.  The median ratio of the
;    the r-band magnitudes synthesized from the spectra and the photo fiber
;    mags in the plugmap is used to reset the zeropoint, if this is desired.
; 
;    The corrected spPlate file is written out in the current directory.  The
;    HDU's which are modified are the flux (HDU 0), the inverse variance
;    (HDU 1), the AND and OR masks (HDU's 2 & 3), and the synthetic 
;    magnitudes (HDU 7).  
;
;    The Header keywords 'FLUXMOD' and 'SMEARMOD' and 'FLUXZPT' are added. 
;    The VERSCOMB keyword (which contains the version of idlspec2d used to 
;    combine the spectra) has the suffix 'FC' added and the keyword 'SMEARUSE' 
;    is set to 'F' (false).  Finally, the BUNIT keyword is modified to contain 
;    the correct flux units of 1E-17 erg/cm^2/s/Ang  (previously it was wrong).
;
; BUGS:
;    One caveat to keep in mind is that setting "remove_smear" is not exactly
;    equivalent to re-reducing the data without smear observations.  This 
;    is because the smear corrections are derived and applied to each exposure 
;    individually before they are combined and there is bound to be some 
;    noise in this process - esp. because the smear observations have low S/N.
;
;
; PROCEDURES CALLED:
;    fibermask_bits()
;    mrdfits()
;    mwrfits
;    plate_flux_recal
;    platesn
;    sxpar()
;    splog
;    undo_smear()
;
; REVISION HISTORY:
;    28-Oct-2002  Written by C. Tremonti (JHU)
;-
;-------------------------------------------------------------------------------

pro spplate_correct, plate, mjd, spectro_data_dir = spectro_data_dir, $
    remove_smear = remove_smear, zeropoint_cor = zeropoint_cor

  ;----------------------------------------------------------------------------
  ; Figure out file names and path
  ;----------------------------------------------------------------------------

  platestr = string(plate,format='(i4.4)')
  mjdstr = string(mjd,format='(i5.5)')

  if NOT keyword_set(spectro_data_dir) then $
    spectro_data_dir = getenv('SPECTRO_DATA')
  if (NOT keyword_set(spectro_data_dir)) then $
     splog, 'Environment variable SPECTRO_DATA must be set!'

  platename = 'spPlate-' + platestr + '-' + mjdstr + '.fits'
  orig_spplate = filepath(platename, root_dir=spectro_data_dir, $
                 subdirectory=platestr)

  ;----------------------------------------------------------------------------
  ; Read back all the HDU's from the plate in question
  ;----------------------------------------------------------------------------

  flux = mrdfits(orig_spplate, 0, hdr)
  invvar = mrdfits(orig_spplate, 1, invhdr)
  andmask = mrdfits(orig_spplate, 2, andhdr)
  ormask = mrdfits(orig_spplate, 3, orhdr)
  disp = mrdfits(orig_spplate, 4, disphdr)
  plugmap = mrdfits(orig_spplate, 5)
  snvec = mrdfits(orig_spplate, 6)
  synthmag = mrdfits(orig_spplate, 7)

  coeff0 = sxpar(hdr, 'COEFF0')
  coeff1 = sxpar(hdr, 'COEFF1')
  npix = sxpar(hdr, 'NAXIS1')
  loglam = coeff0 + coeff1 * lindgen(npix)
  wave = 10.0^loglam

  newflux = flux
  newinvvar = invvar

  ;----------------------------------------------------------------------------
  ; Determine the corrections that need to be applied to remove the effects
  ; of "smear" if this is desired
  ;----------------------------------------------------------------------------

  if keyword_set(remove_smear) then begin

    splog, 'Removing Smear from plate ' + platename 

    objtype = plugmap.objtype

    smearcor = undo_smear(plate, mjd, loglam, objtype, spectro_data_dir = $
                          spectro_data_dir)

    ;------------
    ; Correct the flux & inverse variance

    newflux = flux * smearcor
    newinvvar = invvar  / smearcor^2

    ;-------------
    ; Correct the mask bits to indicate no smear

    newandmask = andmask
    newormask = ormask

    bitval = fibermask_bits('SMEARIMAGE')
    bitset = where((newandmask and bitval) ne 0)
    if bitset[0] ne -1 then newandmask[bitset] = newandmask[bitset] - bitval  
    bitset = where((newormask and bitval) ne 0)
    if bitset[0] ne -1 then newormask[bitset] = newormask[bitset] - bitval  

    bitval = fibermask_bits('SMEARHIGHSN')
    bitset = where((newandmask and bitval) ne 0)
    if bitset[0] ne -1 then newandmask[bitset] = newandmask[bitset] - bitval  
    bitset = where((newormask and bitval) ne 0)
    if bitset[0] ne -1 then newormask[bitset] = newormask[bitset] - bitval  

    bitval = fibermask_bits('SMEARMEDSN')
    bitset = where((newandmask and bitval) ne 0)
    if bitset[0] ne -1 then newandmask[bitset] = newandmask[bitset] - bitval  
    bitset = where((newormask and bitval) ne 0)
    if bitset[0] ne -1 then newormask[bitset] = newormask[bitset] - bitval  

    ;--------------
    ; Fix up header keywords

    sxaddpar, hdr, 'SMEARMOD', 'T', 'Smear correction removed?'

    smearuse = sxpar(hdr, 'SMEARUSE', comment = comment)
    sxaddpar, hdr, 'SMEARUSE', 'F', comment
  endif 

  ;----------------------------------------------------------------------------
  ; Determine the flux correction from the standard stars on each half-plate.
  ; The fluxcor vectors are not saved to file -- if you want to write them to
  ; text then set the keyword "write" to 1.
  ;----------------------------------------------------------------------------

  plate_flux_recal, hdr, newflux, newinvvar, newandmask, plugmap, loglam, $
        /psplot, /stype, writeout = 0, fcor1 = fcor1, fcor2 = fcor2

  ;----------------------------------------------------------------------------
  ; Now apply the flux correction to the flux and invvar arrays -- "fcor" is 
  ; derived from the ratio of (old data / model) so we want to divide 
  ;----------------------------------------------------------------------------

  splog, 'Correcting Flux calibration for plate ' + platename

  for ifiber = 0, 319 do newflux[*,ifiber] = newflux[*,ifiber] / fcor1
  for ifiber = 320, 639 do newflux[*,ifiber] = newflux[*,ifiber] / fcor2

  for ifiber = 0, 319 do newinvvar[*,ifiber] = newinvvar[*,ifiber] * fcor1^2
  for ifiber = 320, 639 do newinvvar[*,ifiber] = newinvvar[*,ifiber] * fcor2^2
  
  sxaddpar, hdr, 'FLUXMOD', 'T', 'Flux Calibration Updated?'

  ;----------------------------------------------------------------------------
  ; Correct flux zero points if this is desired
  ;----------------------------------------------------------------------------

  if keyword_set(zeropoint_cor) AND keyword_set(remove_smear) then begin

    ;-------------------
    ; Compute the flux in the spectra (in ergs/s/cm^2/Hz) at the effective 
    ; wavelength of the r-filter

    flambda2fnu = (wave*wave / 2.99792e18) # replicate(1,640)

    fthru = filter_thru(newflux * flambda2fnu * 1e-17, waveimg=wave, $
                        mask=(newinvvar LE 0), /toair)

    r_spectro_fnu = fthru[*,2]   

    ;------------------
    ; Now compute the same from the photo fiber mag

    r_photo_mag = plugmap.mag[2] + 0.015  ; transform to AB 
    r_photo_fnu = 10.0^(-0.4 * (r_photo_mag + 48.6))

    ;------------------
    ; There ratio of the spectro and phot fluxes is used to set the zeropoint

    ifib1 = indgen(320)
    ifib2 = ifib1 + 320

    fcor1_zpt = median(r_photo_fnu[ifib1] / r_spectro_fnu[ifib1])
    fcor2_zpt = median(r_photo_fnu[ifib2] / r_spectro_fnu[ifib2])

    splog, 'Flux Zeropoint Correction to Spectrograph 1: ' + $
           string(fcor1_zpt, format = '(F6.3)')
    splog, 'Flux Zeropoint Correction to Spectrograph 2: ' + $
           string(fcor2_zpt, format = '(F6.3)')

    newflux[*,ifib1] = newflux[*,ifib1] * fcor1_zpt 
    newflux[*,ifib2] = newflux[*,ifib2] * fcor2_zpt 

    newinvvar[*,ifib1] = newinvvar[*,ifib1] / fcor1_zpt^2 
    newinvvar[*,ifib2] = newinvvar[*,ifib2] / fcor2_zpt^2 
    
    sxaddpar, hdr, 'FLUXZPT', 'T', 'Flux Zeropoint Updated?'

  endif

  ;----------------------------------------------------------------------------
  ; Re-determine the synthetic magnitudes (and make a nice diagnostic plot)
  ;----------------------------------------------------------------------------

   platesn, newflux, newinvvar, newandmask, plugmap, loglam, hdr=hdr, $
            plotfile= 'spSN2d-' + platestr + '-' + mjdstr + '.ps', $
            synthmag=newsynthmag 

  ;----------------------------------------------------------------------------
  ; Update the image header
  ;----------------------------------------------------------------------------

  verscomb = strtrim(sxpar(hdr, 'VERSCOMB', comment = comment), 2)
  sxaddpar, hdr, 'VERSCOMB', verscomb + '_FC', comment

  sxaddpar, hdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang', ' '

  ;----------------------------------------------------------------------------
  ; Write the modified data back to file
  ;----------------------------------------------------------------------------

  mwrfits, float(newflux), platename, hdr, /create
  mwrfits, float(newinvvar), platename, invhdr 
  mwrfits, newandmask, platename, andhdr 
  mwrfits, newormask, platename, orhdr 
  mwrfits, disp, platename, disphdr 
  mwrfits, plugmap, platename 
  mwrfits, snvec, platename 
  mwrfits, newsynthmag, platename 

end
