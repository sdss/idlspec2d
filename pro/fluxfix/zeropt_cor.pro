;+
; NAME: 
;   zeropt_cor
;
; PURPOSE:
;   Compute a zeropoint correction to the flux calibration 
;
; CALLING SEQUENCE:
;   zeropt_cor, hdr, wave, flux, invvar, plugmap, fcor1_zpt = fcor1_zpt, $
;   fcor2_zpt, fcor2_zpt, tsobj_dir = tsobj_dir, fibermags = fibermags
;
; INPUTS:
;    hdr  -- a string array containing the image header
;    wave  -- the wavelength array in Angstroms (not log10) [npix]
;    flux  -- the flux array [npix, nfiber]
;    invvar -- the inverse varriance [npix, nfiber]
;    plugmap -- an array of structures containing the plugmap infor [nfiber]
;    
; OUTPUTS:
;    fcor1_zpt -- zero point correction to be applied to spectrograph 1
;                 (scalar).  Multiply by this numer.
;    fcor2_ztp -- zero point correction to be applied to spectrograph 2
;    tsobj_dir -- directory containing tsObj files.  If set the tsObj 
;                 file with the same mapname is used as the source of the 
;                 r-band photo mag
;    fibermags -- ugriz fiber magnitudes from the tsObj file 
;                          
; COMMENTS:
;   This program computes a zeropoint correction to the spectral flux from 
;   the ratio of the r-band flux in the spectra and the r-band flux measured
;   by the photo fiber mag.  By default the photo fiber mag comes from the
;   plugmap, but if the "tsobj_dir" keyword is set, a tsObjfile is sought
;   matching the mapname of the plate.  If a tsObj file is found then the
;   r-band magnitudes from this file are used instead of those from the
;   plugmap.  This allows the spectra to be calibrated against the most 
;   recent Photo data.  This program produces two scalar variables "fcor_zpt1"  
;   and "fcor_zpt2" that contain the correction which should be applied to 
;   the spectral flux a follows:
;   newflux[*,0:319] = flux[*,0:319] * fcor_zpt1
;   newflux[*,320:639] = flux[*,320:639] * fcor_zpt2
;   
; BUGS:
;   The tsObj files are assumed to be in the format of those produced by
;   Jeff Munn.  Specifically, the tsObjname is the mapname, and the Photo5.3
;   data is in HDU #2
;
; PROCEDURES CALLED:
;    mrdfits()
;    sxpar()
;    splog
;
; REVISION HISTORY:
;    21-Nov-2002  Written by C. Tremonti (JHU)
;-
;-------------------------------------------------------------------------------

pro zeropt_cor, hdr, wave, flux, invvar, plugmap, fcor1_zpt = fcor1_zpt, $
    fcor2_zpt, fcor2_zpt, tsobj_dir = tsobj_dir

  ;----------------------------------------------------------------------------
  ; Get info from image header
  ;----------------------------------------------------------------------------

  plateid = string(sxpar(hdr, 'PLATEID'), format = '(I4.4)')
  mjd = string(sxpar(hdr, 'MJD'), format = '(I5)')
  platename = 'spPlate-' + plateid + '-' + mjd

  mapname = strtrim(sxpar(hdr, 'NAME'), 2)

  ;----------------------------------------------------------------------------
  ; Determine if we are using the plugmap or tsobj for the photo mags
  ;----------------------------------------------------------------------------

  ; By default use the plugmap
  r_photo_mag = plugmap.mag[2]

  ; If the tsobj_dir is set look for a tsobj file corresponding to the mapname
  if keyword_set (tsobj_dir) then begin
     tsobjname = filepath('tsObj-' + mapname + '.fit', root_dir = tsobj_dir)
     if file_test(tsobjname) then begin

        ; Files from Jeff Munn have target data in HDU1 and photo5.3 in HDU2
        tsobj = mrdfits(tsobjname, 2)

        fibermags = tsobj.fibercounts
        r_photo_mag = tsobj.fibercounts[2]

        splog, 'Using r-band magnitudes from tsObj-' + mapname + $
               ' for ' + platename
     endif else begin
        splog, 'tsObj-' + mapname + ' not found for ' + platename
        splog, 'Using plugmap file for zeropoint calibration instead'
     endelse
  endif

  ;----------------------------------------------------------------------------
  ; Compute the flux in the spectra (in ergs/s/cm^2/Hz) at the effective 
  ; wavelength of the r-filter
  ;----------------------------------------------------------------------------

  flambda2fnu = (wave*wave / 2.99792e18) # replicate(1,640)

  fthru = filter_thru(flux * flambda2fnu * 1e-17, waveimg=wave, $
                      mask=(invvar LE 0), /toair)

  r_spectro_fnu = fthru[*,2]   

  ;----------------------------------------------------------------------------
  ; Now compute the same from the photo fiber mag
  ;----------------------------------------------------------------------------

  r_photo_mag = r_photo_mag + 0.015  ; transform to AB 
  r_photo_fnu = 10.0^(-0.4 * (r_photo_mag + 48.6))

  ;----------------------------------------------------------------------------
  ; There ratio of the spectro and photo fluxes is used to set the zeropoint
  ;----------------------------------------------------------------------------

  ifib1 = indgen(320)
  ifib2 = ifib1 + 320

  fcor1_zpt = median(r_photo_fnu[ifib1] / r_spectro_fnu[ifib1])
  fcor2_zpt = median(r_photo_fnu[ifib2] / r_spectro_fnu[ifib2])

  splog, 'Flux Zeropoint Correction to ' + platename + ' Spectrograph 1: ' + $
          string(fcor1_zpt, format = '(F6.3)')
  splog, 'Flux Zeropoint Correction to ' + platename + ' Spectrograph 2: ' + $
          string(fcor2_zpt, format = '(F6.3)')

end
