;+
; NAME:
;   plate_flux_recal
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   plate_flux_recal, hdr, flux, invvar, andmask, plugmap, loglam, $
;         /psplot, /stype, /writeout, fcor1 = fcor1, fcor2 = fcor2
;
; INPUTS:
;   hdr       - header of spPlate file
;   flux      - flux [npix, nfiber]
;   invvar    - inverse variance [npix, nfiber]
;   andmask   - and mask [npix, nfiber]
;   plugmap   - plugmap [nfiber]
;   loglam    - log10 of wavelength vector [npix]
; 
; OPTIONAL INPUT KEYWORDS
;   psplot      - generates a postscript plot named 
;                 sphotocorr-$PLATE-$MJD.ps
;   stype       - if set calls 'stype_standard' to spectral type the 
;                 standard stars.  If not set is assumed that the standards
;                 have already been typed and files of the format 
;                 spStd-$PLATE-$MJD-$SPECTROGRAPH.fits' exist
;   writeout    - if set, the flux correction vectors are saved as text files
;                 named sphotocorr-$PLATE-$MJD-$SPECTROGRAPH.txt
;  
; OPTIONAL OUTPUT KEYWORDS
;   fcor1       - flux correction vector for spectrograph 1
;   fcor2       - flux correction vector for spectrograph 2
;
; OUTPUT:  
;   If 'stype' is call A FITS binary table containing information on which  
;   KURUCZ models best fit the standards is produced.  If 'psplot' is called
;   A postscript plot showing the flux correction is generated.
;
; COMMENTS:   
;   Flux corrections are generated for one plate.  The flux, inverse variance,
;   andmask, plugmap, header, and wavelength array are given as input.  This
;   is so that the smear corrections can be removed first if desired.  A
;   separate correction is produced for each spectrograph and output in the
;   variables "fcor1" and "fcor2".  The corrections are produced by ratioing  
;   the standard star spectra with stellar atmosphere models of the appropriate 
;   spectral type.  The output vectors are the median smoothed average of 
;   the vectors generated for each of the stars (usually 8 per half plate).  
;   This correction is necessary, because the standard stars observed are not 
;   all identical in type to BD+17 4708, as was first assumed.
;
; BUGS:
;   The flux correction occasionally has abrupt jumps due to the outlier 
;   rejection algorithm in 'flux_standard' operating on a small number of spectra.
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   flux_standard
;   mrdfits()
;   stype_standard
;
; REVISION HISTORY:
;   16-Nov-2002  Written by C. Tremonti
;-
;------------------------------------------------------------------------------

pro plate_flux_recal, hdr, flux, invvar, andmask, plugmap, loglam, $
    psplot = psplot, stype = stype, writeout = writeout, fcor1 = fcor1, $
    fcor2 = fcor2
   
;******************************************************************************

  ;-----------------------------------------------------------------------------
  ; Get info from image header 
  ;-----------------------------------------------------------------------------

  platestr = string(sxpar(hdr, 'PLATEID'), format = '(I4.4)')
  mjdstr = string(sxpar(hdr, 'MJD'), format = '(I5)')

  sn1 = string(sxpar(hdr, 'SPEC1_R'), format='(F5.2)')
  sn2 = string(sxpar(hdr, 'SPEC2_R'), format='(F5.2)')
  seeing = string(sxpar(hdr, 'SEEING50'), format='(F4.2)')
  guiding = string(sxpar(hdr, 'RMSOFF50'), format='(F4.2)')
  airmass = string(sxpar(hdr, 'AIRMASS'), format='(F4.2)')

  title = 'Plate: ' + platestr + '  MJD: ' + mjdstr
  outprefix = 'sphotocorr-' + platestr + '-' + mjdstr + '-'

  print, platestr

  ;-----------------------------------------------------------------------------
  ; Use plugmap to find standard stars
  ;-----------------------------------------------------------------------------

  isphoto = where((strtrim(plugmap.objtype) EQ 'SPECTROPHOTO_STD' OR $
            strtrim(plugmap.objtype) EQ 'REDDEN_STD'), nsphoto)

  stdplug = plugmap[isphoto]

  sfib1 = stdplug[where(stdplug.spectrographid eq 1)].fiberid
  sfib2 = stdplug[where(stdplug.spectrographid eq 2)].fiberid

  ;-----------------------------------------------------------------------------
  ; Trim out spectra of standard stars
  ;-----------------------------------------------------------------------------

  flux1 = flux[*, sfib1 - 1]
  invvar1 = invvar[*, sfib1 - 1]
  andmask1 = andmask[*, sfib1 - 1]
  plug1 = plugmap[sfib1 - 1]

  flux2 = flux[*, sfib2 - 1]
  invvar2 = invvar[*, sfib2 - 1]
  andmask2 = andmask[*, sfib2 - 1]
  plug2 = plugmap[sfib2 - 1]

  ;-----------------------------------------------------------------------------
  ; Spectral type the standards (or read back the file of types if this is 
  ; already done)
  ;-----------------------------------------------------------------------------

  stdfiles = 'spStd-'+platestr+'-'+mjdstr+'-'+['1', '2']+'.fits'

  if keyword_set(stype) then begin
    stype_standard, loglam, flux1, invvar1, plug1, stdfiles[0]
    stype_standard, loglam, flux2, invvar2, plug2, stdfiles[1]
  endif

  stdstars1 = mrdfits(stdfiles[0], 1)
  stdstars2 = mrdfits(stdfiles[1], 1)

  ;-----------------------------------------------------------------------------
  ; Compute flux correction  -- "corvector" contains the corrections derived 
  ; for each star and "fcor" contains the average
  ;-----------------------------------------------------------------------------

  flux_standard, loglam, flux1, invvar1, andmask2, stdfiles[0], $
                 fcor = fcor1, corvector = corvector1, goodv = goodv1

  flux_standard, loglam, flux2, invvar2, andmask2, stdfiles[1], $
                 fcor = fcor2, corvector = corvector2, goodv = goodv2

  ;--------
  ; Normalize the corrections in the r-band

  wave = 10.0^loglam
  nwave = where(wave gt 5680 and wave lt 6680)
  norm1 = median(fcor1[nwave])
  norm2 = median(fcor2[nwave])

  fcor1 =  smooth(fcor1/norm1, 75) 
  fcor2 =  smooth(fcor2/norm2, 75) 

  corvector1 = corvector1 / norm1
  corvector2 = corvector2 / norm2

  ;-----------------------------------------------------------------------------
  ; Save correction as text file 
  ;-----------------------------------------------------------------------------

  if keyword_set(writeout) then begin
    openw, unit, outprefix + '1.txt', /get_lun 
    for i = 0, n_elements(wave) - 1 do $
      printf, unit, format = '(F8.3, F8.3)', wave[i], fcor1[i] 
    free_lun, unit

    openw, unit, outprefix + '2.txt', /get_lun 
    for i = 0, n_elements(wave) - 1 do $
      printf, unit, format = '(F8.2, F8.3)', wave[i], fcor2[i]
    free_lun, unit

  endif

  ;-----------------------------------------------------------------------------
  ; Plot Flux Vectors 
  ;-----------------------------------------------------------------------------

  if keyword_set(psplot) then begin

   ; Open Plot file
   set_plot, 'ps'
   device, filename = 'sphotocorr-' + platestr + '-' + mjdstr + '.ps', /color, $
           /portrait, xsize = 7.5, ysize = 10.0, xoffset=0.25, yoffset=0.5, $
           /inch

   !P.MULTI = [0, 1, 2]

    plot, wave, corvector1[*,0], xr=[3800, 9100], /xs, /ys, /nodata, $
          yr=[0.5, 1.4], xtitle = 'Wavelength (A)', charthick = 2, $
          ytitle = 'Flux Ratio (Calibrated SDSS Data / Model)', $
          title = title + '  SPEC: 1'

    color = strarr(n_elements(stdstars1)) + 'gray'
    color[where(goodv1)] = 'black'
    for j = 0, n_elements(stdstars1) - 1 do $
        djs_oplot, wave, corvector1[*,j], color = color[j], thick = 2, nsum=10
    djs_oplot, wave, fcor1, color='red', thick = 4, nsum=10

    xyouts, 4000, 1.32, 'S/N!U2!N(r) = ' + sn1
    xyouts, 4000, 0.55, 'Seeing = ' + seeing
    xyouts, 5200, 0.55, 'Guiding = '+ guiding
    xyouts, 6400, 0.55, 'Airmass = ' + airmass
    xyouts, 7700, 0.55, 'E(B-V) = ' + string(median([stdstars1.e_bv_sfd, $
            stdstars2.e_bv_sfd]), format = '(F4.2)')

    plot, wave, corvector2[*,0], xr=[3800, 9100], /xs, /ys, /nodata, $
          yr=[0.5, 1.4], xtitle = 'Wavelength (A)', charthick = 2, $
          ytitle = 'Flux Ratio (Calibrated SDSS Data / Model)', $
          title = title + '  SPEC: 2'

    color = strarr(n_elements(stdstars2)) + 'gray'
    color[where(goodv2)] = 'black'
    for j = 0, n_elements(stdstars2) - 1 do $
       djs_oplot, wave, corvector2[*,j], color=color[j], thick = 2, nsum=10  
    djs_oplot, wave, fcor2, color='red', thick = 4, nsum=10

    xyouts, 4000, 1.32, 'S/N!U2!N(r) = ' + sn2

    !P.MULTI = 0
    device, /close
    set_plot, 'x'
  endif

end
