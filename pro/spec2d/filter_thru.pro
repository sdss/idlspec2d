;+
; NAME:
;   filter_thru
;
; PURPOSE:
;   Compute throughput in SDSS filters
;
; CALLING SEQUENCE:
;   res = filter_thru( flux, [waveimg=, wset=, mask=, filter_prefix=, /toair ])
;
; INPUTS:
;   flux       - Flux image [NX,NTRACE]
;
; OPTIONAL INPUTS:
;
; OPTIONAL KEYWORDS:
;   waveimg    - Wavelength image in Angstroms [NX,NTRACE], or this could
;                be a single vector if the wavelength mapping is the same
;                for all traces (note the latter is faster to compute)
;   wset       - Wavelength solution in log-lambda; required if WAVEIMG not set
;   mask       - Linearly interpolate over pixels where MASK is nonzero.
;                [NX,NTRACE]
;   filter_prefix  - Use alternate prefix for filter curves to use
;                    (allowed are sdss or doi) [sdss]
;   toair      - Convert the wavelengths to air from vacuum before computing
;
; OUTPUTS:
;   res        - Integrated response in all 5 SDSS filters, ordered ugriz;
;                dimensions are [NTRACE,5] or [5] if NTRACE=1.
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Needs waveimg to be equally spaced in log lambda (MRB 4.5.01) ???
;
; PROCEDURES CALLED:
;   vactoair
;   djs_maskinterp()
;   readcol
;   traceset2xy
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/sdss_u_atm.dat
;   $IDLSPEC2D_DIR/etc/sdss_g_atm.dat
;   $IDLSPEC2D_DIR/etc/sdss_r_atm.dat
;   $IDLSPEC2D_DIR/etc/sdss_i_atm.dat
;   $IDLSPEC2D_DIR/etc/sdss_z_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_u_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_g_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_r_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_i_atm.dat
;   $IDLSPEC2D_DIR/etc/doi_z_atm.dat
;
; REVISION HISTORY:
;   10-Mar-2000  Written by D. Schlegel, Princeton
;   05-Apr-2001  Modified by Michael Blanton to allow alternate filters
;-
;------------------------------------------------------------------------------
function filter_thru, flux, waveimg=waveimg, wset=wset, mask=mask, norm=norm, $
 filter_prefix=filter_prefix, toair=toair

   dims = size(flux, /dimens)
   nx = dims[0]
   if (N_elements(dims) EQ 1) then ntrace = 1 $
    else ntrace = dims[1]

   if (NOT keyword_set(filter_prefix)) then filter_prefix='sdss'
   ffiles = [filter_prefix+'_u_atm.dat', filter_prefix+'_g_atm.dat', $
             filter_prefix+'_r_atm.dat', filter_prefix+'_i_atm.dat', $
             filter_prefix+'_z_atm.dat']
   nfile = N_elements(ffiles)

   if (ntrace EQ 1) then res = fltarr(1, nfile) $
    else res = fltarr(ntrace, nfile)

   ;----------
   ; Get the wavelength at each pixel in FLUX

   if (NOT keyword_set(waveimg)) then begin
      traceset2xy, wset, pixnorm, logwave
      newwaveimg = 10^logwave
   endif else begin
      newwaveimg = waveimg
   endelse
   
   ;----------
   ; Convert wavelengths from vacuum to air if necessary

   if (keyword_set(ttoair)) then $
    vactoair, newwaveimg

   ;----------
   ; Interpolate over masked or low-S/N pixels in each spectrum

   if (keyword_set(mask)) then $
    flux_interp = djs_maskinterp(flux, mask, iaxis=0, /const)

   ;----------
   ; Integrate over each filter response curve

   for ifile=0, nfile-1 do begin

      filename = filepath(ffiles[ifile], $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
      readcol, filename, fwave, fthru, /silent

      filtimg = 0.0 * flux
      if (size(newwaveimg,/n_dimen) EQ 1) then $
       filtimg[*] = interpol(fthru, fwave, newwaveimg) # replicate(1,ntrace) $
      else $
       filtimg[*] = interpol(fthru, fwave, newwaveimg) 

      if (keyword_set(mask)) then $
       res[*,ifile] = total(flux_interp * filtimg, 1) $
      else $
       res[*,ifile] = total(flux * filtimg, 1)

      if (keyword_set(norm)) then begin
         sumfilt = total(filtimg,1)
         res[*,ifile] = res[*,ifile] / (sumfilt + (sumfilt LE 0))
      endif

   endfor

   return, res
end
;------------------------------------------------------------------------------
