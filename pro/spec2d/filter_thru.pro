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
;                    (allowed are sdss, doi, sdss_jun2001) [sdss_jun2001]
;   toair      - Convert the wavelengths to air from vacuum before computing
;
; OUTPUTS:
;   res        - Integrated response in all 5 SDSS filters, ordered ugriz;
;                dimensions are [NTRACE,5] or [5] if NTRACE=1.
;
; COMMENTS:
;   The filter curve files are assumed to be in $IDLUTILS_DIR/data/filters.
;
; EXAMPLES:
;
; BUGS:
;   Needs waveimg to be equally spaced in log lambda (MRB 4.5.01) ???
;    Now calculates pixel size in log lambda to do correct spectrophotometry
;
; PROCEDURES CALLED:
;   vactoair
;   djs_maskinterp()
;   readcol
;   traceset2xy
;
; REVISION HISTORY:
;   10-Mar-2000  Written by D. Schlegel, Princeton
;   05-Apr-2001  Modified by Michael Blanton to allow alternate filters
;-
;------------------------------------------------------------------------------
function filter_thru, flux, waveimg=waveimg, wset=wset, mask=mask, $
 filter_prefix=filter_prefix, toair=toair

   dims = size(flux, /dimens)
   nx = dims[0]
   if (N_elements(dims) EQ 1) then ntrace = 1 $
    else ntrace = dims[1]

   if (NOT keyword_set(filter_prefix)) then filter_prefix='sdss_jun2001'
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

   if (keyword_set(toair)) then $
    vactoair, newwaveimg

   logwave = alog10(newwaveimg)
   diffx   = (findgen(nx-1) + 0.5) # replicate(1,ntrace)
   diffy   = logwave[1:*,*] - logwave[0:nx-2,*]
   xy2traceset, diffx, diffy, diffset, ncoeff=4, xmin=0, xmax=nx-1., /silent
   traceset2xy, diffset, pixnorm, logdiff

   logdiff = abs(logdiff)

   ;----------
   ; Interpolate over masked or low-S/N pixels in each spectrum

   if (keyword_set(mask)) then $
    flux_interp = djs_maskinterp(flux, mask, iaxis=0, /const)

   ;----------
   ; Integrate over each filter response curve

   for ifile=0, nfile-1 do begin

      filename = filepath(ffiles[ifile], $
       root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
      readcol, filename, fwave, fthru, /silent

      ; Use the Goddard routine LINTERP instead of the IDL built-in INTERPOL.
;      filtimg = 0.0 * flux
;      if (size(newwaveimg,/n_dimen) EQ 1) then $
;       filtimg[*] = interpol(fthru, fwave, newwaveimg) # replicate(1,ntrace) $
;      else $
;       filtimg[*] = interpol(fthru, fwave, newwaveimg) 
      linterp, fwave, fthru, newwaveimg[*], filtimg
      filtimg = filtimg * logdiff

      if (size(newwaveimg,/n_dimen) EQ 1) then $
       filtimg = filtimg # replicate(1,ntrace) $
      else $
       filtimg = reform(filtimg, size(flux,/dimens))


      if (keyword_set(mask)) then $
       res[*,ifile] = total(flux_interp * filtimg, 1) $
      else $
       res[*,ifile] = total(flux * filtimg, 1)

      sumfilt = total(filtimg,1)
      res[*,ifile] = res[*,ifile] / (sumfilt + (sumfilt LE 0))
   endfor

   return, res
end

;------------------------------------------------------------------------------
