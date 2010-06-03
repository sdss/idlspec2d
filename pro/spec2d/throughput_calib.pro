;+
; NAME:
;   throughput_calib
;
; PURPOSE:
;   Calculate spectroscopic throughput with fluxcalib files 
;
; CALLING SEQUENCE:
;  throughput = throughput_calib(file, loglam, $
;    primary=primary, exposure=exposure, pixelsize=pixelsize
;
; INPUTS:
;   file       - full name of fluxcalib fits file
;   loglam     - optional, alog10 wavelengths.
;                throughput is set to zero outside of the wavelength range
;                 of both bspline set and the wavemin,wavemax header cards
;
; OPTIONAL KEYWORDS:
;   primary    - effective primary aperature in cm^2 (default 40000)
;   exposure   - effective exposure time of calibration exposure
;   pixelsize  - physical CCD pixel size (guess alog10 1.0e-4)
;                 this could also be passed as a keyword in the header
;
; OUTPUTS:
;   throughput - percent throughput at each loglam position
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;   loglam = 3.57 + findgen(400)/1000.
;   cd,getenv('BOSS_SPECTRO_REDUX')
;   allfiles = findfile('0*/*calib*fits', count=nf)
;   throughput = fltarr(400, nf)
;   for i=0,nf-1 do throughput[*,i] = throughput_calib(allfiles[i],loglam)
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   26-Jun-2001  Scott Burles (Paris)
;-
;------------------------------------------------------------------------------
function throughput_calib, file, loglam, $
    primary=primary, exposure=exposure, pixelsize=pixelsize

   hdr = headfits(file)
   tt = mrdfits(file,1)

   minwave = sxpar(hdr,'WAVEMIN') ? alog10(sxpar(hdr,'WAVEMIN')) : $
                                    min(tt.fullbkpt)
   maxwave = sxpar(hdr,'WAVEMAX') ? alog10(sxpar(hdr,'WAVEMAX')) : $
                                    max(tt.fullbkpt)

   if NOT keyword_set(loglam) then $
     loglam = findgen(500) / 499.0 * (maxwave - minwave) + minwave 

   throughput = loglam * 0.0

   good = where(loglam GE minwave AND loglam LE maxwave)
   if good[0] EQ -1 then return, throughput

   glam = loglam[good]

   calib = bspline_valu(glam, tt)

   ; calib is in units of counts/pix/(10^-17 ergs/s/cm^2/Ang)
   
   if NOT keyword_set(pixelsize) then $
         pixelsize = alog(10.0) * 1.0e-4 * 10^glam
 
   if NOT keyword_set(primary) then primary = 40000.0  ;primary aperture
   if NOT keyword_set(exposure) then begin
      if sxpar(hdr,'CALTIME') then exposure=sxpar(hdr,'CALTIME') $
      else exposure = 55.0   ;guess at smear exposure
   endif
 
   ; factor is number of photons per 10^-17 erg at glam:
 
   factor = (10^glam/12400.0) / 1.609e5 
   throughput[good] = calib / factor / pixelsize / primary / exposure

   if max(throughput[good]) GT 0.5 then begin
     splog,'Throughput over 50%, likely wrong exposure time'
     splog,'Switching time from '+string(exposure)+' to 900 seconds'
     throughput[good] = throughput[good] * exposure / 900.0
   endif

   return, throughput
end
  
   
