;+
; NAME:
;   spthroughput
;
; PURPOSE:
;   Return the throughput and efficiency for a given plate+camera+exposure
;
; CALLING SEQUENCE:
;  photons_per_flux_per_sec = spthroughput( plate, [ indx ], camname=, $
;   expnum=, [ loglam=, exptime=, airmass=, efficiency=, /median ]
;
; INPUTS:
;   plate      - Plate number
;   camname    - Camera name, 'b1', 'b2', 'r1', or 'r2'
;   expnum     - Exposure number
;
; OPTIONAL INPUTS:
;   indx       - Fiber index number(s) 0-indexed; default to [0,N-1]
;                where N is the number of fibers in one spectrograph
;   loglam     - Wavelengths in vacuum log10-Angstroms; if specified, then this
;                is a vector of the same wavelengths for all fibers
;   median     - If set, and LOGLAM is input, then median-filter all
;                of the fibers at each wavelength
;
; OUTPUTS:
;   photons_per_flux_per_sec - Calibration vector(s) of photons per
;                (10^-17 erg/s/cm^2/Ang); return 0 if the flux-calibration
;                failed for this camera and default values were used
;
; OPTIONAL OUTPUTS:
;   loglam     - Wavelength scale in vacuum log-10-Angstroms for the native
;                pixel scale, if LOGLAM is not specified; this is an
;                array with the same dimensions as PHOTONS_PER_FLUX_PER_SEC
;   exptime    - Exposure time (seconds)
;   airmass    - Airmass for this exposure
;   efficiency - Fractional efficiency using the parameters of the SDSS
;                telescope mirror sizes; return 0 if the flux-calibration
;                failed for this camera and default values were used
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   The fiber flats (FFLAT) should also be included.
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   correct_dlam
;   djs_filepath()
;   mrdfits()
;   spframe_read
;   sxpar()
;
; REVISION HISTORY:
;   29-Mar-2006  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function spthroughput, plate, indx1, camname=camname, expnum=expnum, $
 loglam=loglam, exptime=exptime, airmass=airmass, $
 efficiency=efficiency, median=median

   efficiency = 0 ; default return value

   if (n_elements(plate) NE 1) then $
    message, 'PLATE must be a scalar!'
   if (size(camname,/tname) NE 'STRING') then $
    message, 'CAMNAME is not valid!'
   if (NOT keyword_set(expnum)) then $
    message, 'EXPNUM must be set!'
   if (keyword_set(loglam)) then $
    if (size(loglam,/n_dimen) NE 1) then $
     message, 'LOGLAM must be 1-dimensional vector!'

   outdir = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
    + '/' + string(plate,format='(i4.4)')
   filename1 = djs_filepath('spFrame-' $
    +string(format='(a2,a,i8.8,a)',camname,'-',expnum, $
    '.fits*'), root_dir=outdir)
   filename2 = djs_filepath('spFluxcalib-' $
    +string(format='(a2,a,i8.8,a)',camname,'-',expnum, $
    '.fits*'), root_dir=outdir)
   spframe_read, filename1, loglam=loglam1, hdr=objhdr, wset=wset, $
    superflat=superflat
   if (NOT keyword_set(loglam1)) then return, 0
   calibfac = mrdfits(filename2, 0, /silent)
   if (NOT keyword_set(calibfac)) then return, 0
   flatexp = long(strmid(sxpar(objhdr, 'FLATFILE'),7,8))
   exptime = sxpar(objhdr, 'EXPTIME')
   if (arg_present(airmass)) then airmass = sxpar(objhdr, 'AIRMASS')

   if (n_elements(indx1) GT 0) then indx = indx1 $
    else indx = lindgen(sxpar(objhdr,'NAXIS2'))
   nfiber = n_elements(indx)

   correct_dlam, calibfac, 0, wset, /inverse
   loglam1 = loglam1[*,indx]
   calibfac = superflat[*,indx] * calibfac[*,indx] / exptime
   if (keyword_set(loglam)) then begin
      photons_per_flux_per_sec = fltarr(n_elements(loglam),nfiber)
      minlog = min(loglam1, max=maxlog)
      for ifiber=0L, nfiber-1L do $
       photons_per_flux_per_sec[*,ifiber] = interpol(calibfac[*,ifiber], $
        loglam1[*,ifiber], loglam) * (loglam GE minlog AND loglam LE maxlog)
      if (keyword_set(median) AND nfiber GT 1) then begin
         photons_per_flux_per_sec = djs_median(photons_per_flux_per_sec, 2)
         nfiber = 1
      endif
   endif else begin
      photons_per_flux_per_sec = calibfac
      loglam = loglam1
   endelse

   if (arg_present(efficiency)) then begin
      ; 2.5-m telescope with 1.3-m obscured by secondary baffling
      tel_area = !pi * (125.^2 - 65.0^2) ; cm^2
      wave = 10^loglam
      npix = (size(loglam, /dimens))[0]
      dloglam = abs([loglam[1:npix-1,*]-loglam[0:npix-2,*], $ 
       loglam[npix-1,*]-loglam[npix-2,*]])
      dwave = wave * alog(10.) * dloglam
      efficiency = fltarr(npix,nfiber)
      if (size(photons_per_flux_per_sec, /n_dimen) GT 1) then begin
         wave = rebin(wave, npix, nfiber)
         dwave = rebin(dwave, npix, nfiber)
      endif
      efficiency = 1e17 * photons_per_flux_per_sec / dwave / tel_area $
       * 6.62e-27 * 3.e10 / (wave * 1e-8)
   endif

   return, photons_per_flux_per_sec
end
;------------------------------------------------------------------------------
