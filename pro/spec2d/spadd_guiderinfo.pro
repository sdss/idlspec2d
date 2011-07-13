;+
; NAME:
;   spadd_guiderinfo
;
; PURPOSE:
;   Add seeing and RMS offset header cards by parsing guiderMon-$MJD.par
;
; CALLING SEQUENCE:
;   spadd_guiderinfo, hdr
;
; INPUTS:
;   hdr        - FITS header
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   hdr        - FITS header (modified)
;
; COMMENTS:
;   This routine recalls the beginning and ending exposure time-stamps
;   and collects all seeing and guider offsets available from the 
;   speclog file: guiderMon-mmmmm.par (where mmmmm is MJD).
;
;   It reports the 20,50,80% seeing and RMS guiding deviations
;   (closely resembles -1,0,+1 sigma for normal distributions).
;   The following keywords are added: RMSOFF20, RMSOFF50, RMSOFF80,
;   SEEING20, SEEING50, SEEING80.
;   All quantities are quoted in arcseconds.
;   
;   Results of 0.0 mean entries do not exist or are ill-defined.
;
; EXAMPLES:
;   filename = 'sdR-b2-00003976.fit'
;   hdr = headfits(filename)
;   spadd_guiderinfo, hdr
;
; BUGS:
;
; PROCEDURES CALLED:
;   concat_dir()
;   fileandpath()
;   get_tai
;   splog
;   sxaddpar
;   sxpar()
;   yanny_free
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES
;
; DATA FILES:
;   $SPECLOG_DIR/$MJD/guiderMon-$MJD.par
;
; REVISION HISTORY:
;   28-Jan-2002  Written by S. Burles, MIT
;-
;------------------------------------------------------------------------------
pro spadd_guiderinfo, hdr

   if (n_params() LT 1) then begin
      print, 'Syntax - spadd_guiderinfo, hdr'
      return
   endif

   speclog_dir = getenv('SPECLOG_DIR')
   if (NOT keyword_set(speclog_dir)) then $
    message, 'Must set environment variable SPECLOG_DIR'
   mjd = sxpar(hdr, 'MJD')
   mjdstr = string(mjd, format='(i05.5)')
   plugdir = concat_dir(speclog_dir, mjdstr)

   guidermonfile = filepath('guiderMon-'+mjdstr+'.par', root_dir=plugdir)
   guidermon = yanny_readone(guidermonfile, 'GUIDEOBJ', $
    stnames=stnames, /anonymous)

   if (keyword_set(guidermon)) then begin
         ; Ensure that this guiderMon file has all of the following tag
         ; names, which wasn't the case for the early data.
         if (tag_exist(guidermon, 'timestamp') $
          AND tag_exist(guidermon, 'fwhm') $
          AND tag_exist(guidermon, 'dra') $
          AND tag_exist(guidermon, 'ddec')) then begin
            guidermon = guidermon
         endif else begin
            guidermon = 0
            splog, 'WARNING: Invalid format for guiderMon file ' + guidermonfile
         endelse
   endif else begin
      splog, 'WARNING: Empty guiderMon file ' + guidermonfile
   endelse

   see20 = 0.0
   see50 = 0.0
   see80 = 0.0
   rms20 = 0.0
   rms50 = 0.0
   rms80 = 0.0
   nguide = 0

   if (keyword_set(guidermon)) then begin
; No need to introduce the dependence below upon the plate numbers
; being correct in both the image headers and the guiderMon file.
;      plate = sxpar(hdr, 'PLATEID')
;      ifound = where(guidermon.plateId EQ plate)
;      if (ifound[0] NE -1) then guidermon = guidermon[ifound] 

      get_tai, hdr, tai_beg, tai_mid, tai_end
      ; If this is an old SDSS-I file, then the timstamps are all offset
      ; and need correcting; if BOSS, then no corrections
      taiplate = guidermon.timestamp
      if (tag_exist(guidermon, 'focusoffset') EQ 0) then $
       taiplate += 3506716800.0d
      ifound = where(taiplate GE tai_beg AND taiplate LE tai_end, nfound)

      if (nfound GT 8) then begin   ; at least 8 fibers in the time interval?
         guidermon = guidermon[ifound]

         ; Count the number of guider frames based upon the unique number
         ; of time stamps.
         alltimes = guidermon.timestamp
         alltimes = alltimes[uniq(alltimes, sort(alltimes))]
         nguide = n_elements(alltimes)
         splog, 'Number of guider frames = ', nguide

         seeing  = 0.0
         qgood = (guidermon.fwhm GT 0) AND finite(guidermon.fwhm)
         if (tag_exist(guidermon, 'focusoffset')) then $
          qgood *= (abs(guidermon.focusoffset) LT 100) ; in-focus fibers only
         if (tag_exist(guidermon, 'exists')) then $
          qgood *= (guidermon.exists EQ 'T')
         if (tag_exist(guidermon, 'enabled')) then $
          qgood *= (guidermon.enabled EQ 'T')
         igood = where(qgood, ngood)
         if (ngood GT 0) then begin
            seeing = guidermon[igood].fwhm
            seeing = seeing[sort(seeing)]           
            see20 = seeing[(long(ngood*0.20) - 1) > 0]
            see50 = seeing[(long(ngood*0.50) - 1) > 0]
            see80 = seeing[(long(ngood*0.80) - 1) > 0]
         endif else splog, 'Warning: No non-zero FWHM entries'

         rmsoff = 0.0
         qgood = guidermon.dra NE 0.0 OR guidermon.ddec NE 0.0
         if (tag_exist(guidermon, 'exists')) then $
          qgood *= (guidermon.exists EQ 'T')
         if (tag_exist(guidermon, 'enabled')) then $
          qgood *= (guidermon.enabled EQ 'T')

         ;- Ensure offsets aren't NaN or Inf
         qgood *= finite(guidermon.dra) * finite(guidermon.ddec)

         igood = where(qgood, ngood)
         if (ngood GT 0) then begin
            rmsoff = sqrt(guidermon[igood].dra^2 + guidermon[igood].ddec^2)
            rmsoff = rmsoff[sort(rmsoff)]
            rms20 = rmsoff[(long(ngood*0.20) - 1) > 0]
            rms50 = rmsoff[(long(ngood*0.50) - 1) > 0]
            rms80 = rmsoff[(long(ngood*0.80) - 1) > 0]
         endif else splog, 'Warning: No non-zero DRA and DDEC entries'
      endif else splog, 'Warning: Too few inclusive times ', tai_beg, tai_end
   endif

   sxaddpar, hdr, 'RMSOFF80', rms80, $
    ' 80% RMS offset of guide fibers (arcsec)', after='CAMERAS'
   sxaddpar, hdr, 'RMSOFF50', rms50, $
    ' 50% RMS offset of guide fibers (arcsec)', after='CAMERAS'
   sxaddpar, hdr, 'RMSOFF20', rms20, $
    ' 20% RMS offset of guide fibers (arcsec)', after='CAMERAS'
   sxaddpar, hdr, 'SEEING80', see80, $
    ' 80% seeing during exposure (arcsec)', after='CAMERAS'
   sxaddpar, hdr, 'SEEING50', see50, $
    ' 50% seeing during exposure (arcsec)', after='CAMERAS'
   sxaddpar, hdr, 'SEEING20', see20, $
    ' 20% seeing during exposure (arcsec)', after='CAMERAS'
   sxaddpar, hdr, 'NGUIDE', nguide, $
    ' Number of guider frames during exposure', after='CAMERAS'

   return
end
;------------------------------------------------------------------------------
