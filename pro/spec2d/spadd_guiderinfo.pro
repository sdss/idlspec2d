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
   yanny_read, guidermonfile, pdata, stnames=stnames, /anonymous

   if (keyword_set(pdata)) then begin
      i = (where(stnames EQ 'GUIDEOBJ'))[0]
      if (i NE -1) then begin
         tags = tag_names(*pdata[i])
         ; Ensure that this guiderMon file has all of the following tag
         ; names, which wasn't the case for the early data.
         if ((where(tags EQ 'TIMESTAMP'))[0] NE -1 $
          AND (where(tags EQ 'FWHM'))[0] NE -1 $
          AND (where(tags EQ 'DRA'))[0] NE -1 $
          AND (where(tags EQ 'DDEC'))[0] NE -1) then begin
            guidermon = *pdata[i]
         endif else begin
            splog, 'WARNING: Invalid format for guiderMon file ' + guidermonfile
         endelse
      endif else begin
         splog, 'WARNING: No GUIDEOBJ entries in guiderMon file ' + guidermonfile
      endelse

      yanny_free, pdata
   endif else begin
      splog, 'WARNING: Empty guiderMon file ' + guidermonfile
   endelse

   see20 = 0.0
   see50 = 0.0
   see80 = 0.0
   rms20 = 0.0
   rms50 = 0.0
   rms80 = 0.0

   if (keyword_set(guidermon)) then begin
;      plate = sxpar(hdr, 'PLATEID')
;      found = where(guidermon.plateId EQ plate)
;      if (found[0] NE -1) then begin
;         guidermon = guidermon[found] 

      get_tai, hdr, tai_beg, tai_mid, tai_end
      taiplate = guidermon.timeStamp + 3506716800.0d
      found2 = where(taiplate GE tai_beg AND taiplate LE tai_end, ct2)

      if (ct2 GT 8) then begin   ; at least 8 fibers in the time interval?
         guidermon = guidermon[found2]
         seeing  = 0.0
         good_seeing = where(guidermon.fwhm GT 0.0,ct_seeing)
         if good_seeing[0] NE -1 then begin
            seeing = guidermon[good_seeing].fwhm
            seeing = seeing[sort(seeing)]           
            see20 = seeing[long(ct_seeing*0.2) - 1]
            see50 = seeing[long(ct_seeing*0.5) - 1]
            see80 = seeing[long(ct_seeing*0.8) - 1]
         endif else splog, 'Warning: No non-zero FWHM entries'

         rmsoff = 0.0
         good_rms = where(guidermon.dra NE 0.0 OR $
                          guidermon.ddec NE 0.0,ct_rms)
         if good_rms[0] NE -1 then begin
            rmsoff = sqrt(guidermon[good_rms].dra^2 +$
                          guidermon[good_rms].ddec^2)

            rmsoff = rmsoff[sort(rmsoff)]
            rms20 = rmsoff[long(ct_rms*0.2) - 1]
            rms50 = rmsoff[long(ct_rms*0.5) - 1]
            rms80 = rmsoff[long(ct_rms*0.8) - 1]
         endif else splog, 'Warning: No non-zero DRA and DDEC entries'
      endif else splog, 'Warning: No inclusive times ', tai_beg, tai_end
   endif

   sxaddpar, hdr, 'RMSOFF80', rms80, $
       ' 80% RMS offset of guide fibers (arcsec)', after='CAMERAS'
   sxaddpar, hdr, 'RMSOFF50', rms50, $
       ' 50% RMS offset of guide fibers (arcsec)', after='CAMERAS'
   sxaddpar, hdr, 'RMSOFF20', rms20, $
       ' 20% RMS offset of guide fibers (arcsec)', after='CAMERAS'
   sxaddpar, hdr, 'SEEING80', see80, ' 80% seeing during exposure (arcsec)', $
                             after='CAMERAS'
   sxaddpar, hdr, 'SEEING50', see50, ' 50% seeing during exposure (arcsec)', $
                              after='CAMERAS'
   sxaddpar, hdr, 'SEEING20', see20, ' 20% seeing during exposure (arcsec)', $
                             after='CAMERAS'

   return
end
;------------------------------------------------------------------------------
