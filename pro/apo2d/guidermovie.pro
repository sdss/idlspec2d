;+
; NAME:
;   guidermovie
;
; PURPOSE:
;   Plot the guider images for an entire night using ATV.
;
; CALLING SEQUENCE:
;   guidermovie, [ mjd=, _EXTRA=KeywordsForATV ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - MJD number; if not set, then select the most recent MJD
;                in the $RAWDATA_DIR directory.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $RAWDATA_DIR/$MJD/guider/gimg*.fits.gz
;
; PROCEDURES CALLED:
;   atv
;   atvxyouts
;   djs_filepath()
;   fileandpath()
;   mrdfits()
;   splog
;
; REVISION HISTORY:
;   13-May-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro guidermovie, mjd=mjd, _EXTRA=KeywordsForATV

   ;----------
   ; If MJD is not specified, then find the most recent MJD for output files

   rawdata_dir = getenv('RAWDATA_DIR')
   if (NOT keyword_set(rawdata_dir)) then $
    rawdata_dir = '/data/spectro'

   if (NOT keyword_set(mjd)) then begin
      mjdlist = get_mjd_dir(rawdata_dir, mjstart=1, mjend=99999, mjd='?????')
      mjd = (reverse(mjdlist[sort(mjdlist)]))[0]
      splog, 'Selecting MJD=', mjd, ' (override this with MJD keyword)'
   endif

   ;----------
   ; Get file list for guider images (g-zipped only)

   gfiles = findfile( djs_filepath('gimg*fits.gz', root_dir=rawdata_dir, $
    subdirectory=[string(mjd,format='(i5.5)'),'guider']), count=nfile )
   if (nfile EQ 0) then begin
      splog, 'No files found
      return
   endif

   ;----------
   ; Read the dark image

   dark = mrdfits('~jeg/Data/guider/darks/30s/SUMDKBIN.FIT')

   gfiles = gfiles[sort(gfiles)]
   ifile = 0
   cc = 'L'

   while (cc NE 'Q') do begin

      if (cc EQ 'L') then begin
         ifile = (ifile + 1) < (nfile - 1)
         c1 = strupcase(get_kbrd(0))
         if (keyword_set(c1)) then cc = c1
      endif

      if (cc NE 'L') then begin
         print, 'Press b=back one frame'
         print, '      f=forward one frame'
         print, '      l=loop mode'
         print, '      q=quit (and enter interactive mode for this plot)'
         cc = strupcase(get_kbrd(1))
         case cc of
            'B': ifile = (ifile - 1) > 0
            'F': ifile = (ifile + 1) < (nfile - 1)
            else:
         endcase
      endif

      img = mrdfits(gfiles[ifile], 0, hdr, /silent) + 32768.
      utstring = sxpar(hdr, 'UTTIME')
      ww = strsplit(utstring, ': ', /extract)
      hr = long(ww[0])
      min = long(ww[1])
      sec = long(ww[2])
      if (ww[3] EQ 'PM') then hr = hr + 12
      utstring = string(mjd, hr, min, sec, $
       format='(i5.5,"  ",i2.2,":",i2.2,":",i2.2," UT")')

      splog, ifile, nfile, fileandpath(gfiles[ifile]), utstring, $
       format='("File ",i4," of ",i4,": ",a,2x,a)'

      ; Dark-subtract
      if (total(size(dark,/dimens) - size(img,/dimens)) EQ 0) then begin
         res = linfit(dark[*,0:50], img[*,0:50])
         img = img - res[0] - res[1] * dark
      endif

      atv, img, _EXTRA=KeywordsForATV
      atvxyouts, 0, 0, utstring, color='red', charsize=5

   endwhile
end
;------------------------------------------------------------------------------
