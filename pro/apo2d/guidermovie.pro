;+
; NAME:
;   guidermovie
;
; PURPOSE:
;   Plot the guider images for an entire night using ATV.
;
; CALLING SEQUENCE:
;   guidermovie, [ mjd=, plate=, expnum=, _EXTRA=KeywordsForATV ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - MJD number; if not set, then select the most recent MJD
;                in the $RAWDATA_DIR directory.
;   plate      - Plate number; if specified, then select all exposures
;                within the time frame for the science+smear exposures
;                for this plate during the night in question.
;   expnum     - Exposure numbers (for the gimg*.fits guider images)
;                as an array.
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
;   $RAWDATA_DIR/$MJD/sdR-b1-????????.fit.gz
;
; PROCEDURES CALLED:
;   atv
;   atvxyouts
;   djs_filepath()
;   fileandpath()
;   get_tai
;   jdcnv
;   mrdfits()
;   splog
;   sdsshead()
;
; REVISION HISTORY:
;   13-May-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro guidermovie, mjd=mjd, plate=plate, expnum=expnum, _EXTRA=KeywordsForATV

   quiet = !quiet
   !quiet = 1

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
   mjdstr = string(mjd, format='(i5.5)')

   ;----------
   ; Get file list for guider images (g-zipped only)

   gfiles = findfile( djs_filepath('gimg*fits.gz', root_dir=rawdata_dir, $
    subdirectory=[string(mjd,format='(i5.5)'),'guider']), count=nfile )
   if (nfile EQ 0) then begin
      splog, 'No files found
      return
   endif

   ;----------
   ; Trim to only specified exposure numbers (and update NFILE variable)

   if (keyword_set(expnum)) then begin
      qkeep = bytarr(nfile)
      for ifile=0, nfile-1 do begin
         thisexp = long( strmid(fileandpath(gfiles[ifile]),5,4) )
         if ((where(expnum EQ thisexp))[0] NE -1) then qkeep[ifile] = 1b
      endfor
      ikeep = where(qkeep, nfile)
      if (nfile EQ 0) then begin
         splog, 'No files matching EXPNUM'
         return
      endif
      gfiles = gfiles[ikeep]
   endif

   ;----------
   ; Read the headers for all the guider files

   print, 'Reading guider image headers...'
   monthname = ['','Jan','Feb','Mar','Apr','May','Jun','Jul', $
    'Aug','Sep','Oct','Nov','Dec']
   datestring = strarr(nfile)
   timearray = dblarr(nfile)
   for ifile=0, nfile-1 do begin
      hdr = headfits(gfiles[ifile])

      utdate = sxpar(hdr,'UTDATE')
      ww = strsplit(utdate, ', ', /extract)
      year = long(ww[3])
      month = (where(ww[1] EQ monthname))[0] > 0
      date = long(ww[2])

      uttime = sxpar(hdr,'UTTIME')
      ww = strsplit(uttime, ': ', /extract)
      hr = long(ww[0])
      min = long(ww[1])
      sec = long(ww[2])

      datestring[ifile] = string(year, month, date, hr, min, sec, $
       format='(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2," Z")')
      jdcnv, year, month, date, hr + min/60. + sec/3600., thisjd
      timearray[ifile] = (thisjd - 2400000.5D) * (24.D*3600.D)

      print, format='("File ",i5," of ",i5,a1,$)', $
       ifile+1, nfile, string(13b)
   endfor

   ;----------
   ; Sort these files by date+time

   isort = sort(datestring)
   gfiles = gfiles[isort]
   datestring = datestring[isort]
   timearray = timearray[isort]

   ;----------
   ; If PLATE is specified, then read the headers of all the sdR files,
   ; determine which times span the science+smear exposures, and trim
   ; the guider images to that list.

   if (keyword_set(plate)) then begin

      ;----------
      ; Set input directory for sdR files

      rawdata_dir = getenv('RAWDATA_DIR')
      if (NOT keyword_set(rawdata_dir)) then $
       rawdata_dir = '/data/spectro'

      ;----------
      ; Find sdR files for b1 camera only for this MJD (g-zipped only!)

      sdrname = findfile(filepath('sdR-b1-????????.fit.gz', $
       root_dir=rawdata_dir, subdir=mjdstr), count=nsdr)

      if (nsdr EQ 0) then begin
         print, 'No sdR files found for MJD=' + mjdstr
         return
      endif

      print, 'Reading sdR image headers...'
      for isdr=0, nsdr-1 do begin
         qdone = fits_wait(sdrname[isdr], deltat=1, tmax=1, /header_only)
         if (qdone) then begin
            thishdr = sdsshead(sdrname[isdr])
            thisplate = long(sxpar(thishdr,'PLATEID'))
            thisflavor = strtrim(sxpar(thishdr,'FLAVOR'),2)
            if (thisplate EQ plate AND (thisflavor EQ 'science' $
             OR thisflavor EQ 'smear')) then begin
               get_tai, thishdr, tai_beg, tai_mid, tai_end
               if (NOT keyword_set(tai_first)) then begin
                  tai_first = tai_beg
                  tai_last = tai_end
               endif else begin
                  tai_first = tai_first < tai_beg
                  tai_last = tai_last > tai_end
               endelse
            endif
         endif
         print, format='("File ",i5," of ",i5,a1,$)', $
          isdr+1, nsdr, string(13b)
      endfor

      if (NOT keyword_set(tai_first)) then begin
         print, 'No science/smear exposures for this plate', plate
         return
      endif
      ikeep = where(timearray GE tai_first AND timearray LE tai_last, nfile)
      if (nfile EQ 0) then begin
         print, 'No guider images during timestamps for this plate', plate
         return
      endif
      print, 'Trimming to ', nfile, ' guider frames for plate ', plate
      gfiles = gfiles[ikeep]
      datestring = datestring[ikeep]
      timearray = timearray[ikeep]
   end

   ;----------
   ; Read the dark image

   dark = mrdfits('~jeg/Data/guider/darks/30s/SUMDKBIN.FIT')


   titlestring = 'MJD ' + mjdstr
   if (keyword_set(plate)) then $
    titlestring = 'Plate ' + string(plate,format='(i4)') + ' ' + titlestring

   ifile = 0
   cc = 'F'
   lastfile = -1

   while (strupcase(cc) NE 'Q') do begin

      if (cc EQ 'F') then begin
         ifile = ifile + 1
         if (ifile GE nfile) then begin
            ifile = nfile - 1
            cc = 's'
         endif
         c1 = get_kbrd(0)
         if (keyword_set(c1)) then cc = c1
      endif else if (cc EQ 'B') then begin
         ifile = ifile - 1
         if (ifile LT 0) then begin
            ifile = 0
            cc = 's'
         endif
         c1 = get_kbrd(0)
         if (keyword_set(c1)) then cc = c1
      endif else begin
         print, 'Press b=back one frame'
         print, '      f=forward one frame'
         print, '      B=loop backward'
         print, '      F=loop forward'
         print, '      s=stop'
         print, '      q=quit (and enter interactive mode for this plot)'
         cc = get_kbrd(1)
         case cc of
            'b': ifile = (ifile - 1) > 0
            'f': ifile = (ifile + 1) < (nfile - 1)
            else:
         endcase
      endelse

      ; Only display this image if it is different from the last image displayed
      if (ifile NE lastfile) then begin
         img = mrdfits(gfiles[ifile], 0, /silent) + 32768.

         splog, ifile, nfile, fileandpath(gfiles[ifile]), utstring, $
          format='("File ",i4," of ",i4,": ",a,2x,a)'

         ; Dark-subtract
         if (total(size(dark,/dimens) - size(img,/dimens)) EQ 0) then begin
            res = linfit(dark[*,0:50], img[*,0:50])
            img = img - res[0] - res[1] * dark
         endif

         atv, img, _EXTRA=KeywordsForATV
         atvxyouts, 0, (size(img,/dimens))[1], titlestring, $
          color='red', charsize=4
         atvxyouts, 0, 0, datestring[ifile], $
          color='red', charsize=4

         lastfile = ifile
      endif else begin
         wait, 0.5 ; Don't sit around burning CPU cycles.
      endelse

   endwhile
end
;------------------------------------------------------------------------------
