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
;                in the $BOSS_SPECTRO_DATA directory.
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
;   Display the horrendous guider images for plate 889 on MJD=52346:
;     IDL> guidermovie, mjd=52346, plate=889
;
; BUGS:
;
; DATA FILES:
;   $BOSS_SPECTRO_DATA/$MJD/guider/gimg*.fits.gz
;   $BOSS_SPECTRO_DATA/$MJD/sdR-b1-????????.fit.gz
;
; PROCEDURES CALLED:
;   atv
;   atvxyouts
;   djs_filepath()
;   fileandpath()
;   findopfile()
;   get_tai
;   jdcnv
;   mrdfits()
;   splog
;   sdsshead()
;   sxpar()
;
; REVISION HISTORY:
;   13-May-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro guidermovie, mjd=mjd, plate=plate, expnum=expnum, _EXTRA=KeywordsForATV

   common gmovie, com_mjd, com_gfiles, com_datestring, com_timearray, $
    com_plate, com_mjdplate, com_tai_first, com_tai_last

   quiet = !quiet
   !quiet = 1
   if (NOT keyword_set(com_mjd)) then com_mjd = 0
   if (NOT keyword_set(com_plate)) then com_plate = 0
   if (NOT keyword_set(com_mjdplate)) then com_mjdplate = 0

   ;----------
   ; If MJD is not specified, then find the most recent MJD for output files

   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   if (NOT keyword_set(rawdata_dir)) then $
    rawdata_dir = '/data/spectro'

   if (NOT keyword_set(mjd)) then begin
      mjdlist = get_mjd_dir(rawdata_dir, mjstart=1, mjend=99999, mjd='?????')
      mjd = (reverse(mjdlist[sort(mjdlist)]))[0]
      splog, 'Selecting MJD=', mjd, ' (override this with MJD keyword)'
   endif
   mjdstr = string(mjd, format='(i5.5)')

   if (com_mjd EQ mjd) then begin
      ;----------
      ; If the last call to this procedure was for the same MJD, then use
      ; the results cached in these common variables.

      print, 'Using cached values for guider image header info'
      mjd = com_mjd
      gfiles = com_gfiles
      datestring = com_datestring
      timearray = com_timearray
      nfile = n_elements(gfiles)
   endif else begin
      ;----------
      ; Get file list for guider images (g-zipped only)

      gfiles = findfile( djs_filepath('gimg*fits.gz', root_dir=rawdata_dir, $
       subdirectory=[string(mjd,format='(i5.5)'),'guider']), count=nfile )
      if (nfile EQ 0) then begin
         splog, 'No files found
         return
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

         ; The old format for the date was, for ex, "Sun, Apr 1, 2001".
         ; The new format for the date is, for ex, "2002/07/08".
         utdate = sxpar(hdr,'UTDATE')
         utdate = repstr(utdate, ',', ' ') ; Replace commas w/whitespace
         ww = strsplit(utdate, ', ', /extract)
         if (n_elements(ww) GE 3) then begin
            year = long(ww[3])
            month = (where(ww[1] EQ monthname))[0] > 0
            date = long(ww[2])
         endif else begin
            ww = strsplit(utdate, '/', /extract)
            year = long(ww[0])
            month = long(ww[1])
            date = long(ww[2])
         endelse

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

      com_mjd = mjd
      com_gfiles = gfiles
      com_datestring = datestring
      com_timearray = timearray
   endelse

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

      if (plate EQ com_plate AND mjd EQ com_mjdplate) then begin
         ;----------
         ; If the last call to this procedure was for the same PLATE and MJD,
         ; then use the results cached in these common variables.

         print, 'Using cached values for sdR header info'
         tai_first = com_tai_first
         tai_last = com_tai_last
      endif else begin
         ;----------
         ; Set input directory for sdR files

         rawdata_dir = getenv('BOSS_SPECTRO_DATA')
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

         tai_first = 0
         tai_last = 0
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

         com_plate = plate
         com_mjdplate = mjd
         com_tai_first = tai_first
         com_tai_last = tai_last
      endelse

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
   ; Read the dark image (select only from the 30-sec darks)

   darkdir = filepath('',root_dir=getenv('IDLSPEC2D_DIR'), subdir='etc')
   darkfile = findopfile('dark-30-*.fits*', mjd[0], darkdir)
   if (keyword_set(darkfile)) then $
    dark = mrdfits(filepath(darkfile, root_dir=darkdir), /fscale)

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
            cc = ' '
         endif
      endif else if (cc EQ 'B') then begin
         ifile = ifile - 1
         if (ifile LT 0) then begin
            ifile = 0
            cc = ' '
         endif
      endif else begin
         print, 'Press b=back one frame'
         print, '      f=forward one frame'
         print, '      B=loop backward'
         print, '      F=loop forward'
         print, '      s=select frame index'
         print, '      q=quit (and enter interactive mode for this plot)'
         print, '      <any other key>=stop'
         case cc of
            'b': ifile = (ifile - 1) > 0
            'f': ifile = (ifile + 1) < (nfile - 1)
            's': begin
               read, ifile, prompt='Enter frame index number: '
               ifile = (ifile > 0) < (nfile - 1)
               cc = ' '
               end
            else:
         endcase
         cc = ' ' ; Clear this command (since it's not a looping command)
      endelse

      ; Only display this image if it is different from the last image displayed
      if (ifile NE lastfile) then begin
         img = mrdfits(gfiles[ifile], 0, hdr, /silent) + 32768.
         exptime = sxpar(hdr, 'EXPTIME')

         splog, ifile+1, nfile, fileandpath(gfiles[ifile]), utstring, $
          format='("File ",i4," of ",i4,": ",a,2x,a)'

         ; Dark-subtract, scaling the dark current from the 30-sec dark
         qsamesize = total(size(img,/dimens) EQ size(dark,/dimens)) EQ 2
         if (keyword_set(exptime)*qsamesize) then $
          img -= (exptime/30.) * dark

         atv, img, _EXTRA=KeywordsForATV
         atvxyouts, 0, (size(img,/dimens))[1], titlestring, $
          color='red', charsize=4
         atvxyouts, 0, 0, datestring[ifile], $
          color='red', charsize=4

         lastfile = ifile
      endif else begin
         wait, 0.5 ; Don't sit around burning CPU cycles.
      endelse

      ;----------
      ; If we're looping, then just read whatever is in the keyboard buffer,
      ; otherwise sit and wait for a keystroke.

      if (cc EQ 'F' OR cc EQ 'B') then begin
         c1 = get_kbrd(0)
         if (keyword_set(c1)) then cc = c1
      endif else begin
         cc = get_kbrd(1) ; Sit and wait for keystroke
      endelse

   endwhile
end
;------------------------------------------------------------------------------
