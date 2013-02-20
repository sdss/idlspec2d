;+
; NAME:
;   sos_logs_concat
;
; PURPOSE:
;   Concatenate all the SOS log files into a single file
;
; CALLING SEQUENCE:
;   sos_logs_concat, [ mjd=, outfile= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - Search for all files in $SPECTROLOGS_DIR/MJD;
;                default to '5????'
;   outfile    - Output file name; default to 'logfile-all.fits'
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   All the SOS output files in $SPECTROLOG_DIR/$MJD/logfile-$MJD.fits
;   are concatenated together into a single file.  The HDUs are:
;     HDU #1: bias and dark frames
;     HDU #2: flats
;     HDU #3: arcs
;     HDU #4: science frames
;   The exact format of each data structure is the same as the first one read.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;  18-Feb-2013  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro sos_logs_concat, mjd=mjd1, outfile=outfile1

   spectrolog_dir = getenv('SPECTROLOG_DIR')

   if (keyword_set(mjd1)) then mjdstr = strtrim(mjd1,2) $
    else mjdstr = '5????'
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile = 'logfile-all.fits'

   mjdlist = ''
   for i=0, n_elements(mjdstr)-1 do begin
      spawn, '\ls -d '+spectrolog_dir + '/' + mjdstr[i], mjdlist1
      if (keyword_set(mjdlist1)) then $
       mjdlist = keyword_set(mjdlist) ? [mjdlist,mjdlist1] : mjdlist1
   endfor
   if (NOT keyword_set(mjdlist)) then begin
      print, 'No files found'
      return
   endif
   nmjd = n_elements(mjdlist)
   nhdu = 4

   ; Find if a logfile exists in each MJD, and count the number
   ; of entries in each HDU
   filelist = 0
   thislist1 = create_struct('FILENAME', '', 'NUM', lonarr(nhdu))
   for i=0L, nmjd-1L do begin
      mjdstr = string(fileandpath(mjdlist[i]),format='(i5.5)')
      thisname = filepath('logfile-'+mjdstr+'.fits', $
       root_dir=mjdlist[i])
      thisname = findfile(thisname, count=ct)
      if (ct GT 0) then begin
         thislist1.filename = thisname
         for j=0, nhdu-1 do $
          thislist1.num[j] = sxpar(headfits(thisname,exten=j+1), 'NAXIS2')
         if (NOT keyword_set(filelist)) then filelist = thislist1 $
          else filelist = [filelist, thislist1]
      endif
   endfor

   ; Read all data and copy into one data structure per HDU
   nfile = n_elements(filelist)
   count = lonarr(nhdu)
   outdat = ptrarr(nhdu)
   for i=0L, nfile-1L do begin
      print, 'Reading file ', i, ' of ', nfile, ': ' + filelist[i].filename
      for j=0, nhdu-1 do begin
         if (filelist[i].num[j] GT 0) then begin
            dat1 = mrdfits(filelist[i].filename,j+1)
            if (NOT keyword_set(outdat[j])) then begin
               blankdat = dat1[0]
               struct_assign, {junk:0}, blankdat
               outdat[j] = ptr_new( replicate(blankdat,total(filelist.num[j])) )
            endif
;            (*outdat[j])[count[j]+lindgen(filelist[i].num[j])] = dat1
            copy_struct_inx, dat1, *outdat[j], $
             index_to=count[j]+lindgen(filelist[i].num[j])
            count[j] += filelist[i].num[j]
         endif
      endfor
   endfor

   ; Write output file
   print, 'Writing file '+outfile
   mwrfits, 0, outfile, /create
   for j=0, nhdu-1 do $
    mwrfits, *outdat[j], outfile

   return
end
;------------------------------------------------------------------------------
