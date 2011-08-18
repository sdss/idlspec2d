;+
; NAME:
;   apoheader
;
; PURPOSE:
;   Print the subset of interesting header keywords from the raw sdR files.
;
; CALLING SEQUENCE:
;   apoheader, expnum, [ mjd= ]
;
; INPUTS:
;   expnum     - Exposure number
;
; OPTIONAL INPUTS:
;   mjd        - Optionally specify the MJD for the subdirectory in which
;                to search for the file.  This will greatly speed things up.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   Print out all FITS header keywords of interest for exposure # 14728
;   (assuming that exposure is still on disk somewhere under $BOSS_SPECTRO_DATA):
;     IDL> apoheader, 14728
;
; BUGS:
;
; PROCEDURES CALLED:
;   fileandpath()
;   fits_wait
;   sdsshead()
;   sxpar()
;
; REVISION HISTORY:
;   24-Apr-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro apoheader, expnum, mjd=mjd

   if (n_params() LT 1 OR n_elements(expnum) NE 1) then begin
      print, 'Syntax - apoheader, expnum'
      return
   end

   quiet = !quiet
   !quiet = 1

   ;----------
   ; Set input directory for sdR files

   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   if (NOT keyword_set(rawdata_dir)) then $
    rawdata_dir = '/data/spectro'

   ;----------
   ; Sanity checks on EXPNUM

   expnum = long(expnum)
   if (expnum LE 0 OR expnum GT 99999999L) then begin
      print, 'EXPNUM must be a number between 1 and 99999999'
      !quiet = quiet
      return
   endif

   ;----------

   if (keyword_set(mjd)) then subdir = strtrim(string(mjd),2) $
    else subdir = '*'

   ;----------
   ; Read the headers for the 4 files.

   camname = ['b1','r1','b2','r2']
   ncam = n_elements(camname)
   filename = strarr(ncam)
   phdr = ptrarr(ncam)
   for icam=0, ncam-1 do begin
      fileroot = string(camname[icam], expnum, format='("sdR-",a2,"-",i8.8)')
      filename[icam] = (findfile(filepath(fileroot+'.fit*', $
       root_dir=rawdata_dir, subdir=subdir), count=ct))[0]
      if (ct GT 0) then begin
         qdone = fits_wait(filename[icam], deltat=1, tmax=1, /header_only)
         if (qdone) then begin
;            thishdr = headfits(filename[icam])
            thishdr = sdsshead(filename[icam], /do_lock)
            sxaddpar, thishdr, '(END-BEG', $
             sxpar(thishdr,'TAI-END') - sxpar(thishdr,'TAI-BEG')
            phdr[icam] = ptr_new(thishdr)
         endif
      endif
   endfor

   ;----------

   cardname = [ $
    'EXPOSURE' , $
    'CAMERAS'  , $
    'FLAVOR'   , $
    'MJD'      , $
    'PLATEID'  , $
    'NAME'     , $
    'EXPTIME'  , $
    '(END-BEG' , $
    'TAI-BEG'  , $
    'TAI-END'  , $
    'TAI'      , $
    'FFS'      , $
    'NE'       , $
    'HGCD'     , $
    'OBSCOMM'  , $
    'QUALITY'  , $
; These are informational keywords...
    'TILEID'   , $
    'CARTID'   , $
    'RA'       , $
    'DEC'      , $
    'RADEG'    , $
    'DECDEG'   , $
    'AIRTEMP'  ]

   print, fileandpath(filename), format='(10x,4(x,a16))'
   print, format='(10x,4(" ----------------"))'
   sval = strarr(ncam)
   for icard=0, n_elements(cardname)-1 do begin
      for icam=0, ncam-1 do begin
         if (keyword_set(phdr[icam])) then begin
            thisval = sxpar(*phdr[icam], cardname[icard])
            if (cardname[icard] EQ 'TAI' $
             OR cardname[icard] EQ 'TAI-BEG' $
             OR cardname[icard] EQ 'TAI-END') then $
             format='(e15.9)' $
            else format=''
            sval[icam] = strtrim(string(thisval,format=format))
         endif else begin
            sval[icam] = ''
         endelse
      endfor
      print, cardname[icard], sval, format='(a8,2x,4a17)'
   endfor

   !quiet = quiet

   return
end
;------------------------------------------------------------------------------
