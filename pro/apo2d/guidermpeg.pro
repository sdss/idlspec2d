;+
; NAME:
;   guidermpeg
;
; PURPOSE:
;   Create an MPEG guider movie.
;   Plot the guider images for an entire night.
;
; CALLING SEQUENCE:
;   guidermovie, [ mjd=, expnum=, _EXTRA=KeywordsForATV ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - MJD number; if not set, then select the most recent MJD
;                in the $BOSS_SPECTRO_DATA directory.
;   expnum     - Exposure numbers as an array
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   Make a movie of the horrible guider images while the slit-head
;   became unlatched while observing plate 889 on MJD 52346:
;     IDL> guidermpeg, mjd=52346, expnum=805+lindgen(25)
;
; BUGS:
;
; DATA FILES:
;   $BOSS_SPECTRO_DATA/$MJD/guider/gimg*.fits.gz
;
; PROCEDURES CALLED:
;   djs_filepath()
;   fileandpath()
;   findopfile()
;   mrdfits()
;   splog
;   sxpar()
;
; REVISION HISTORY:
;   13-May-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro guidermpeg, mjd=mjd, expnum=expnum

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

   ;----------
   ; Get file list for guider images (g-zipped only)

   gfiles = findfile( djs_filepath('gimg*fits.gz', root_dir=rawdata_dir, $
    subdirectory=[string(mjd,format='(i5.5)'),'guider']), count=nfile )
   if (nfile EQ 0) then begin
      splog, 'No files found
      return
   endif

   ;----------
   ; Trim to only specified exposure numbers

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
   ; Read the dark image (select only from the 30-sec darks)

   darkdir = filepath('',root_dir=getenv('IDLSPEC2D_DIR'), subdir='etc')
   darkfile = findopfile('dark-30-*.fits*', mjd[0], darkdir)
   if (keyword_set(darkfile)) then $
    dark = mrdfits(filepath(darkfile, root_dir=darkdir), /fscale)

   ;----------
   ; Constuct the output byte array

   xsize = 320L
   ysize = 240L
   bytarr = bytarr(xsize,ysize,nfile)
   npix = xsize * ysize
   xoff = 40
   yoff = 16

   ;----------
   ; Loop through all images

   for ifile=0, nfile-1 do begin

      img = mrdfits(gfiles[ifile], 0, hdr, /silent) + 32768
      exptime = sxpar(hdr, 'EXPTIME')

      ; Dark-subtract, scaling the dark current from the 30-sec dark
      qsamesize = total(size(img,/dimens) EQ size(dark,/dimens)) EQ 2
      if (keyword_set(exptime)*qsamesize) then $
       img -= (exptime/30.) * dark

      ; Trim image to XSIZE,YSIZE
      img = img[xoff:xoff+xsize-1,yoff:yoff+ysize-1]

      ; Rescale to a byte image
      isort = sort(img)
      min = img[isort[long(0.25*npix)]] ; 25-th percentile
      max = img[isort[long(0.995*npix)]] ; 99.5-th percentile
      bytarr[*,*,ifile] = bytscl(img, min=min, max=max)
   endfor

   ppmtompeg, bytarr, 'junk.mpeg', tmpdir='./'

end
;------------------------------------------------------------------------------
