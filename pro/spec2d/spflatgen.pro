;+
; NAME:
;   spflatgen
;
; PURPOSE:
;   Wrapper for SPFLATTEN2 for generating pixel-to-pixel flat-fields.
;
; CALLING SEQUENCE:
;   spflatgen, [ mjd=, expnum=, expstart=, expend=, $
;    indir=, outdir=, docam= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - If INDIR not set, then look for files in $BOSS_SPECTRO_DATA/MJD.
;   expnum     - If set, then use these exposure numbers
;   expstart   - If set, then only use exposure numbers >= EXPSTART
;   expend     - If set, then only use exposure numbers >= EXPEND
;   indir      - Look for input files in this directory; default to current
;                directory if neither MJD or INDIR are set.
;   outdir     - Output directory; default to same as INDIR.
;   docam      - Camera names; default to all cameras: ['b1', 'b2', 'r1', 'r2']
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine looks for flats and arcs in a given night that are logged
;   (according to the headers) to be spectroscopic pixel flats.  We expect
;   at least 7 flats in a sequence plus one or more arcs.  If this is not
;   true for a given camera, then no pixel flats are generated for that camera.
;
;   Trigger a failure if there are more than 25 biases for a given camera,
;   since that is probably too many to load into memory.
;
;   Four FITS files are produced, one for each camera:
;     pixflat-MJD-CAMERA.fits
;   where MJD is the 5-digit modified Julian date, and CAMERA
;   is 'b1', 'b2', 'r1', and 'r2'.
;
;   The header contains the number of files used in each camera (NEXP),
;   and an identifier for each of those files (EXPID*).  Those identifiers
;   contain the camera name, MJD, flat exposure number, and arc exposure 
;   number, all dash-separated.
;
; EXAMPLES:
;   The following nights look like they contain a flat sequence:
;     51686 51781 51809 51852 51893 51950 51978 52010 52038 52069
;     52161 62193 52216 52245 52276 52304 52333 52363 52392 52423
;     52454 52514 52551 52573 52602 52633 52655 52689 52719 52747
;     52777 52808 52868 52899 52956 52983 53014 53044 53073 53100
;     53132 53162 53220 53257 ...???
;     53610 53633
;     53649 (15 dithers)
;     53651 (15 dithers)
;   Generate one of these sets of flats with:
;     spflatgen, mjd=51781, outdir='.'
;
;   The initial set of pixel flats were generated using six of the
;   flats taken on MJD=51441 from three different plates: exposure
;   numbers 1351,1352,1362,1363,1372,1373 with exp #1350 as the
;   fiducial flat and #1347 as the arc.
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_filepath()
;   sdsshead()
;   spflatten2
;   splog
;   sxpar()
;
; REVISION HISTORY:
;   06-Jul-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro spflatgen, mjd=mjd, expnum=expnum, expstart=expstart, expend=expend, $
 indir=indir, outdir=outdir, docam=docam

   if (keyword_set(docam)) then camnames = docam $
    else camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   ;----------
   ; Find all file names in the directory corresponding to this MJD

   if (keyword_set(mjd) AND NOT keyword_set(indir)) then begin
      indir = filepath('', root_dir=getenv('BOSS_SPECTRO_DATA'), $
       subdirectory=string(mjd,format='(i5.5)'))
   endif
   if (keyword_set(indir) AND NOT keyword_set(outdir)) then $
    outdir = indir

   files = findfile(djs_filepath('sdR-*.fit*',root_dir=indir), count=nfile)
   if (nfile EQ 0) then begin
      splog, 'No files found.'
      return
   endif

   ;----------
   ; Trim to exposure numbers specified by EXPNUM,EXPSTART,EXPEND

   if (keyword_set(expnum) OR keyword_set(expstart) OR keyword_set(expend)) $
    then begin
      qkeep = bytarr(nfile) + 1b
      exposure = long( strmid(fileandpath(files),7,8) )
      if (keyword_set(expnum)) then begin
         for ifile=0, nfile-1 do $
          qkeep[ifile] = (total(exposure[ifile] EQ long(expnum)) NE 0)
      endif
      if (keyword_set(expstart)) then $
       qkeep = qkeep AND (exposure GE long(expstart))
      if (keyword_set(expend)) then $
       qkeep = qkeep AND (exposure LE long(expend))
      ikeep = where(qkeep, nfile)
      if (nfile EQ 0) then begin
         splog, 'No files found matching EXPNUM,EXPSTART,EXPEND.'
         return
      endif
      splog, 'Trimming to ', nfile, ' files specified by EXPNUM,EXPSTART,EXPEND.'
      files = files[ikeep]
   endif

   ;----------
   ; Parse the FITS headers

   cameras = strarr(nfile)
   obscomm = strarr(nfile)
   quality = strarr(nfile)
   mjdarr = lonarr(nfile)
   exposure = lonarr(nfile)

   print, 'Reading FITS headers...', format='(a,$)'
   for ifile=0, nfile-1 do begin
      hdr = sdsshead(files[ifile])
      cameras[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
      obscomm[ifile] = strtrim(sxpar(hdr, 'OBSCOMM'),2)
      quality[ifile] = strtrim(sxpar(hdr, 'QUALITY'),2)
      mjdarr[ifile] = sxpar(hdr, 'MJD')
      exposure[ifile] = sxpar(hdr, 'EXPOSURE')
      print, '.', format='(a,$)'
   endfor
   print

   for icam=0, ncam-1 do begin
      ; Select flats and arcs for this camera
      iflats = where(cameras EQ camnames[icam] $
       AND quality NE 'bad' $
       AND obscomm EQ '{dithered flats-flat}', nflat)
      iarcs = where(cameras EQ camnames[icam] $
       AND quality NE 'bad' $
       AND obscomm EQ '{dithered flats-arc}', narc)

      ; Re-sort each list by exposure number
      if (nflat GT 0) then $
       iflats = iflats[ sort(exposure[iflats]) ]
      if (narc GT 0) then $
       iarcs = iarcs[ sort(exposure[iarcs]) ]

      ; Trim the flats to only the list of contiguous flats,
      ; according to their exposure numbers
      if (nflat GT 0) then $
       iflats = iflats[ where(exposure[iflats] - exposure[iflats[0]] $
        - lindgen(nflat) EQ 0, nflat) ]

      if (nflat GT 25) then begin
         splog, 'Too many flats (' + string(nbias) + ') to load into memory'
         splog, 'Skipping pixflat generation for camera ' + camnames[icam]
      endif else if (nflat GE 7 AND narc GE 1) then begin
         pixflatname = 'pixflat-' + string(mjdarr[iflats[0]],format='(i5.5)') $
          + '-' + camnames[icam] + '.fits'
         splog, 'Generating pixel flat ' + pixflatname
         splog, 'Output directory ' + outdir
         spflatten2, files[iflats[(nflat-1)/2]], files[iarcs[0]], files[iflats], $
          outfile=pixflatname, outdir=outdir
      endif else begin
         splog, 'Expected sequence of at least 7 flats + 1 arc, got ' $
          + strtrim(string(nflat),2) + ' + ' + strtrim(string(narc),2)
         splog, 'Skipping pixflat generation for camera ' + camnames[icam]
      endelse
   endfor

   return
end
;------------------------------------------------------------------------------
