;+
; NAME:
;   spflatgen
;
; PURPOSE:
;   Wrapper for SPFLATTEN2 for generating pixel-to-pixel flat-fields.
;
; CALLING SEQUENCE:
;   spflatgen, [ mjd=, indir=, outdir= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - If INDIR not set, then look for files in $RAWDATA_DIR/MJD.
;   indir      - Look for input files in this directory; default to current
;                directory if neither MJD or INDIR are set.
;   outdir     - Output directory; default to same as INDIR.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine looks for flats and arcs in a given night that are logged
;   (according to the headers) to be spectroscopic pixel flats.  We expect
;   exactly 7 flats in a sequence plus one or more arcs.  If there are not
;   exactly 7 flats and at least one arc for a given camera, then no pixel
;   flats are generated for that camera.
;
;   Four FITS files are produced, one for each camera:
;     pixflat-MJD-CAMERA.fits
;   where MJD is the 5-digit modified Julian date, and CAMERA
;   is 'b1', 'b2', 'r1', and 'r2'.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_filepath()
;   sdssproc
;   spflatten2
;   splog
;   sxpar()
;
; REVISION HISTORY:
;   06-Jul-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro spflatgen, mjd=mjd, indir=indir, outdir=outdir

   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   if (keyword_set(mjd) AND NOT keyword_set(indir)) then begin
      indir = filepath('', root_dir=getenv('RAWDATA_DIR'), $
       subdirectory=string(mjd,format='(i5.5)'))
   endif
   if (keyword_set(indir) AND NOT keyword_set(outdir)) then $
    outdir = indir

   files = findfile(djs_filepath('sdR-*.fit*',root_dir=indir), count=nfile)
   if (nfile EQ 0) then begin
      splog, 'No files found.'
      return
   endif

   cameras = strarr(nfile)
   obscomm = strarr(nfile)
   mjdarr = lonarr(nfile)

   print, 'Reading FITS headers...', format='(a,$)'
   for ifile=0, nfile-1 do begin
      sdssproc, files[ifile], hdr=hdr
      cameras[ifile] = sxpar(hdr, 'CAMERAS')
      obscomm[ifile] = sxpar(hdr, 'OBSCOMM')
      mjdarr[ifile] = sxpar(hdr, 'MJD')
      print, '.', format='(a,$)'
   endfor
   print

   for icam=0, ncam-1 do begin
      iflats = where(cameras EQ camnames[icam] $
       AND obscomm EQ 'dithered flats-flat', nflat)
      iarcs = where(cameras EQ camnames[icam] $
       AND obscomm EQ 'dithered flats-arc', narc)
      if (nflat EQ 7 AND narc GE 1) then begin
         pixflatname = 'pixflat-' + camnames[icam] $
          + '-' + string(mjdarr[iflats[0]],format='(i5.5)') + '.fits'
         spflatten2, files[iflats[4]], files[iarcs[0]], files[iflats], $
          pixflatname
      endif else begin
         splog, 'Expected 7 flats + some arcs, got ' $
          + strtrim(string(nflat),2) + ' + ' + strtrim(string(narc),2)
         splog, 'Skipping pixflat generation for camera ' + camnames[icam]
      endelse
   endfor

   return
end
;------------------------------------------------------------------------------
