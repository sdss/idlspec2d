;+
; NAME:
;   spflatgen
;
; PURPOSE:
;   Wrapper for SPFLATTEN2 for generating pixel-to-pixel flat-fields.
;
; CALLING SEQUENCE:
;   spflatgen, [ mjd=, indir=, outdir=, docam= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - If INDIR not set, then look for files in $RAWDATA_DIR/MJD.
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
;   Four FITS files are produced, one for each camera:
;     pixflat-MJD-CAMERA.fits
;   where MJD is the 5-digit modified Julian date, and CAMERA
;   is 'b1', 'b2', 'r1', and 'r2'.
;
; EXAMPLES:
;   The following nights look like they contain a flat sequence:
;     51686 51781 51809 51852 51893 51950 51978 52010 52038 52069
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
pro spflatgen, mjd=mjd, indir=indir, outdir=outdir, docam=docam

   if (keyword_set(docam)) then camnames = docam $
    else camnames = ['b1', 'b2', 'r1', 'r2']
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
   exposure = lonarr(nfile)

   print, 'Reading FITS headers...', format='(a,$)'
   for ifile=0, nfile-1 do begin
      hdr = sdsshead(files[ifile])
      cameras[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
      obscomm[ifile] = strtrim(sxpar(hdr, 'OBSCOMM'),2)
      mjdarr[ifile] = sxpar(hdr, 'MJD')
      exposure[ifile] = sxpar(hdr, 'EXPOSURE')
      print, '.', format='(a,$)'
   endfor
   print

   for icam=0, ncam-1 do begin
      ; Select flats and arcs for this camera
      iflats = where(cameras EQ camnames[icam] $
       AND obscomm EQ '{dithered flats-flat}', nflat)
      iarcs = where(cameras EQ camnames[icam] $
       AND obscomm EQ '{dithered flats-arc}', narc)

      ; Re-sort each list by exposure number
      iflats = iflats[ sort(exposure[iflats]) ]
      iarcs = iarcs[ sort(exposure[iarcs]) ]

      ; Trim the flats to only the list of contiguous flats,
      ; according to their exposure numbers
      iflats = iflats[ where(exposure[iflats] - exposure[iflats[0]] $
       - lindgen(nflat) EQ 0, nflat) ]

      if (nflat GE 7 AND narc GE 1) then begin
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
