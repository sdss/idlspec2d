;+
; NAME:
;   spbiasgen
;
; PURPOSE:
;   Routine to generate mean biases for a night.
;
; CALLING SEQUENCE:
;   spbiasgen, [ mjd=, expnum=, expstart=, expend=, timesep=, $
;    indir=, outdir=, docam=, sigrej=, maxiter=, /noproc ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - If INDIR not set, then look for files in $BOSS_SPECTRO_DATA/MJD.
;   expnum     - If set, then use these exposure numbers
;   expstart   - If set, then only use exposure numbers >= EXPSTART
;   expend     - If set, then only use exposure numbers >= EXPEND
;   timesep    - Discard any bias that isn't within TIMESEP after another bias;
;                default to 300 sec.  The time is obtained from TAI-BEG in the
;                FITS header.  Note that this always discards the first bias in
;                any sequence.
;   indir      - Look for input files in this directory; default to current
;                directory if neither MJD or INDIR are set.
;   outdir     - Output directory; default to same as INDIR.
;   docam      - Camera names; default to all cameras: ['b1', 'b2', 'r1', 'r2']
;   sigrej     - Sigma rejection level; default to 1, 1, 1.1, 1.3, 1.6 or 1.9
;                for 1,2,3,4,5 or 6 flats.  For more then 6 flats, default
;                to 2.0.
;   maxiter    - Number of rejection iterations; default to 3.
;   noproc     - If set, then use MRDFITS to process images, otherwise use
;                SDSSPROC and remove the overscan regions
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine looks for biases in a given night that are logged as such.
;   Discard the first bias in any sequence, and generate a mean bias from
;   all the others for each camera.  There must be at least 3 frames.
;   A sigma-clipped mean is computed for each pixel.
;
;   Trigger a failure if there are more than 25 biases for a given camera,
;   since that is probably too many to load into memory.
;
;   Four FITS files are produced, one for each camera:
;     pixbias-MJD-CAMERA.fits
;   where MJD is the 5-digit modified Julian date, and CAMERA
;   is 'b1', 'b2', 'r1', and 'r2'.
;
;   The header contains the number of files used in each camera (NEXP),
;   and an identifier for each of those files (EXPID*).  Those identifiers
;   contain the camera name, MJD, and exposure number, all dash-separated.
;
; EXAMPLES:
;   The following SDSS-I nights probably contain a bias sequence:
;     51686 51781 51809 51852 51893 51950 51978 52010 52038 52069
;     52245 52276 52305 52333 52363 52392 52423 52454 52514 52551
;     52573 52602 52633 52655 52689 52718 52747 (<- should declare some bad!)
;     52808 52868 52899 52956 52983 53014 53044 53073 53100 53132
;     53162 53220 53257 ...???
;     53610 53633
;   (Those are the nights with a pixel flat sequence.)  There are other
;   nights with many biases that may be useful, but may have been done
;   for other testing purposes:
;     52084 52085 52121 52229 52294 52550 52572 52601 52629
;     52657 52660 52661 52717 52807 52863 52864 52865 52926
;   Generate one of these sets of biases with:
;     spbiasgen, mjd=51893, expstart=7607, expend=7616, outdir='.'
;
;   For BOSS:
;     spbiasgen,mjd=55268,expstart=111620,expend=111645,outdir='.'
;     spbiasgen,mjd=55271,expstart=111881,expend=111905,outdir='.'
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_avsigclip()
;   djs_filepath()
;   fileandpath()
;   get_tai
;   sdsshead()
;   sdssproc
;   splog
;   sxaddpar
;   sxpar()
;   writefits
;
; INTERNAL SUPPORT ROUTINES:
;   spbiasgen1
;
; REVISION HISTORY:
;   30-Aug-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro spbiasgen1, files, outfile=outfile, outdir=outdir, $
 sigrej=sigrej, maxiter=maxiter, noproc=noproc

   nfile = n_elements(files)

   if (NOT keyword_set(sigrej)) then begin
      if (nfile LE 2) then sigrej = 1.0 $ ; Irrelevant for only 1 or 2 files
       else if (nfile EQ 3) then sigrej = 1.1 $
       else if (nfile EQ 4) then sigrej = 1.3 $
       else if (nfile EQ 5) then sigrej = 1.6 $
       else if (nfile EQ 6) then sigrej = 1.9 $
       else sigrej = 2.0
   endif
   if (NOT keyword_set(maxiter)) then maxiter = 3

   for ifile=0, nfile-1 do begin
      splog, 'Reading file #', ifile+1, ' of ', nfile
      if (keyword_set(noproc)) then begin
         thisimg = mrdfits(files[ifile], 0, hdr, /fscale)
      endif else begin
         sdssproc, files[ifile], thisimg, thisivar, hdr=hdr, $
          bcfile='opBC-empty.par'
      endelse
      if (ifile EQ 0) then begin
         hdr0 = hdr
         imgarr = make_array(dimension=[size(thisimg,/dimens),nfile], /float)
         inmask = make_array(dimension=[size(thisimg,/dimens),nfile], /byte)
         sxaddpar, hdr0, 'NEXP', nfile, 'Number of exposures in this file', $
          before='EXPTIME'
      endif
      imgarr[*,*,ifile] = thisimg
      if (keyword_set(thisivar)) then $
       inmask[*,*,ifile] = thisivar LE 0

      sxaddpar, hdr0, string(ifile+1,format='("EXPID",i2.2)'), $
       string( sxpar(hdr,'CAMERAS'), sxpar(hdr,'EXPOSURE'), $
        format='(a2,"-",i8.8)'), $
        'ID string for exposure ' + strtrim(ifile+1,2), before='EXPTIME'
   endfor

   mnimg = djs_avsigclip(imgarr, sigrej=sigrej, maxiter=maxiter, $
    inmask=inmask)

   writefits, djs_filepath(outfile, root_dir=outdir), mnimg, hdr0

   return
end
;------------------------------------------------------------------------------
pro spbiasgen, mjd=mjd, expnum=expnum, expstart=expstart, expend=expend, $
 timesep=timesep, indir=indir, outdir=outdir, docam=docam, noproc=noproc

   if (NOT keyword_set(timesep)) then timesep = 300
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
   flavor = strarr(nfile)
   mjdarr = lonarr(nfile)
   exposure = lonarr(nfile)
   taitime = dblarr(nfile)

   print, 'Reading FITS headers...', format='(a,$)'
   for ifile=0, nfile-1 do begin
      hdr = sdsshead(files[ifile])
      cameras[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
      flavor[ifile] = strtrim(sxpar(hdr, 'FLAVOR'),2)
      mjdarr[ifile] = sxpar(hdr, 'MJD')
      exposure[ifile] = sxpar(hdr, 'EXPOSURE')
      get_tai, hdr, tai_beg, tai_mid, tai_end
      taitime[ifile] = tai_beg
      print, '.', format='(a,$)'
   endfor
   print

   for icam=0, ncam-1 do begin
      ; Select biases for this camera
      ibias = where(cameras EQ camnames[icam] $
       AND (strtrim(flavor,2) EQ 'bias' OR strtrim(flavor,2) EQ 'dark'), nbias)

      ; Discard any bias that isn't within TIMESEP seconds after
      ; another bias (this always discards the first one in any sequence).
      jbias = ibias ; This is a copy of those indices, putting -1 for bad ones
      for i=0, nbias-1 do begin
         junk = where(taitime[ibias[i]] GE taitime[ibias] $
                  AND taitime[ibias[i]] LE taitime[ibias] + timesep, nmatch)
         if (nmatch LE 1) then jbias[i] = -1 ; One match is the frame itself.
      endfor
      j = where(jbias NE -1, nbias)
      if (nbias GT 25) then begin
         splog, 'Too many biases (' + string(nbias) + ') to load into memory'
         splog, 'Skipping pixbias generation for camera ' + camnames[icam]
      endif else if (nbias GE 3) then begin
         jbias = jbias[j]

         pixbiasname = 'pixbias-' + string(mjdarr[jbias[0]],format='(i5.5)') $
          + '-' + camnames[icam] + '.fits'
         splog, 'Generating pixel bias ' + pixbiasname
         splog, 'Output directory ' + outdir
         spbiasgen1, files[jbias], outfile=pixbiasname, outdir=outdir, $
          sigrej=sigrej, maxiter=maxiter, noproc=noproc

      endif else begin
         splog, 'Expected at least 3 biases (not including 1st), got ' $
          + strtrim(string(nbias),2)
         splog, 'Skipping pixbias generation for camera ' + camnames[icam]
      endelse
   endfor

   return
end
;------------------------------------------------------------------------------
