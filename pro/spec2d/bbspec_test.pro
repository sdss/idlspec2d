;+
; NAME:
;   bbspec_test
;
; PURPOSE:
;   Run bbspec 2D extraction code as an after-burner to existing reductions
;
; CALLING SEQUENCE:
;   bbspec_test, scifile, [ outfile=, /clobber, /batch, _EXTRA= ]
;
; INPUTS:
;   scifile    - spFrame file name
;
; OPTIONAL INPUTS:
;   outfile    - Output FITS file with 2D image model; default to
;                'ymodel-test.fits'
;   clobber    - If set, then clobber any existing PSF and re-generate it
;   batch      - If set, then batch each bundle of 20 fibers to a different
;                PBS job; FRANGE keyword cannot also be set
;   _EXTRA     - Keywords for BBSPEC_EXTRACT, such as FRANGE,YRANGE
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine calls the bbspec PSF-construction and 2D extraction
;   code as an afterburner to the pipeline, assuming that there already
;   exists an spFrame, spArc, spFlat.  Existing spBasisPSF files are not
;   over-written unless /CLOBBER is set.
;
;   If the spBasisPSF file cotains only the first N fibers, then only
;   those N fibers are extracted.  If /BATCH is set, then N must be
;   a multiple of 20 since the code is then called in explicit batches
;   of 20 fibers.
;
;   The output file contains 3 HDUs with the 2-D model image,
;   extracted fluxes, and extracted inverse variances.
;
; EXAMPLES:
;
; BUGS:
;
; REVISION HISTORY:
;   24-May-2011  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function bbspec_test_batch, scifile, _EXTRA=Extra

   fq = "'"
   idlcmd = 'bbspec_test,'+fq+scifile+fq
   jobid = bbspec_batch(idlcmd, _EXTRA=Extra)

   return, jobid
end
;------------------------------------------------------------------------------
pro bbspec_test, scifile, outfile=outfile1, clobber=clobber, batch=batch, $
 _EXTRA=Extra

   t0 = systime(1)

   if (n_params() NE 1) then $
    message, 'Wrong number of parameters'
   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   if (keyword_set(batch) AND keyword_set(Extra)) then begin
      if (tag_exist(Extra,'FRANGE')) then $
       message, 'Cannot specify both BATCH and FRANGE'
   endif

   hdr = headfits(scifile)
   if (NOT keyword_set(hdr)) then $
    message, 'Error reading file '+scifile
   mjdstr = string(sxpar(hdr, 'MJD'),format='(i5.5)')
   indir = concat_dir(rawdata_dir, mjdstr)

   rawfile = 'sdR-'+strmid(sxpar(hdr, 'FILENAME'),4,11)+'.fit*'
   splog, 'Reading raw science image'
   sdssproc, rawfile, image, invvar, indir=indir, $
    /applybias, /applypixflat, /applycrosstalk, minflat=0.8, maxflat=1.2
   if (NOT keyword_set(image)) then $
    message, 'Unable to find raw science frame'

   arcstr = strmid(sxpar(hdr, 'ARCFILE'),4,11)
   flatstr = strmid(sxpar(hdr, 'FLATFILE'),4,11)
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile = 'ymodel-test.fits'

   arcfile = (findfile('spArc-'+arcstr+'.fits*', count=ct))[0]
   if (ct EQ 0) then $
    message, 'Unable to find spArc file'
   flatfile = (findfile('spFlat-'+flatstr+'.fits*', count=ct))[0]
   if (ct EQ 0) then $
    message, 'Unable to find spFlat file'

   basisfile = 'spBasisPSF-*-' + strmid(sxpar(hdr, 'ARCFILE'),4,11) + '.fits'
   junk = (findfile(basisfile+'*', count=ct))[0]
   splog, 'Selecting file '+basisfile
   if (ct EQ 0 OR keyword_set(clobber)) then begin
      splog, 'Generating sdProc file for arc'
      arcname = 'sdR-'+arcstr+'.fit'
      sdssproc, arcname, indir=indir, /outfile, $
       /applybias, /applypixflat, /applycrosstalk
      splog, 'Generating spBasisPSF file '+basisfile
      pyfile = djs_filepath('make-my-psf.py', root_dir=getenv('BBSPEC_DIR'), $
       subdir='examples')
      cmd = 'python '+pyfile+' '+arcstr+' '+flatstr
      spawn, cmd, res, errcode
      if (keyword_set(errcode)) then begin
         splog, errcode
         message, 'Error calling '+cmd
      endif
   endif else begin
      splog, 'Use existing PSF file '+basisfile
   endelse

   splog, 'Reading existing traceset '+flatfile
   xset = mrdfits(flatfile, 1)
   traceset2xy, xset, xx, ximg

   splog, 'Running 2D extraction'
   if (NOT keyword_set(batch)) then begin
      bbspec_extract, image, invvar, flux, fluxivar, basisfile=basisfile, $
       ximg=ximg, ymodel=bb_ymodel, _EXTRA=Extra
   endif else begin
      nfiber = sxpar(headfits(basisfile),'NAXIS2')
      njob = nfiber / 20
      jobid = lonarr(njob)
      splog, 'Batching ', njob, ' jobs for ', nfiber, ' fibers in PSF file'
      tmproot = 'tmp-'+strmid(scifile,8,11)+'-' $
       +string(lindgen(njob),format='(i2.2)')
      tmpoutfile = tmproot+'.fits'
      for i=0, njob-1 do $
       if (file_test(tmpoutfile[i])) then file_delete, tmpoutfile[i]
      for i=0, njob-1 do $
       jobid[i] = bbspec_test_batch(scifile, _EXTRA=Extra, $
        frange=[i*20,i*20+19], outfile=tmpoutfile[i], tmproot=tmproot[i]+'-')
      bbspec_batch_wait, jobid
      bb_ymodel = 0
      flux = 0
      fluxivar = 0
      for i=0, njob-1 do bb_ymodel += mrdfits(tmpoutfile[i],0)
      for i=0, njob-1 do flux += mrdfits(tmpoutfile[i],1)
      for i=0, njob-1 do fluxivar += mrdfits(tmpoutfile[i],2)
   endelse

   splog, 'Writing file '+outfile
   mwrfits, bb_ymodel, outfile, /create
   mwrfits, flux, outfile
   mwrfits, fluxivar, outfile

   splog, 'Elapsed time = ', systime(1)-t0, ' sec'

   return
end
;------------------------------------------------------------------------------
