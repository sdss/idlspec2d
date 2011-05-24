;+
; NAME:
;   bbspec_test
;
; PURPOSE:;
;   Run bbspec 2D extraction code as an after-burner to existing reductions
;
; CALLING SEQUENCE:
;   bbspec_test, scifile, [ outfile=, /clobber, _EXTRA= ]
;
; INPUTS:
;   scifile    - spFrame file name
;
; OPTIONAL INPUTS:
;   outfile    - Output FITS file with 2D image model; default to
;                'ymodel-test.fits'
;   clobber    - If set, then clobber any existing PSF and re-generate it
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
pro bbspec_test, scifile, outfile=outfile1, clobber=clobber, _EXTRA=Extra

   if (n_params() NE 1) then $
    message, 'Wrong number of parameters'

   image = mrdfits(scifile, 0, hdr)
   if (NOT keyword_set(image)) then $
    message, 'Error reading file '+scifile
   invvar = mrdfits(scifile, 1)
   arcstr = strmid(sxpar(hdr, 'ARCFILE'),4,11)
   flatstr = strmid(sxpar(hdr, 'FLATFILE'),4,11)
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile = 'ymodel-test.fits'

   arcfile = (findfile('spArc'+arcstr+'.fits*', count=ct))[0]
   if (ct EQ 0) then $
    message, 'Unable to find spArc file'
   flatfile = (findfile('spFlat'+flatstr+'.fits*', count=ct))[0]
   if (ct EQ 0) then $
    message, 'Unable to find spFlat file'

   basisfile = 'spBasisPSF-' + strmid(sxpar(hdr, 'ARCFILE'),4,11) + '.fits'
   junk = findfile(basisfile, count=ct)
   if (ct EQ 0 OR keyword_set(clobber)) then begin
      splog, 'Generating sdProc file for arc'
      arcname = 'sdR-'+arcstr+'.fit'
      rawdata_dir = getenv('RAWDATA_DIR')
      mjdstr = sxpar(hdr, 'MJD')
      indir = concat_dir(rawdata_dir, mjdstr)
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
   endif

   splog, 'Reading existing traceset '+flatfile
   xset = mrdfits(flatfile, 1)
   traceset2xy, xset, xx, ximg

   splog, 'Running 2D extraction'
   bbspec_extract, image, invvar, flux, fluxivar, basisfile=basisfile, $
    ximg=ximg, ymodel=bb_ymodel, _EXTRA=Extra

   splog, 'Writing file '+outfile
   mwrfits, bb_ymodel, outfile, /create
   mwrfits, flux, outfile
   mwrfits, fluxivar, outfile

   return
end
;------------------------------------------------------------------------------
