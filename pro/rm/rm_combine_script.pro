;+
; NAME:
;   rm_combine_script
;
; PURPOSE:
;   Script to process epochs with the xyfit custom flux calibration
;
; CALLING SEQUENCE:
;
; INPUTS:
;   planfile   - Name(s) of output plan file
;
; OPTIONAL INPUTS:
;   run2d      - Name of the run2d
;   finaldir   - Additional subdirectory for output
;   xyfit      - Compute 2d flux corrections in the xy focal plane
;   bscore     - Fraction of best exposure score to use as a threshold for discarding exposures
;   minsn2     - Minimum S/N^2 to include science frame in coadd; default
;                to 0 to only include those with S/N > 0.
;                Note that all exposures with a score less than 0.2 times
;                the score of the best exposure are discarded; for those
;                purposes, the score used is the worst of all 4 cameras.
;
;
; Optional Keywords:
;   MWM_fluxer    - Utilize MWM optional settings (ie gaia reddening and different S/N cuts)
;   nofcorr       - Skip the step to generate and use the spFluxcorr* files
;   nodist        - Skip the step to generate and use the spFluxdistort* files
;   radec_coadd   - Coadd using ra-dec matching rather then catalogID matching
;   no_reject     - Turns off rejection in the coadding
;   onestep_coadd - Legacy algorithm for coadd. Coadding blue+red and all exposures
;                    at the the same time.
;   epoch         - Epoch Coadd flag for input and outputs
;   legacy        - Flag for Pre-SDSSV 2 Spectrograph data at APO
;   plates        - Flat for SDSSV 1 Spectrograph plate data at APO
;   loaddesi      - Load the DESI (JG) models for fluxing
;   skipfluxing   - Skip the step to generate spFluxcalib* files
;   skipfcorr     - Skip creation of flux-correction vectors and use prexisting spFluxcorr* files
;
; OUTPUT:
;
; COMMENTS:
; EXAMPLES:
;
; BUGS:
;   This routine spawns the Unix command 'mkdir'.
;
; PROCEDURES CALLED:
;   get_field_dir
;   djs_filepath
;   rm_spcombine_v5
;
;
;------------------------------------------------------------------------------

pro rm_combine_script, planfile, run2d=run2d,skipfluxing=skipfluxing, skipfcorr=skipfcorr, $
     nofcorr=nofcorr,nodist=nodist, finaldir=finaldir,xyfit=xyfit, $
     loaddesi=loaddesi,legacy=legacy,plates=plates,minsn2=minsn2,bscore=bscore,$
     MWM_fluxer=MWM_fluxer, radec_coadd=radec_coadd, no_reject=no_reject, $
     onestep_coadd=onestep_coadd, epoch=epoch

RESOLVE_ROUTINE,'sdss_maskbits',/EITHER,/SKIP_EXISTING, /quiet
RESOLVE_ALL, /SKIP_EXISTING, /quiet, /CONTINUE_ON_ERROR;, class='COMMON'
CPU, TPOOL_NTHREADS = 1

; first determine the proper topdir
if not keyword_set(epoch) then begin
    fieldstr = (strsplit(repstr(repstr(planfile,'spPlancomb',''),'.par', ''),'-',/extract))[0]
endif else begin
    fieldstr = (strsplit(repstr(repstr(planfile,'spPlancombepoch',''),'.par', ''),'-',/extract))[0]
endelse
if not keyword_set(finaldir) then finaldir = '';'recalib/' ; 'recalib/test20/'
if not keyword_set(xyfit) then xyfit = 1L

if ~keyword_set(run2d) then run2d = getenv('RUN2D')
for i=0L, n_elements(planfile) - 1L do begin
   topdir = get_field_dir(getenv('BOSS_SPECTRO_REDUX'), run2d, fieldstr[i])
   if keyword_set(epoch) then begin
        topdir = djs_filepath('epoch', root_dir = topdir)
   endif
   rm_spcombine_v5, planfile[i],finaldir=finaldir,xyfit=xyfit, topdir=topdir, $
     skipfluxing=skipfluxing, skipfcorr=skipfcorr, nofcorr=nofcorr, $ 
     nodist=nodist, loaddesi=loaddesi, legacy=legacy,plates=plates, $
     minsn2=minsn2, bscore=bscore, MWM_fluxer=MWM_fluxer, epoch=epoch, $
     radec_coadd=radec_coadd, no_reject=no_reject, onestep_coadd=onestep_coadd

endfor

end
