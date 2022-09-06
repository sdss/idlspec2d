; Script to re-process epochs with the xyfit custom flux calibration
; Example: rm_combine_script, planfile
; use run2d='v5_7_1' for 2014 data and run2d='v5_10_10' for eBOSS data (2015-)
pro rm_combine_script, planfile, run2d=run2d,skipfluxing=skipfluxing, skipfcorr=skipfcorr, $
     nofcorr=nofcorr,nodist=nodist, method=method, finaldir=finaldir,xyfit=xyfit, $
     loaddesi=loaddesi,legacy=legacy,plates=plates,minsn2=minsn2,bscore=bscore,$
     MWM_fluxer=MWM_fluxer, radec_coadd=radec_coadd, no_reject=no_reject, $
     onestep_coadd=onestep_coadd ;,special=special

RESOLVE_ROUTINE,'sdss_maskbits',/EITHER,/SKIP_EXISTING, /quiet
RESOLVE_ALL, /SKIP_EXISTING, /quiet, /CONTINUE_ON_ERROR;, class='COMMON'
CPU, TPOOL_NTHREADS = 1
;if n_elements(planfile) eq 0 then $

; first determine the proper topdir
fieldstr = (strsplit(repstr(repstr(planfile,'spPlancomb',''),'.par', ''),'-',/extract))[0]
   
if not keyword_set(finaldir) then finaldir = '';'recalib/' ; 'recalib/test20/'
if not keyword_set(xyfit) then xyfit = 1L

if ~keyword_set(run2d) then run2d = getenv('RUN2D')
for i=0L, n_elements(planfile) - 1L do begin
   if keyword_set(legacy) or keyword_set(plates) then begin
       topdir = getenv('BOSS_SPECTRO_REDUX') + '/' + run2d + '/' + fieldstr[i] + '/'
   endif else begin
       topdir = getenv('BOSS_SPECTRO_REDUX') + '/' + run2d + '/' + fieldstr[i] + '/'
   endelse
   rm_spcombine_v5, planfile[i],finaldir=finaldir,xyfit=xyfit, topdir=topdir, $
     skipfluxing=skipfluxing, skipfcorr=skipfcorr, nofcorr=nofcorr, $ 
     nodist=nodist, loaddesi=loaddesi, legacy=legacy,plates=plates, $
     minsn2=minsn2, bscore=bscore, MWM_fluxer=MWM_fluxer, $
     radec_coadd=radec_coadd, no_reject=no_reject, onestep_coadd=onestep_coadd

endfor

end
