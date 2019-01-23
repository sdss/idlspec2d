;+
; NAME
;
; PURPOSE:
;   output the xpos and ypos focal positions for each target
;

pro rm_output_focal_pos,plate,mjd,outdir=outdir,run2d=run2d,run1d=run1d,calibdir=calibdir

   if ~keyword_set(outdir) then outdir='/data3/quasar/yshen/ftp/bossredux/v5_7_1/ascii_spec/'

   ; set spPlate path
   if ~keyword_set(run2d) then run2d=getenv('RUN2D')
   if ~keyword_set(run1d) then run1d=getenv('RUN1D')
   if ~keyword_Set(calibdir) then calibdir='wh_skysub/'
 
   ; get target list
   target_file=getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   info=mrdfits(target_file,2,/silent)
   mean_mjd=info.mean_mjd

   fiber_all=fibermap.fiberid
   plate_all=(fibermap.plate)[*,0]
   mjd_all=(fibermap.mjd)[*,0]

   ind=where(plate_all eq plate and mjd_all eq mjd)
   fiber=fiber_all[ind,*]

   rm_readspec,plate,fiber,mjd=mjd,plugmap=plugmap,calibdir=calibdir

   ; extended source flag
   eflag = lonarr(1000L)
   ind=where(fibermap.objc_type eq 3)
   eflag[ind]=1L

   fmt='(i3.3, " ", i4.4, " ", f10.6, " ", f10.6, " ", f6.4, " ", A8, " ", i0, " ", f9.4, " ", f9.4)'

   outfile=outdir+string(plate,format='(i4.4)')+'-'+string(mjd,format='(i5.5)')+'_xypos.dat'
   openw, lun, outfile, /get_lun
   printf,lun, '# RMID  FIBERID RA  DEC  ZFINAL   SOURCE_TYPE  EXTENDED? XFOCAL  YFOCAL'
   for i=0L, 999L do begin
      printf, lun, format=fmt, i,fiber[i], fibermap[i].ra,fibermap[i].dec,fibermap[i].zfinal, $
          fibermap[i].sourcetype, eflag[i], plugmap[i].xfocal,plugmap[i].yfocal
   endfor
   close, lun
   free_lun, lun

end
