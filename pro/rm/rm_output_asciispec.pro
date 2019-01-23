;+
; PURPOSE:
;   Output ascii spectra of the spPlate files
;
; Keywords:
;
;   calibdir=''            ; BOSS pipeline reduction
;   calibdir='recalib/'    ; xyfit reduction
;   calibdir='wh_skysub/'  ; xyfit reduction + WH05 sky subtraction
;   epoch=[33,38]          ; only output epochs33 to epoch38
;---------------------------------

pro rm_output_asciispec, calibdir=calibdir,outdir=outdir,clobber=clobber $
     ,run2d=run2d,run1d=run1d,epoch=epoch, pt_corr=pt_corr, rmid=rmid


   ; set run version
   if not keyword_set(run2d) then run2d=getenv('RUN2D')
   if not keyword_set(run1d) then run1d=getenv('RUN1D')

   ; path for spPlate* files
   if not keyword_Set(calibdir) then calibdir='wh_skysub/'

   ; path for output ascii files
   if not keyword_Set(outdir) then outdir='/data3/yshen/ftp/sdssrm/collab/bossredux/v5_7_1/ascii_spec/wh_skysub/'
        ;'/data3/quasar/yshen/ftp/bossredux/v5_7_1/ascii_spec/wh_skysub/'

   target_file=getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   info=mrdfits(target_file,2,/silent)
   mean_mjd=info.mean_mjd
   mjdall = fibermap[0].mjd

   if ~keyword_set(epoch) then begin
     ; get all the available epochs
     ind = where(fibermap[0].plate gt 0, nnn)
     epoch=[1,nnn]
   endif 
   nnn=epoch[1]-epoch[0]+1

   fmt='(f10.3, e12.4, e12.4,i12,i12)'
   if ~keyword_set(rmid) then rmid=indgen(1000L)
   nobj = n_elements(rmid)
   for i=0L, nobj - 1 do begin
      specdir = outdir + 'RMID_' + string(rmid[i],format='(i3.3)') + '/'
      if file_test(specdir) eq 0 then spawn, 'mkdir ' + specdir
      if keyword_set(pt_corr) then begin ; apply the pt corrections from PrepSpec
         pt=rm_get_pt(rmid[i],err=pt_err,mjd_all=mjdall[epoch[0]-1:epoch[1]-1],topdir=prepspec_dir)
      endif   
      for jep=epoch[0]-1, epoch[1]-1 do begin

         plate=(fibermap[rmid[i]].plate)[jep]
         fiber=(fibermap[rmid[i]].fiberid)[jep]
         mjd=(fibermap[rmid[i]].mjd)[jep]
         platestr=string(plate,format='(i4.4)')
         fiberstr=string(fiber,format='(i4.4)')
         mjdstr=string(mjd,format='(i5.5)')
         specfile=specdir+'RMID_' + string(rmid[i],format='(i3.3)') + '_EP' + string(jep+1,format='(i2.2)') + $
            '_' + platestr + '-' + mjdstr + '-' + fiberstr + '.dat'
         gzipfile=specfile+'.gz'
         if file_test(gzipfile) eq 0 or keyword_set(clobber) then begin
            rm_readspec,plate,fiber,mjd=mjd,wave=wave,flux=flux,flerr=flerr,andmask=andmask,ormask=ormask,$
               calibdir=calibdir,run2d=run2d
            if keyword_set(pt_corr) then begin ; apply the pt corrections
               if pt[jep] gt 1d-2 and pt[jep] lt 100 and pt_err[jep] gt 0 then begin
                  flux = flux/pt[jep]
                  ind_good = where(flerr gt 1d-5)
                  if ind_good[0] ne -1 then $
                    flerr[ind_good] = sqrt( (flerr[ind_good]/pt[jep])^2 + (flux[ind_good]*pt_err[jep]/pt[jep]^2)^2 )
               endif
            endif

            npix=n_elements(wave)
            openw, lun, specfile,/get_lun
            printf,lun,'# Mean_MJD='+string(mean_mjd[jep],format='(f9.3)')
            printf,lun,'# RMID='+string(rmid[i],format='(i3.3)') + ' SOURCE_TYPE=' + fibermap[rmid[i]].sourcetype
            printf,lun,'# ZFINAL='+string(fibermap[rmid[i]].zfinal,format='(f6.4)')
            printf,lun,'# EPOCH='+string(jep+1,format='(i2)')
            printf,lun,'# PLATE='+platestr + ' MJD='+mjdstr + ' FIBERID='+fiberstr
            printf,lun,'# Wave[vaccum]  Flux  Flux_Err  ANDMASK  ORMASK'
            printf,lun,'# [Ang]         [10^(-17) erg/cm2/s/Ang]'
            for jj=0L, npix-1 do printf,lun,wave[jj],flux[jj],flerr[jj],andmask[jj],ormask[jj],format=fmt
            close, lun
            free_lun, lun
            spawn, 'gzip -f ' + specfile
         endif
      endfor
      
      splog, 'finished obj: ', i+1, '/', nobj

   endfor

end
