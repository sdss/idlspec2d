; a tiny fraction of the Hbeta and Halpha fits were bad due to sensitivity to initial guesses; 
; ~1% of objects were affected (both measurements and error estimation)
; find these fits and re-fit with the updated rm_qsofit routine

pro rm_fix_badfits, istart=istart

   ; read in the master file
   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)

   badfits_file = getenv('IDLRM_DIR') + '/misc/badfits.fits'
   tmp=file_test(badfits_file)
   if tmp eq 0 then begin ; find all the bad fits and compile a list
      epoch_id=indgen(18) 
      rm_get_rms_prop,result=result,epoch_id=epoch_id
      ; find all the objects with Hbeta or Halpha fits
      ; ind = where(result.hbeta_ngood gt 0 or result.halpha_ngood gt 0,nobj) 
      ; all 849 objects
      ind = indgen(849) 
      nobj = n_elements(ind)
      nepoch=n_elements(epoch_id)
      rm_ID = result[ind].rm_id
      plate = (fibermap[rm_id].plate)[epoch_id, *]
      fiber = (fibermap[rm_id].fiberid)[epoch_id,*]
      mjd = (fibermap[rm_id].mjd)[epoch_id,*]

      obj_refit = {rm_ID:-1L, epoch_badfit:lonarr(nepoch), refit:0L}
      obj_refit = replicate(obj_refit, nobj)
      obj_refit.rm_ID = rm_ID

      ntot=0L
      for i=0L, nobj-1 do begin
          epoch_badfit = lonarr(nepoch)
          for j=0L, nepoch - 1 do begin
             topdir = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
           + '/' + string(plate[j,i],format='(i4.4)') + '/qsofit/fits/' 
             fitsfile=topdir+string(plate[j,i],format='(i4.4)') + '-' + $
                string(mjd[j,i],format='(i5.5)') + '-' + $
                string(fiber[j,i],format='(i4.4)') + '.fits'
             if file_test(fitsfile) eq 1 then begin
                para = mrdfits(fitsfile, 1, /silent)
                ;indd = where(strmatch(para.linename,'Halpha*') or strmatch(para.linename,'Hbeta*'))
                redchi2=para.line_redchi2
                if max(redchi2) gt 5. then begin 
                   obj_refit[i].refit = 1L
                   epoch_badfit[j] = 1L
                   ntot = ntot + 1L
                   splog, 'Object refit:', obj_refit[i].rm_id,plate[j,i],mjd[j,i],fiber[j,i], ntot
                endif
             endif
          endfor
          obj_refit[i].epoch_badfit = epoch_badfit
      endfor
   
      mwrfits, obj_refit, badfits_file, /create
   endif

   ; start re-fit these bad fits
   obj_refit=mrdfits(badfits_file,1)
   ; keep only those that need refit
   ind = where(obj_refit.refit eq 1)
   obj_refit = obj_refit[ind]

   nnn=n_elements(obj_refit)

   ; Setup error handeling
   errfile = getenv('IDLRM_DIR') + '/misc/' + 'refit.err.out.mjd' + string(current_mjd(), format='(f9.3)')
   openw, lun, errfile, /get_lun

   if n_elements(istart) eq 0 then istart=0L
   for i=istart, nnn - 1L do begin

      rm_id = obj_refit[i].rm_id
      plate = fibermap[rm_id].plate
      fiber = fibermap[rm_id].fiberid
      mjd = fibermap[rm_id].mjd
      ind = where(plate gt 0, nfit)
      plate=plate[ind] & fiber=fiber[ind] & mjd=mjd[ind]

      ; the idea is to re-fit all epochs of this particular object
      for j=0L, nfit - 1 do begin

         ; Estabilish error handlers
         Catch, error_status
         if error_status ne 0 then begin
            ;print, 'Error index: ', error_status
            printf, lun, 'rm_ID plate mjd fiber'
            printf, lun, rm_id, plate[j], mjd[j], fiber[j]
            printf, lun, 'Error message: ', !error_state.msg
            printf, lun, ' '
            catch, /cancel
            continue
         endif

         topdir = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
           + '/' + string(plate[j],format='(i4.4)') + '/'
         outdir = topdir + 'qsofit/'
         errdir = outdir + 'err/'
         rm_fit1fiber,plate[j],fiber[j],mjd[j], $
            ra=fibermap[rm_id].ra,dec=fibermap[rm_id].dec, $
            zfit=fibermap[rm_id].zfinal,rm_indd=rm_id, $
            outdir=outdir, /fits,/psplot, errdir=errdir,ntrial=50L

      endfor
   endfor

   ; close error log
   close, lun
   free_lun, lun

end
