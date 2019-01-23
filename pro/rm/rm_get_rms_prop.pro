;+
; NAME:
;   rm_get_rms_prop
; 
; PURPOSE:
;   Get the fractional RMS variations from the qsofit results 
;
;
; CALLING SEQUENCE: 
;   rm_get_rms_prop, [plate,mjd], result=result, [epoch_id=[0,3,4,5,7,8,9,10,11,12]  ]
;
; INPUTS:
; 
; OPTIONAL INPUTS:
;   plate      -  Array of plate numbers of the epochs used to compute RMS
;   mjd        -  Array of MJD of the epochs used to compute RMS
;   epoch_id   -  Array of indices of epochs used to compute RMS (start from 0)
;   rm_ID      -  Indices of objects to compute RMS (same order in fibermap);
;                 default rm_ID=0:848
;
; OUTPUTS:
;   result     -  fractional RMS variation and median errors for selected lines/continuum;
;                 .*_ngood: number of good epochs in computing the RMS
;-----------------
pro rm_get_rms_prop,plate=plate,mjd=mjd,rm_ID,result=result,epoch_id=epoch_id,silent=silent

   if n_elements(plate) ne n_elements(mjd) then begin
      splog, 'Plate and MJD must have the same number of elements. Return.'
      return
   endif

   if n_elements(plate) gt 0 and n_elements(epoch_id) gt 0 then begin
      splog, 'Plate and epoch_id cannot be set at the same time. Return.'
      return
   endif

   topdir='/data3/quasar/yshen/spectro/bossredux/'  + getenv('RUN2D') + '/'

   ; read in the master file
   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   plate_all=fibermap[0].plate
   mjd_all=fibermap[0].mjd
   ind=where(plate_all gt 0)
   plate_all=plate_all[ind]
   mjd_all=mjd_all[ind]
   fiber_all=(fibermap.fiberid)[ind,0:848]

   ; default is to use all epochs
   if n_elements(plate) eq 0 then begin 
      plate=plate_all
      mjd=mjd_all
      if keyword_set(epoch_id) then begin
         plate=plate[epoch_id] & mjd=mjd[epoch_id]
      endif
   endif 
   if n_elements(rm_ID) eq 0 then rm_ID=indgen(849)
   nobj=n_elements(rm_ID)

   ; determin which epochs are included in the RMS computation
   flag=lonarr(n_elements(plate_all))
   for i=0, n_elements(plate)-1 do begin
      ind=where(plate_all eq plate[i] and mjd_all eq mjd[i])
      if ind[0] ne -1 then flag[i] = 1
   endfor
   ind_epoch = where(flag eq 1, nepoch)
   ; override the index of epochs to use
   if keyword_set(epoch_ID) then ind_epoch = epoch_ID
   fiber = fiber_all[*, rm_ID]
   fiber = fiber[ind_epoch,*]

   struct={rm_ID:-1L,plate:plate,mjd:mjd,hbeta_rms:0.D, hbeta_mederr:-1.D, hbeta_ngood:0L, $
             hbeta_LC:dblarr(nepoch),hbeta_err:dblarr(nepoch)-1,hbeta_goodflag:lonarr(nepoch), $
             mgii_rms:0.d, mgii_mederr:-1.D, mgii_ngood:0L, $
             mgii_LC:dblarr(nepoch),mgii_err:dblarr(nepoch)-1,mgii_goodflag:lonarr(nepoch), $
             oiii_rms:0.D, oiii_mederr:-1.D, oiii_ngood:0L, $
             oiii_LC:dblarr(nepoch),oiii_err:dblarr(nepoch)-1,oiii_goodflag:lonarr(nepoch), $
             halpha_rms:0.D, halpha_mederr:-1.D, halpha_ngood:0L, $
             halpha_LC:dblarr(nepoch),halpha_err:dblarr(nepoch)-1,halpha_goodflag:lonarr(nepoch), $
             ciii_rms:0.D, ciii_mederr:-1.D, ciii_ngood:0L, $
             ciii_LC:dblarr(nepoch),ciii_err:dblarr(nepoch)-1,ciii_goodflag:lonarr(nepoch), $
             civ_rms:0.D, civ_mederr:-1.D, civ_ngood:0L, $
             civ_LC:dblarr(nepoch),civ_err:dblarr(nepoch)-1, civ_goodflag:lonarr(nepoch), $
             L1350_rms:0.D, L1350_mederr:-1.D, L1350_ngood:0L, $
             L1350_LC:dblarr(nepoch),L1350_err:dblarr(nepoch)-1, L1350_goodflag:lonarr(nepoch), $ 
             L1700_rms:0.D, L1700_mederr:-1.D, L1700_ngood:0L, $
             L1700_LC:dblarr(nepoch),L1700_err:dblarr(nepoch)-1, L1700_goodflag:lonarr(nepoch), $
             L3000_rms:0.D, L3000_mederr:-1.D, L3000_ngood:0L, $
             L3000_LC:dblarr(nepoch),L3000_err:dblarr(nepoch)-1, L3000_goodflag:lonarr(nepoch), $
             L5100_rms:0.D, L5100_mederr:-1.D, L5100_ngood:0L, $
             L5100_LC:dblarr(nepoch),L5100_err:dblarr(nepoch)-1, L5100_goodflag:lonarr(nepoch) $
                   }
   ;
   result=replicate(struct,nobj)
   result.rm_ID = rm_ID

   ; Hbeta
   arr1 = dblarr(nepoch,nobj) & err1 = dblarr(nepoch,nobj) & redchi2_1=dblarr(nepoch,nobj)
   ; MgII
   arr2 = dblarr(nepoch,nobj) & err2 = dblarr(nepoch,nobj) & redchi2_2=dblarr(nepoch,nobj)
   ; OIII
   arr3 = dblarr(nepoch,nobj) & err3 = dblarr(nepoch,nobj) & redchi2_3=dblarr(nepoch,nobj)
   ; L3000
   arr4 = dblarr(nepoch,nobj) & err4 = dblarr(nepoch,nobj) & redchi2_4=dblarr(nepoch,nobj)
   ; L1350
   arr5 = dblarr(nepoch,nobj) & err5 = dblarr(nepoch,nobj) & redchi2_5=dblarr(nepoch,nobj)
   ; L1700
   arr6 = dblarr(nepoch,nobj) & err6 = dblarr(nepoch,nobj) & redchi2_6=dblarr(nepoch,nobj)
   ; L5100
   arr7 = dblarr(nepoch,nobj) & err7 = dblarr(nepoch,nobj) & redchi2_7=dblarr(nepoch,nobj)
   ; Halpha
   arr8 = dblarr(nepoch,nobj) & err8 = dblarr(nepoch,nobj) & redchi2_8=dblarr(nepoch,nobj)
   ; CIII_all
   arr9 = dblarr(nepoch,nobj) & err9 = dblarr(nepoch,nobj) & redchi2_9=dblarr(nepoch,nobj)
   ; CIV
   arr10 = dblarr(nepoch,nobj) & err10 = dblarr(nepoch,nobj) & redchi2_10=dblarr(nepoch,nobj)

   for i_ep=0, nepoch-1 do begin
   
      platestr=string(plate[i_ep],format='(i4.4)')
      mjdstr=string(mjd[i_ep],format='(i5.5)')
      fitsfile=topdir+platestr+'/qsofit/qso_prop-' + platestr + '-' + mjdstr + '.fits'
      qsofit = mrdfits(fitsfile,1,/silent)
      qsofit = qsofit[rm_ID]
         
      ; do it for Hbeta
      arr1[i_ep,*] = (qsofit.hbeta)[2,*] ; this is logL
      err1[i_ep,*] = (qsofit.hbeta_err)[2,*]  ; this is err in logL
      redchi2_1[i_ep,*] = qsofit.hbeta_redchi2
      ; do it for MgII
      arr2[i_ep,*] = (qsofit.mgii)[2,*] ; this is logL
      err2[i_ep,*] = (qsofit.mgii_err)[2,*]  ; this is err in logL
      redchi2_2[i_ep,*] = qsofit.mgii_redchi2
      ; do it for oiii
      arr3[i_ep,*] = (qsofit.oiii5007)[2,*] ; this is logL
      err3[i_ep,*] = (qsofit.oiii5007_err)[2,*]  ; this is err in logL
      redchi2_3[i_ep,*]=qsofit.oiii5007_redchi2
      ; do it for L3000
      arr4[i_ep,*] = (qsofit.logL3000)[*] ; this is logL
      err4[i_ep,*] = (qsofit.logL3000_err)[*]  ; this is err in logL
      redchi2_4[i_ep,*]=qsofit.conti_redchi2
      ; do it for L1350
      arr5[i_ep,*] = (qsofit.logL1350)[*] ; this is logL
      err5[i_ep,*] = (qsofit.logL1350_err)[*]  ; this is err in logL
      redchi2_5[i_ep,*]=qsofit.conti_redchi2
      ; L1700
      arr6[i_ep,*] = (qsofit.logL1700)[*] ; this is logL
      err6[i_ep,*] = (qsofit.logL1700_err)[*]  ; this is err in logL
      redchi2_6[i_ep,*]=qsofit.conti_redchi2
      ; L5100
      arr7[i_ep,*] = (qsofit.logL5100)[*] ; this is logL
      err7[i_ep,*] = (qsofit.logL5100_err)[*]  ; this is err in logL
      redchi2_7[i_ep,*]=qsofit.conti_redchi2
      ; halpha
      arr8[i_ep,*] = (qsofit.halpha)[2,*] ; this is logL
      err8[i_ep,*] = (qsofit.halpha_err)[2,*]  ; this is err in logL
      redchi2_8[i_ep,*]=qsofit.halpha_redchi2
      ; ciii_all
      arr9[i_ep,*] = (qsofit.ciii_all)[2,*] ; this is logL
      err9[i_ep,*] = (qsofit.ciii_all_err)[2,*]  ; this is err in logL
      redchi2_9[i_ep,*]=qsofit.ciii_all_redchi2
      ; civ
      arr10[i_ep,*] = (qsofit.civ)[2,*] ; this is logL
      err10[i_ep,*] = (qsofit.civ_err)[2,*]  ; this is err in logL
      redchi2_10[i_ep,*]=qsofit.civ_redchi2


   endfor
   if not keyword_set(silent) then splog, 'Finished reading all epochs'

   ;message, 'stop'
   maxredchi2=100. ; maximum value of redchi2 to deem a good fit epoch

   for iobj=0,nobj-1 do begin

      ; for hbeta
      arr = arr1[*, iobj] & err = err1[*,iobj] & redchi2=redchi2_1[*,iobj]
      ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
      if nnn gt 0 then result[iobj].hbeta_goodflag[ind] = 1
      result[iobj].hbeta_ngood=nnn
      result[iobj].hbeta_lc = arr
      result[iobj].hbeta_err = err
      if nnn gt 1 then begin
         arr = arr[ind] & err = err[ind]
         ; covert to linear flux scale
         result[iobj].hbeta_rms=stddev(10.D^arr)/mean(10.D^arr)
         result[iobj].hbeta_mederr=median(err)*alog(10.D)
      endif

      ; for MgII
      arr = arr2[*, iobj] & err = err2[*,iobj] & redchi2=redchi2_2[*,iobj]
      ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
      if nnn gt 0 then result[iobj].mgii_goodflag[ind] = 1
      result[iobj].mgii_ngood=nnn
      result[iobj].mgii_lc = arr
      result[iobj].mgii_err = err
      if nnn gt 1 then begin
         arr = arr[ind] & err = err[ind]
         ; covert to linear flux scale
         result[iobj].mgii_rms=stddev(10.D^arr)/mean(10.D^arr)
         result[iobj].mgii_mederr=median(err)*alog(10.D)
      endif

      ; for OIII
      arr = arr3[*, iobj] & err = err3[*,iobj] & redchi2=redchi2_3[*,iobj]
      ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
      if nnn gt 0 then result[iobj].oiii_goodflag[ind] = 1
      result[iobj].oiii_ngood=nnn
      result[iobj].oiii_lc = arr
      result[iobj].oiii_err = err
      if nnn gt 1 then begin
         arr = arr[ind] & err = err[ind]
         ; covert to linear flux scale
         result[iobj].oiii_rms=stddev(10.D^arr)/mean(10.D^arr)
         result[iobj].oiii_mederr=median(err)*alog(10.D)
      endif

      ; for L3000
      arr = arr4[*, iobj] & err = err4[*,iobj] & redchi2=redchi2_4[*,iobj]
      ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
      if nnn gt 0 then result[iobj].L3000_goodflag[ind] = 1
      result[iobj].L3000_ngood=nnn
      result[iobj].L3000_lc = arr
      result[iobj].L3000_err = err
      if nnn gt 1 then begin
         arr = arr[ind] & err = err[ind]
         ; covert to linear flux scale
         result[iobj].L3000_rms=stddev(10.D^arr)/mean(10.D^arr)
         result[iobj].L3000_mederr=median(err)*alog(10.D)
      endif

      ; L1350
      arr = arr5[*, iobj] & err = err5[*,iobj] & redchi2=redchi2_5[*,iobj]
      ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
      if nnn gt 0 then result[iobj].L1350_goodflag[ind] = 1
      result[iobj].L1350_ngood=nnn
      result[iobj].L1350_lc = arr
      result[iobj].L1350_err = err
      if nnn gt 1 then begin
         arr = arr[ind] & err = err[ind]
         ; covert to linear flux scale
         result[iobj].L1350_rms=stddev(10.D^arr)/mean(10.D^arr)
         result[iobj].L1350_mederr=median(err)*alog(10.D)
      endif

      ; L1700
      arr = arr6[*, iobj] & err = err6[*,iobj] & redchi2=redchi2_6[*,iobj]
      ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
      if nnn gt 0 then result[iobj].L1700_goodflag[ind] = 1
      result[iobj].L1700_ngood=nnn
      result[iobj].L1700_lc = arr
      result[iobj].L1700_err = err
      if nnn gt 1 then begin
         arr = arr[ind] & err = err[ind]
         ; covert to linear flux scale
         result[iobj].L1700_rms=stddev(10.D^arr)/mean(10.D^arr)
         result[iobj].L1700_mederr=median(err)*alog(10.D)
      endif

     ; L5100
     arr = arr7[*, iobj] & err = err7[*,iobj] & redchi2=redchi2_7[*,iobj]
     ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
     if nnn gt 0 then result[iobj].L5100_goodflag[ind] = 1
     result[iobj].L5100_ngood=nnn
     result[iobj].L5100_lc = arr
     result[iobj].L5100_err = err
     if nnn gt 1 then begin
        arr = arr[ind] & err = err[ind]
        ; covert to linear flux scale
        result[iobj].L5100_rms=stddev(10.D^arr)/mean(10.D^arr)
        result[iobj].L5100_mederr=median(err)*alog(10.D)
     endif

     ; halpha
     arr = arr8[*, iobj] & err = err8[*,iobj] & redchi2=redchi2_8[*,iobj]
     ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
     if nnn gt 0 then result[iobj].halpha_goodflag[ind] = 1
     result[iobj].halpha_ngood=nnn
     result[iobj].halpha_lc = arr
     result[iobj].halpha_err = err
     if nnn gt 1 then begin
        arr = arr[ind] & err = err[ind]
        ; covert to linear flux scale
        result[iobj].halpha_rms=stddev(10.D^arr)/mean(10.D^arr)
        result[iobj].halpha_mederr=median(err)*alog(10.D)
     endif

     ; ciii_all
     arr = arr9[*, iobj] & err = err9[*,iobj] & redchi2=redchi2_9[*,iobj]
     ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
     if nnn gt 0 then result[iobj].ciii_goodflag[ind] = 1
     result[iobj].ciii_ngood=nnn
     result[iobj].ciii_lc = arr
     result[iobj].ciii_err = err
     if nnn gt 1 then begin
        arr = arr[ind] & err = err[ind]
        ; covert to linear flux scale
        result[iobj].ciii_rms=stddev(10.D^arr)/mean(10.D^arr)
        result[iobj].ciii_mederr=median(err)*alog(10.D)
     endif
 
     ; civ
     arr = arr10[*, iobj] & err = err10[*,iobj] & redchi2=redchi2_10[*,iobj]
     ind = where(arr gt 0 and err gt 0 and redchi2 lt maxredchi2,nnn)
     if nnn gt 0 then result[iobj].civ_goodflag[ind] = 1
     result[iobj].civ_ngood=nnn
     result[iobj].civ_lc = arr
     result[iobj].civ_err = err
     if nnn gt 1 then begin
        arr = arr[ind] & err = err[ind]
        ; covert to linear flux scale
        result[iobj].civ_rms=stddev(10.D^arr)/mean(10.D^arr)
        result[iobj].civ_mederr=median(err)*alog(10.D)
     endif


   endfor

end

