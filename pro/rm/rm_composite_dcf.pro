
pro rm_coadd_dcf, dcf_data, tau_dcf=tau_dcf, dcf_coadd=dcf_coadd, edcf_coadd=edcf_coadd, $
     indp_npair_coadd=indp_npair_coadd, line=line

    tau_dcf=dcf_data[0].tau_dcf

    ngrid=n_elements(tau_dcf)

    ; get the DCF for all objects
    tagname=tag_names(dcf_data[0])
    ind=where(strmatch(tagname,  strupcase(line+'_DCF')) )
    dcf=dcf_data.(ind[0])
    ind=where(strmatch(tagname,  strupcase(line+'_EDCF')) )
    edcf=dcf_data.(ind[0])
    ind=where(strmatch(tagname,  strupcase(line+'_indp_npair')) )
    indp_npair=dcf_data.(ind[0])

    ; compute the invvar
    invvar=0.*edcf > 0
    ind=where(edcf gt 0)
    invvar[ind] = 1./ (edcf[ind])^2

    ; force problematic dcf and invvar to be zero
    ind = where(finite(dcf) ne 1)
    dcf[ind]=0. & invvar[ind]=0. 
    ;message, 'stop'

    ; get the total # of independent pairs in each tau bin
    indp_npair_coadd=total(indp_npair,2)

    ; get the invvar weighted average DCF
    sum1 = total(dcf*invvar, 2, /double)
    sum2 = total(invvar, 2, /double)
    ind = where(sum2 gt 0)
    dcf_coadd=dblarr(ngrid)
    dcf_coadd[ind] = sum1[ind]/sum2[ind]

    ; bootstrap to get the uncertainty in the coadded DCF

end

pro rm_composite_dcf, epoch_id=epoch_id

   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   ; keep only RM targets
   target=fibermap[0:848]

   ; set-up the epochs of the LC
   if not keyword_set(epoch_id) then epoch_id=indgen(18)
   ; remove epoch 3 and 7 with the lowest sn
   ;mask_indx=[2,6]
   ;epoch_id[mask_indx]=-1
   ind=where(epoch_id ge 0)
   epoch_id=epoch_id[ind]
   maxep=string(max(epoch_id+1), format='(i2.2)')

   ; get all the LC data
   rm_get_rms_prop,result=LC_data_all, /silent,epoch_id=epoch_id

   nobj=n_elements(target)

   ; read in more target information
   file='/data3/quasar/yshen/work/composite_lag/target_info.fits'
   info=mrdfits(file,1)
    
   tau_dcf=-50.+findgen(41)*5
   ;tau_dcf=-50. + findgen(21)*10.
   ngrid=n_elements(tau_dcf)
   result={rm_ID:-1L, tau_obs:0.D,logL5100:0.D,imag:0.D, Mi_z2:0.D,z:0.D,objc_type:0L,tau_dcf:tau_dcf, $
           hbeta_dcf:dblarr(ngrid),hbeta_edcf:dblarr(ngrid)-1.,hbeta_indp_npair:lonarr(ngrid), $
           mgii_dcf:dblarr(ngrid),mgii_edcf:dblarr(ngrid)-1.,mgii_indp_npair:lonarr(ngrid)}
   result=replicate(result,nobj)
   result.rm_ID=indgen(849)
   result.tau_obs=info.tau_obs
   result.logL5100=info.logL5100
   result.imag=info.imag
   result.Mi_z2=info.Mi_z2
   result.z=info.zfinal
   result.objc_type=info.objc_type

   tagname = tag_names(LC_data_all[0])
   tagname_result=tag_names(result[0])
   mjdlist = LC_data_all[0].mjd
   sncut=0.
   line=['Hbeta','MgII']
   conti=['L1350','L1700','L3000', 'L5100']
   nline=n_elements(line)
   for i=0L, nobj-1 do begin
      
      LC_data = LC_data_all[i]

      ; 
      for jj=0L, nline - 1 do begin
         ind = where(strmatch(tagname,  strupcase(line[jj])+'*') )
         ngood = LC_data.(ind[2])
         LC_arr = LC_data.(ind[3])
         LC_err = LC_data.(ind[4])
         goodflag = LC_data.(ind[5])
         if ngood gt 0 then begin
            ind_good = where(goodflag eq 1)
            LC_arr = 10.D^LC_arr[ind_good]
            LC_err = LC_err[ind_good] ; NB, err in dex
            LC_err = LC_err*alog(10.D)*LC_arr
            mjd_good = mjdlist[ind_good]
            ; remove noisy measurements?
            ind_goodsn = where(LC_arr/LC_err gt sncut)  ; 20% flux measurements, which is too precise
            cadence_id = mjd_good[ind_goodsn]
            lc_line=LC_arr[ind_goodsn] & err2=LC_err[ind_goodsn]
         
            ; get continuum LC, starting from the shortest wavelength
            flag0=0 & iconti=0
            while flag0 eq 0 do begin
               ind = where(strmatch(tagname, strupcase(conti[iconti])+'*') )
               ngood = LC_data.(ind[2])
               LC_arr = LC_data.(ind[3])
               LC_err = LC_data.(ind[4])
               goodflag = LC_data.(ind[5])
               if ngood gt 0 then begin
                  ;ind_good = where(LC_arr gt 0 and LC_err gt 0)
                  ind_good = where(goodflag eq 1)
                  LC_arr = 10.D^LC_arr[ind_good]
                  LC_err = LC_err[ind_good] ; NB, err in dex
                  LC_err = LC_err*alog(10.D)*LC_arr
                  mjd_good = mjdlist[ind_good]
                  ; remove noisy measurements?
                  ind_goodsn = where(LC_arr/LC_err gt sncut) ; 20% flux measurements
                  cadence_id_c = mjd_good[ind_goodsn]
                  lc_conti=LC_arr[ind_goodsn] & err1=LC_err[ind_goodsn]
                  ; we have found the continuum LC
                  flag0=1L
               endif
               iconti=iconti+1
            endwhile

            if ngood gt 0 and flag0 eq 1 then begin
               rm_dcf,cadence_id_c,lc_conti,cadence_id,lc_line,err1=median(err1),err2=median(err2),$
                  DCF=DCF,eDCF=eDCF,tau_grid=tau_dcf,indp_npair=indp_npair
               
               ; assign the DCF
               ind1=where(strmatch(tagname_result,  strupcase(line[jj])+'_DCF') )
               result[i].(ind1[0])=dcf
               ind1=where(strmatch(tagname_result,  strupcase(line[jj])+'_EDCF') )
               result[i].(ind1[0])=edcf
               ind1=where(strmatch(tagname_result,  strupcase(line[jj])+'_INDP_NPAIR') )
               result[i].(ind1[0])=indp_npair
            endif
         endif
      endfor
      splog, 'Finished obj: ', i

   endfor
   
   ; output the results
   outfile = '/data3/quasar/yshen/work/composite_lag/dcf_data_maxep_'+maxep+'.fits'
   mwrfits,result,outfile,/create

end
