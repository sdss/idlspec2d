; coadd the zdcf for a given list of objects and a specific line

pro coadd_zdcf, rmid, line=line, result=result, coadd=coadd, diag=diag

   file='/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
   target=mrdfits(file,1)

   nobj=n_elements(rmid)
   flag=lonarr(nobj)

   ; where the ZDCF results are stored
   topdir='/data2/jli184/ZDCF_ACBFJ/zdcf_data/'

   tt = temporary(result)
   for i=0, nobj-1 do begin

      rmtag='rm'+string(rmid[i],format='(i3.3)')
      file=topdir+rmtag+'_'+line+'_zdcf.dcf'

      if file_test(file) eq 1 then begin
         flag[i]=1
         readcol,file,format='d,d,d, d,d,d', tau, etau_lo, etau_hi, dcf, edcf_lo, edcf_hi, /silent
         if n_elements(tau_arr) eq 0 then tau_arr=tau

         edcf = 0.5*(edcf_lo + edcf_hi)

         tmp={rmid:rmid[i], tau:tau, etau_lo:etau_lo, etau_hi: etau_hi, dcf:dcf, edcf_lo:edcf_lo, edcf_hi:edcf_hi, $
                edcf:edcf }
         if n_elements(result) eq 0 then result=tmp else result=[result, tmp]

      endif

   endfor

   ; now do the coadd with different weights
   nnn=n_elements(result)
   nbin=n_elements(result[0].tau)
 
   ; 1) simple median
   med_dcf = median( result.dcf, dim=2 )
   err_dcf = stddev( result.dcf, dim=2) / sqrt(nnn)

   ; 2) inverse variance weights
   weights = 0.*result.edcf
   ind=where(result.edcf gt 1d-6)
   weights[ind] = 1.D/( (result.edcf)[ind])^2

   mean_w = total( result.dcf*weights, 2, /double)/total(weights, 2, /double)
   err_w = sqrt(total( (result.edcf*weights)^2, 2, /double)) / total(weights, 2, /double)

   ; 3) modified inverse variance weights where those extremely large weights (>5 stddev) are replaced with the median value
   weights2=weights
   for i=0L, nbin - 1 do begin
      sigma=stddev(weights2[i, *])
      indd = where( weights2[i, *] - median(weights2[i, *]) gt 5.*sigma)
      if indd[0] ne -1 then weights2[i, indd] = median(weights2[i, *])
   endfor
   mean_w2 = total( result.dcf*weights2, 2, /double)/total(weights2, 2, /double)
   err_w2 = sqrt(total( (result.edcf*weights2)^2, 2, /double)) / total(weights2, 2, /double)

   ; 4) use the median error in zdcf for each objects as the weight, so that the weight distribution is the same across tau bins
   weights3=weights
   for i=0L, nnn - 1 do begin
      weights3[*,i] = median(weights3[*,i])
   endfor
   mean_w3 = total( result.dcf*weights3, 2, /double)/total(weights3, 2, /double)
   err_w3 = sqrt(total( (result.edcf*weights3)^2, 2, /double)) / total(weights3, 2, /double)
 
   coadd={rmid:result.rmid, line:line, weights:weights, weights2:weights2,weights3:weights3,tau:result[0].tau, etau_lo:result[0].etau_lo, etau_hi:result[0].etau_hi, $
        med_dcf:med_dcf, edcf_med:err_dcf, wmea_dcf:mean_w, edcf_wmea:err_w, wmea_dcf2:mean_w2,edcf_wmea2:err_w2, $
        wmea_dcf3:mean_w3,edcf_wmea3:err_w3}

   ; make diagnostic plots
   target=target[coadd.rmid]
   if keyword_set(diag) then begin
      LoadCT, 0
      img1=result.dcf & img2=result.dcf/result.edcf
      weight_1d=reform(weights3[0,*])
      ind_sort=sort(weight_1d)
      ind_sort=sort(target.zfinal)      

      p = [0.02, 0.1, 0.98, 0.95]
      window, 1, retain=2
      range=[-1.,1.]
      cgimage, img1[*, ind_sort], position=p,/erase, title='ZDCF', minvalue=range[0], maxvalue=range[1]
      xyouts, 0.5, 0.955, /norm, line+'  ZDCF'
      cgColorbar, Position=[p[0], p[1]-0.05, p[2], p[1]-0.01], range=range

      window, 2, retain=2, title='idl 1'
      range=[-10, 30.]
      cgimage, img2[*, ind_sort], position=p,/erase, minvalue=range[0], maxvalue=range[1]
      xyouts, 0.5, 0.955, /norm, line+'  ZDCF/err'
      cgColorbar, Position=[p[0], p[1]-0.05, p[2], p[1]-0.01], range=range

      window, 3, retain=2, title='weights'
      loadct, 4
      range=[0,1]
      cgimage, weights3[*, ind_sort]/max(weights3), position=p,/erase, minvalue=range[0], maxvalue=range[1]
      xyouts, 0.5, 0.955, /norm, line+'  Weight'
      cgColorbar, Position=[p[0], p[1]-0.05, p[2], p[1]-0.01], range=range

      ;message, 'stop'
   endif


end
