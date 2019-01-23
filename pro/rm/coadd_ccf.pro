; coadd the ccf for a given list of objects and a specific line

pro coadd_ccf, rmid, line=line, result=result

   nobj=n_elements(rmid)
   flag=lonarr(nobj)

   topdir='/data3/yshen/ftp/sdssrm/collab/prepspec/ACBFJ/'

   for i=0, nobj-1 do begin

      rmtag='rm'+string(rmid[i],format='(i3.3)')
      file=topdir+rmtag+'/ccf_'+line
      if file_test(file) eq 1 then begin
         flag[i]=1
         readcol,file,format='d,d', tau, ccf, /silent
         if n_elements(tau_arr) eq 0 then tau_arr=tau

         if n_elements(ccf_arr) eq 0 then ccf_arr = ccf $
           else ccf_arr=[ [ccf_arr], [ccf] ] 

      endif

   endfor
   ind_coadd=where(flag eq 1, ngood)

   ccf_coadd=median(ccf_arr, dim=2)
   ccf_coadd_err=stddev(ccf_arr, dim=2)/sqrt(double(ngood))

   result={line:line, tau:tau_arr, ccf_arr:ccf_arr, rmid:rmid[ind_coadd], ccf_coadd:ccf_coadd, $
             ccf_coadd_err:ccf_coadd_err}


end
