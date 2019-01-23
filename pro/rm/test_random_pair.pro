;+
; NAME:
;  test_random_pair
;
; PURPOSE:
;  perform batch tests for PrepSpec light curves
;  use shuffled random continuum/line LC pairs
;  set /sym_range to search the lag in symmetric range of lags, this however will 
;  increase the uncertainty in lag centroid with MC method, and rendering the lag
;  undetected. 
;
pro test_random_pair,rmid,topdir=topdir, use_epochid=use_epochid, $
     dcf=dcf, ep_rej=ep_rej, ccf_tau=ccf_tau, sym_range=sym_range,outfile=outfile, $
     random_lc=random_lc

   if n_elements(rmid) eq 0 then begin ; default is the first 100 objects
       file = getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
       target=mrdfits(file,1)
       target=target[0:848]
       ind=sort(target.zfinal)
       rmid=ind[0:99]
   endif

   if ~keyword_Set(topdir) then topdir='/data3/quasar/yshen/work/lags/prepspec/ACBFJ/'

   dt=2.
   if n_elements(sym_range) eq 0 then sym_range=1 ; default symmetric search range
   if keyword_set(sym_range) then ccf_tau = -70. + (findgen(132./dt + 5L)*dt)
   if ~keyword_set(ccf_tau) then ccf_tau = -20. + (findgen(100./dt + 6L)*dt)

   ; open an outfile to output the measurements
   if ~keyword_Set(outfile) then outfile1=topdir + 'ccf_random_pairing' else $
     outfile1=topdir + outfile
   openw, outlun, outfile1, /get_lun
   printf, outlun, '# Tau search range: ',min(ccf_tau),max(ccf_tau)
   if n_elements(ep_rej) gt 0 then printf, outlun, '# Rejected epochs: ', ep_rej
   printf, outlun, '# RMID line  linechi2  lag_mea  peak_sig'
   fmt='(i3.3, " ", a4, " ", e8.2, " ", 3(f5.1, " "), e9.3)'

   nep_tot=32 ; this is the current # of spec epochs
   if n_elements(use_epochid) eq 0 then begin ; assign usable epochs
      if n_elements(ep_rej) gt 0 then begin
         epochid = indgen(nep_tot)
         flag = lonarr(nep_tot)
         flag[ep_rej-1]=1L
         ind_keep=where(flag eq 0, complement=ind_bad)
         use_epochid=ind_keep
      endif
   endif

   ; read in more target information
   file='/data3/quasar/yshen/work/composite_lag/target_info.fits'
   info=mrdfits(file,1)

   if n_elements(rmid) eq 0 then rmid=lindgen(849)
   nobj=n_elements(rmid)

   suffix=['hb_t','mg2_t']
   line=['hb','mg2']
   nline=n_elements(suffix)
   for i=0L, nobj-1 do begin

     rmtag='rm'+string(rmid[i],format='(i3.3)')

     ; first get 5100 LC
     lcfile=topdir+rmtag+'/'+rmtag+'_c*.dat'
     lcfile_found=file_search(lcfile,count=nfile)

     if nfile gt 0 then begin
        ; choose the best continuum LC
        conti_tag=strmid(lcfile_found,7,4,/reverse)
        ;print, conti_tag
        ind_use=where(conti_tag eq '5100')
       
        if ind_use[0] eq -1 then ind_use=where(conti_tag eq $
           string(median(conti_tag),format='(i0)' ))

        lcfile_use=lcfile_found[ind_use]
        readcol,lcfile_use,format='d,d,d',t_conti,f_conti,e_conti,/silent
        ; remove bad epochs
        if keyword_Set(use_epochid) then begin
           t_conti=t_conti[use_epochid]
           f_conti=f_conti[use_epochid]
           e_conti=e_conti[use_epochid]
        endif
        ind_good = where(e_conti gt 0)
        if ind_good[0] ne -1 then begin
           t_conti=t_conti[ind_good]
           f_conti=f_conti[ind_good]
           e_conti=e_conti[ind_good]
        endif

        ; now do each BLR LC
        ind_other = where(rmid ne rmid[i], nother)
        for jj=0L, nother - 1 do begin

          rmtag_other = 'rm'+string(rmid[ind_other[jj]],format='(i3.3)')
          for j=0,nline-1 do begin

          ; find the line LC for the other 99 objects

            lcfile=topdir+rmtag_other+'/'+rmtag_other+'_'+suffix[j]+'.dat'
            if file_test(lcfile) eq 1 then begin
               readcol,lcfile,format='d,d,d',t_line,f_line,e_line,/silent
               ; remove bad epochs
               if keyword_Set(use_epochid) then begin
                 t_line=t_line[use_epochid]
                 f_line=f_line[use_epochid]
                 e_line=e_line[use_epochid]
               endif
               ind_good = where(e_line gt 0)
               if ind_good[0] ne -1 then begin
                 t_line=t_line[ind_good]
                 f_line=f_line[ind_good]
                 e_line=e_line[ind_good]
               endif

               if keyword_set(random_lc) then begin ; shuffle the line LC epochs to test if a correlation still exists
                 nep_line=n_elements(t_line)
                 Nran = randomu(seed, nep_line)
                 ind_ran=sort(Nran)
                 f_line=f_line[ind_ran] & e_line=e_line[ind_ran] 
               endif

               tag_plot='c'+conti_tag[ind_use]+'-'+suffix[j]
            
               ; compute the line chi^2
               ; I believe the mean of the line LC has already been subtracted in prepspec
               line_chi2=total( f_line^2/e_line^2, /double )

               rm_lag_diag,t_conti,f_conti,e_conti,t_line,f_line,e_line,dt=dt,tag=tag_plot, $
                /noplot, tau_exp=info[rmid[i]].tau_obs,dcf=dcf, ccf_tau=ccf_tau, $
                lag_mea=lag_mea,peak_sig=peak_sig,ntrial=100,/mirror,/dosub

               printf, outlun, rmid[ind_other[jj]], line[j], line_chi2, $
                 lag_mea, peak_sig, format=fmt

            endif
          endfor
        endfor
      endif

   endfor

   close, outlun
   free_lun, outlun

end
