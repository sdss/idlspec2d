; Plot the lag for a list of objects
; ep_rej is the epoch number to be rejected, i.e., ep_rej=7 is the lowest s/n epoch

pro rm_plot_lag, rmid, topdir=topdir,dcf=dcf

if ~keyword_Set(topdir) then topdir='/data3/quasar/yshen/work/lags/prepspec/ACBFJ/'

; rmid=[]
; call rm_plot_lag_one one-by-one

file1=topdir+'/output/lc_all' ;'/data3/quasar/yshen/work/lags/prepspec/output/paper/lc_all'
openw,lun1,file1,/get_lun
printf,lun1,'#RMID line redshift MJD-50000 f_conti e_conti f_line e_line bad?'
file2=topdir+'/output/tau_all' ; '/data3/quasar/yshen/work/lags/prepspec/output/paper/tau_all'
openw,lun2,file2,/get_lun
printf, lun2, '#RMID line redshift tau_rest peak_sig sig_rms sig_mean fwhm_rms fwhm_mean logVP'

rm_plot_lag_one2,775,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,tau_input=[19.18,4.34,-12.79]*(1.+ 0.1725)
rm_plot_lag_one2,769,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7],tau_input=[21.54,5.80,-7.67]*(1. + 0.1871)
rm_plot_lag_one2,840,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,tau_input=[10.93, 20.93,-6.57]*(1. + 0.2439)

;rm_plot_lag_one2,768,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,vrange=[-15000.,10000.]  ; Hb blends with HeII, but it is a good lag detection
rm_plot_lag_one2,272,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf, ep_rej=[7],tau_input=[21.91,7.90,-10.36]*(1. + 0.2628)
rm_plot_lag_one2,320,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf, tau_input=[29.58,2.50,-15.72]*(1. + 0.2647)
;rm_plot_lag_one2,377,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7] ; hb marginally detected
;rm_plot_lag_one2,160,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7],yrange_ccf=[0.5,1], tau_input=[10.08,11.33,-16.18]*(1. + 0.3593)
rm_plot_lag_one2,789,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7],tau_input=[17.15,2.71,-2.73]*(1. + 0.4253)

rm_plot_lag_one2,191,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,tau_input=[23.25,2.72,-11.16]*(1. + 0.4418)
;rm_plot_lag_one2,101,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf ; good
rm_plot_lag_one2,101,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,tau_input=[36.68,10.40,-4.76]*(1. + 0.4581)
;rm_plot_lag_one2,229,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf
rm_plot_lag_one2,229,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,tau_input=[32.33,12.93,-5.31]*(1. + 0.4696)
;rm_plot_lag_one2,371,line='ha',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf
rm_plot_lag_one2,645,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7],yrange_ccf=[0.6,1], tau_input=[14.22,6.55,-8.07]*(1. + 0.4738)

rm_plot_lag_one2,767,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7],tau_input=[25.09,2.00,-2.59]*(1. + 0.5266)

rm_plot_lag_one2,694,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf, tau_input=[14.14,12.93,-9.46]*(1. + 0.5324)
;rm_plot_lag_one2,301,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7]
;rm_plot_lag_one2,519,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[3,7] ; good, HeII/Hbeta blends in the rms model
rm_plot_lag_one2,267,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[3,7],tau_input=[18.57,7.15,-3.80]*(1. + 0.5872)
;rm_plot_lag_one2,457,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[3,7] ; OK
rm_plot_lag_one2,457,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[3,7],tau_input=[29.07,3.64,-8.78]*(1. + 0.6037)

;rm_plot_lag_one2,634,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf ; var hb barely detected
;rm_plot_lag_one2,551,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7] ; var mgii barely detected
;rm_plot_lag_one2,632,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7] ; rising red wing
;rm_plot_lag_one2,510,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7]
;rm_plot_lag_one2,589,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[3,7]
rm_plot_lag_one2,589,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[3,7],tau_input=[34.05,6.72,-12.03]*(1. + 0.7510)
;rm_plot_lag_one2,762,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf,ep_rej=[7]
;;rm_plot_lag_one,797,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[7,29,30,31,32]
;rm_plot_lag_one2,272,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf
;rm_plot_lag_one2,320,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf
;rm_plot_lag_one2,252,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[7],dcf=dcf
;rm_plot_lag_one2,377,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[7],dcf=dcf
;rm_plot_lag_one2,160,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[7],dcf=dcf
;rm_plot_lag_one2,101,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf
;rm_plot_lag_one2,518,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf
;rm_plot_lag_one2,229,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf
;rm_plot_lag_one2,371,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf
;rm_plot_lag_one2,767,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[7],dcf=dcf
;rm_plot_lag_one2,694,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,dcf=dcf
;rm_plot_lag_one2,519,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[3,7],dcf=dcf
;rm_plot_lag_one2,457,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[7],dcf=dcf
;rm_plot_lag_one2,551,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[7],dcf=dcf
;rm_plot_lag_one2,762,line='hb',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[7],dcf=dcf
;rm_plot_lag_one2,140,line='mg2',lun_tau=lun2,lun_lc=lun1,topdir=topdir,ep_rej=[7],dcf=dcf

close,lun1,lun2
free_lun, lun1,lun2

end

pro rm_plot_lag_one,rmid,topdir=topdir,ep_rej=ep_rej,line=line,conti=conti, $
   dt=dt,ccf_tau=ccf_tau,mirror=mirror,xrange_ccf=xrange_ccf,yrange_ccf=yrange_ccf, $
   dcf=dcf, lun_tau=lun_tau, lun_lc=lun_lc

  if ~keyword_set(topdir) then topdir='/data3/quasar/yshen/work/lags/prepspec/'
  if ~keyword_set(conti) then conti='c5100'
  if ~keyword_set(line) then line='hb'

  common rm_info_block, rm_info
  if n_elements(rm_info) eq 0 then begin
    target_file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
    rm_info=mrdfits(target_file,1,/silent)
  endif

  rmtag='rm'+string(rmid,format='(i3.3)')
  ; get continuum LC
  lcfile=topdir+rmtag+'/'+rmtag+'_'+conti+'.dat'
  readcol,lcfile,format='d,d,d',t_conti,f_conti,e_conti,/silent
  nep_c=n_elements(t_conti)
  ; get line LC
  lcfile=topdir+rmtag+'/'+rmtag+'_'+line+'_t.dat'
  readcol,lcfile,format='d,d,d',t_line,f_line,e_line,/silent
  nep_l=n_elements(t_line)
  ; get rms sigma
  sigfile=topdir+rmtag+'/'+rmtag+'_vvblr.dat'
  readcol,sigfile,format='d,x,x,x,d,d',wave_line,fwhm_rms_all, sig_rms_all,skip=1L
  if line eq 'hb' then wave_s=4861.327
  if line eq 'mg2' then wave_s=2798.000
  ind_line=where(abs(wave_line - wave_s) le 1.)
  sig_rms=sig_rms_all[ind_line]

  ; remove bad epoch if required
  flag=lonarr(nep_c)
  t_conti_all=t_conti & f_conti_all=f_conti & e_conti_all=e_conti
  t_line_all=t_line & f_line_all=f_line & e_line_all=e_line
  nep_c_all=n_elements(t_conti_all)
  nep_l_all=n_elements(t_line_all)
  if n_elements(ep_rej) ne 0 then begin
     flag[ep_rej-1]=1L
     ind=where(flag eq 0, complement=indd)
     if indd[0] ne -1 then begin ; get the bad epochs
       t_conti_bad=t_conti[indd] & f_conti_bad=f_conti[indd] & e_conti_bad=e_conti[indd]
       t_line_bad=t_line[indd] & f_line_bad=f_line[indd] & e_line_bad=e_line[indd]
     endif
     t_conti=t_conti[ind] & f_conti=f_conti[ind] & e_conti=e_conti[ind]
     t_line=t_line[ind] & f_line=f_line[ind] & e_line=e_line[ind]
     nep_c=n_elements(t_conti)
     nep_l=n_elements(t_line)
  endif

  ; output the LCs
  outdir=topdir+'output/'
  if ~keyword_set(lun_lc) then begin
    fmt='(f11.3, " ", e14.7, " ", e14.7, " ", i0)'
    lcfile=outdir+rmtag+'_'+conti+'.dat'
    openw,lun,lcfile,/get_lun
    for i=0L,nep_c_all-1 do printf,lun,t_conti_all[i],f_conti_all[i],e_conti_all[i], $
      flag[i],format=fmt
    close,lun
    free_lun,lun
    lcfile=outdir+rmtag+'_'+line+'.dat'
    openw,lun,lcfile,/get_lun
    for i=0L,nep_l_all-1 do printf,lun,t_line_all[i],f_line_all[i],e_line_all[i], $
      flag[i],format=fmt
    close,lun
    free_lun,lun
  endif else begin
    fmt='(i3.3, " ", f11.3, 2(" ", e14.7, " ", e14.7), " ", i0)'
    for i=0L,nep_c_all-1 do printf,lun_lc,rmid,t_conti_all[i],f_conti_all[i],e_conti_all[i],$
      f_line_all[i],e_line_all[i],flag[i],format=fmt
  endelse

  ; now compute the CCF
  if ~keyword_set(dt) then dt=2.
  if ~keyword_set(ccf_tau) then ccf_tau = -10. + (findgen(70./dt + 5L)*dt)
  ; run xcorr_interp once to get the cent_tau and peak sig
  ccf0 = xcorr_interp(t_line,f_line,t_conti,f_conti,tau=ccf_tau,ndata=ndata0 $
            , /nogap,cent_tau=cent_tau0, peak_tau=peak_tau0 $
            , mirror=mirror,peak_sig=peak_sig)
  ; run again to get error estimates
  ccf = xcorr_interp(t_line,f_line,t_conti,f_conti,tau=ccf_tau,ndata=ndata $
            , /nogap,err1=e_line,err2=e_conti, /bootstrap,cent_tau=cent_tau, peak_tau=peak_tau $
            , mirror=mirror,ntrial=500L)

  ; now generate a plot
  ; the map between line names and what to plot
  line1=['blr','ha','hb','mg2','he2']
  figfile=outdir+rmtag+'_'+line+'_lag.eps'
  charsize=1.5 & ticklen=0.05
  begplot,name=figfile,/cmyk,/encap,xsize=12.5, ysize=4
  line2=textoidl(['blr','H\alpha', 'H\beta', 'MgII', 'HeII'])
  indx_line=where(line1 eq line)
  pos11=[0.08,0.55,0.49,0.91] & pos2=[0.57, 0.16, 0.98, 0.91]
  pos12=[0.08,0.16,0.49,0.52]
  xrange=[min(t_conti),max(t_conti)]
  yrange=[min(f_conti),max(f_conti)]
  tag='J'+radec2sdssname(rm_info[rmid].ra,rm_info[rmid].dec)
  plot,[0],[0],xrange=xrange,yrange=yrange, $
    title=tag, pos=pos11, charsize=charsize,xticklen=ticklen, $
    xtickname=replicate(' ', 6L)
  oploterror,t_conti,f_conti,e_conti,psym=symcat(9)
  if n_elements(t_conti_bad) ne 0 then oploterror,[t_conti_bad],[f_conti_bad],[e_conti_bad],psym=symcat(9), $
    color=cgcolor('red'),errcolor=cgcolor('red')
  xyouts, pos11[2] + 0.01, pos11[3] - 0.5*(pos11[3]-pos11[1]), conti, /norm,charsize=charsize, $
    orien=90, align=0.5
  yrange=[min(f_line),max(f_line)]
  plot,[0],[0],xrange=xrange,yrange=yrange,xtitle='MJD-50000 [days]', $
    pos=pos12,charsize=charsize,xticklen=ticklen,xtickformat='(i0)',/noerase
  oploterror,t_line,f_line,e_line,psym=symcat(9),color=cgcolor('blue')
  if n_elements(t_line_bad) ne 0 then oploterror,[t_line_bad],[f_line_bad],[e_line_bad],psym=symcat(9), $
    color=cgcolor('red'),errcolor=cgcolor('red')
  xyouts, pos12[2] + 0.01, pos12[3] - 0.5*(pos12[3]-pos12[1]), line2[indx_line], /norm,charsize=charsize,color=cgcolor('blue'), orien=90, align=0.5
  xyouts, 0.02, 0.5, 'Relative Flux [arbitrary units]',/norm,orientation=90,align=0.5,charsize=charsize
  ; now plot the CCF
  tag='RMID='+string(rmid,format='(i3.3)') + ', z='+string(rm_info[rmid].zfinal,format='(f6.4)')
  if ~keyword_set(xrange_ccf) then xrange_ccf=[-20,80]
  if ~keyword_set(yrange_ccf) then yrange_ccf=[-0.5,1]
  plot,[0],[0],xrange=xrange_ccf,yrange=yrange_ccf,xtitle='Observed Time Delay [days]',ytitle='CCF', $
    pos=pos2,charsize=charsize,xticklen=ticklen*0.5,title=tag,/noerase
  oplot, ccf_tau, ccf0
  tau_mea = [cent_tau[0], quantile_1d(0.16, cent_tau), quantile_1d(0.84, cent_tau) ]
  oplot, [tau_mea[0],tau_mea[0]],yrange,color=cgcolor('blue')
  oplot, [tau_mea[1],tau_mea[1]],yrange,color=cgcolor('blue'),line=2
  oplot, [tau_mea[2],tau_mea[2]],yrange,color=cgcolor('blue'),line=2
  xyouts, pos2[0] + 0.02, pos2[3] - 0.08, /norm, 'peak sig='+string(peak_sig,format='(f0.3)'),charsize=charsize

  print, tau_mea
  fmt='(i3.3, " ", f6.4, " ", 3(f5.1, " "), f0.3, " ", i0, " ", f6.3)'
  tau_rest=tau_mea/(1. + rm_info[rmid].zfinal)
  tau_rest[1]=tau_rest[1]-tau_rest[0]
  tau_rest[2]=tau_rest[2]-tau_rest[0]
  ; compute a virial product (VP)
  logvp=alog10(3d10*tau_rest[0]*24.*3600.*(sig_rms*1d5)^2/6.6742d-8/1.989d33)
  if keyword_set(lun_tau) then printf, lun_tau,rmid,rm_info[rmid].zfinal, $
     tau_rest, peak_sig, round(sig_rms),logvp, format=fmt
  endplot

  ; now estimate the FWHM/dispersion from the rms broad-line spectra

end

; this is a slightly different version, which outputs 3 vertical panels of plots
pro rm_plot_lag_one2,rmid,topdir=topdir,ep_rej=ep_rej,line=line,conti=conti, $
   dt=dt,ccf_tau=ccf_tau,mirror=mirror,xrange_ccf=xrange_ccf,yrange_ccf=yrange_ccf, $
   dcf=dcf, lun_tau=lun_tau, lun_lc=lun_lc, vrange=vrange, wide=wide,tau_input=tau_input, $
   fl_err_min=fl_err_min,fc_err_min=fc_err_min, add_fl_const=add_fl_const, minerr=minerr

  if ~keyword_set(topdir) then topdir='/data3/quasar/yshen/work/lags/prepspec/'
  if ~keyword_set(conti) then conti='c5100'
  if ~keyword_set(line) then line='hb'


  if n_elements(minerr) eq 0 then minerr=0.03 ; default is to inflate LC errors to 3% fractional

  if n_elements(wide) eq 0 then wide=1

  common rm_info_block, rm_info
  if n_elements(rm_info) eq 0 then begin
    target_file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
    rm_info=mrdfits(target_file,1,/silent)
  endif

  rmtag='rm'+string(rmid,format='(i3.3)')
  ; get continuum LC
  lcfile=topdir+rmtag+'/'+rmtag+'_'+conti+'.dat'
  readcol,lcfile,format='d,d,d',t_conti,f_conti,e_conti,/silent
  nep_c=n_elements(t_conti)
  ; get line LC
  lcfile=topdir+rmtag+'/'+rmtag+'_'+line+'_t.dat'
  readcol,lcfile,format='d,d,d',t_line,f_line,e_line,/silent
  nep_l=n_elements(t_line)
  ; get rms sigma
  sigfile=topdir+rmtag+'/'+rmtag+'_vvblr.dat'
  readcol,sigfile,format='d,d,x,d,d,d,d',wave_line,Jrms_all, fwhm_rms_all,efwhm_all, sig_rms_all, esig_all,skip=1L
  if line eq 'hb' then wave_s=4862.68
  if line eq 'mg2' then wave_s=2798.75
  if line eq 'ha' then wave_s=6564.61
  if line eq 'he2_4686' then wave_s=4687.02
  ind_line=where(abs(wave_line - wave_s) le 1.)
  if ind_line[0] ne -1 then sig_rms=sig_rms_all[ind_line[0]] else sig_rms=0
  if ind_line[0] ne -1 then esig_rms=esig_all[ind_line[1]] else esig_rms=-1
  if ind_line[0] ne -1 then fwhm_rms=fwhm_rms_all[ind_line[0]] else fwhm_rms=0
  if ind_line[0] ne -1 then efwhm_rms=efwhm_all[ind_line[1]] else efwhm_rms=-1
  if ind_line[0] ne -1 then Jrms=Jrms_all[ind_line[0]] else Jrms=-1

  ; for avg spectrum
  avgfile=topdir+rmtag+'/'+rmtag+'_vblr.dat'
  readcol,avgfile,format='d,d,x,d,d,d,d',wave_line,Javg_all,fwhm_mean_all,efwhm_mean_all,sig_mean_all, esig_mean_all
  ind_line=where(abs(wave_line - wave_s) le 1.)
  if ind_line[0] ne -1 then Javg=Javg_all[ind_line[0]] else Javg=-1
  if ind_line[0] ne -1 then sig_mean=sig_mean_all[ind_line[0]] else sig_mean=0
  if ind_line[0] ne -1 then esig_mean=esig_mean_all[ind_line[1]] else esig_mean=-1
  if ind_line[0] ne -1 then fwhm_mean=fwhm_mean_all[ind_line[0]] else fwhm_mean=0
  if ind_line[0] ne -1 then efwhm_mean=efwhm_mean_all[ind_line[1]] else efwhm_mean=-1
  ; print, fwhm_mean, efwhm_mean
  ; pause

  ; convert the line LC to real units, i.e., adding mean flux and scale the amplitude
  ; in 1d-16 erg/s/cm2
  f_line = (Javg + Jrms * f_line)*1d16
  e_line = (Jrms * e_line)*1d16

  ; remove bad epoch if required
  flag=lonarr(nep_c)
  t_conti_all=t_conti & f_conti_all=f_conti & e_conti_all=e_conti
  t_line_all=t_line & f_line_all=f_line & e_line_all=e_line
  nep_c_all=n_elements(t_conti_all)
  nep_l_all=n_elements(t_line_all)
  if n_elements(ep_rej) ne 0 then begin
     flag[ep_rej-1]=1L
     ind=where(flag eq 0, complement=indd)
     if indd[0] ne -1 then begin ; get the bad epochs
       t_conti_bad=t_conti[indd] & f_conti_bad=f_conti[indd] & e_conti_bad=e_conti[indd]
       t_line_bad=t_line[indd] & f_line_bad=f_line[indd] & e_line_bad=e_line[indd]
     endif
     t_conti=t_conti[ind] & f_conti=f_conti[ind] & e_conti=e_conti[ind]
     t_line=t_line[ind] & f_line=f_line[ind] & e_line=e_line[ind]
     nep_c=n_elements(t_conti)
     nep_l=n_elements(t_line)
  endif
 ;if keyword_set(add_fl_const) then f_line=f_line+10.*stddev(f_line) ; assuming line variability is 10% of mean flux
  if keyword_set(minerr) then begin ; using inflated error on the LC, i.e, minerr=0.03 of the median error
      e_conti = e_conti > (minerr*median(f_conti) )
      if n_elements(ep_rej) ne 0 then e_conti_bad = e_conti_bad > (minerr*median(f_conti) )
      ;e_conti_all = e_conti_all > (minerr*median(f_conti) ) ; keep the original LC errors
      e_line = e_line > (minerr*median(f_line) )
      if n_elements(ep_rej) ne 0 then e_line_bad = e_line_bad > (minerr*median(f_line) )
      ;e_line_all = e_line_all > (minerr*median(f_line) ) ; keep the original LC errors
  endif


  ; output the LCs
  outdir=topdir+'output/'
  if ~keyword_set(lun_lc) then begin
    fmt='(f11.3, " ", e14.7, " ", e14.7, " ", i0)'
    lcfile=outdir+rmtag+'_'+conti+'.dat'
    openw,lun,lcfile,/get_lun
    for i=0L,nep_c_all-1 do printf,lun,t_conti_all[i],f_conti_all[i],e_conti_all[i], $
      flag[i],format=fmt
    close,lun
    free_lun,lun
    lcfile=outdir+rmtag+'_'+line+'.dat'
    openw,lun,lcfile,/get_lun
    for i=0L,nep_l_all-1 do printf,lun,t_line_all[i],f_line_all[i],e_line_all[i], $
      flag[i],format=fmt
    close,lun
    free_lun,lun
  endif else begin
    fmt='(i3.3, " ", a8, " ", f6.4, " ", f11.3, 2(" ", e14.7, " ", e14.7), " ", i0)'
    for i=0L,nep_c_all-1 do printf,lun_lc,rmid, line, rm_info[rmid].zfinal, t_conti_all[i],f_conti_all[i],e_conti_all[i],$
      f_line_all[i],e_line_all[i],flag[i],format=fmt
  endelse

  ; now compute the CCF
  if ~keyword_set(dt) then dt=2.
  if ~keyword_set(ccf_tau) then ccf_tau = -20. + (findgen(100./dt + 6L)*dt)
  ccf_tau_wide=-90. + (findgen(160./dt + 11L)*dt)
  ; -70. + (findgen(140./dt + 5L)*dt) ;; -10. + (findgen(70./dt + 5L)*dt)
  ; run xcorr_interp once to get the cent_tau and peak sig
  ccf0 = xcorr_interp(t_line,f_line,t_conti,f_conti,tau=ccf_tau,ndata=ndata0 $
            , /nogap,cent_tau=cent_tau0, peak_tau=peak_tau0 $
            , mirror=mirror,peak_sig=peak_sig)
  ccf0_wide=xcorr_interp(t_line,f_line,t_conti,f_conti,tau=ccf_tau_wide $
            , /nogap, mirror=mirror)
  acf0_wide=xcorr_interp(t_conti,f_conti,t_conti,f_conti,tau=ccf_tau_wide $
            , /nogap, /mirror)
  ; run again to get error estimates
  if ~keyword_set(tau_input) then begin ; 
    ccf = xcorr_interp(t_line,f_line,t_conti,f_conti,tau=ccf_tau,ndata=ndata $
            , /nogap,err1=e_line,err2=e_conti, /bootstrap,cent_tau=cent_tau, peak_tau=peak_tau $
            , mirror=mirror,ntrial=5000L,f1_err_min=fl_err_min,f2_err_min=fc_err_min)
  endif

  ; now generate a plot
  ; the map between line names and what to plot
  line1=['blr','ha','hb','mg2','he2_4686']
  figfile=outdir+rmtag+'_'+line+'_lag.eps'
  charsize=1.5 & ticklen=0.08
  begplot,name=figfile,/cmyk  ;,/encap  ,xsize=12.5, ysize=4
  line2=textoidl(['blr','H\alpha', 'H\beta', 'MgII', 'HeII'])
  indx_line=where(line1 eq line)
  ;pos11=[0.08,0.55,0.49,0.91] & pos2=[0.57, 0.16, 0.98, 0.91]
  ;pos12=[0.08,0.16,0.49,0.52]
  ;pos1 = [0.1, 0.70, 0.97, 0.965]
  pos11 = [0.13, 0.84, 0.96, 0.96]
  pos12 = [0.13, 0.71, 0.96, 0.83]
  pos2 = [0.13, 0.39, 0.96, 0.63]
  pos3 = [0.13, 0.07, 0.96, 0.31]

  xrange=[min(t_conti),max(t_conti)]
  yrange=[min(f_conti),max(f_conti)]
  symsize=0.6

  tag='J'+radec2sdssname(rm_info[rmid].ra,rm_info[rmid].dec) + $
    ', RMID='+string(rmid,format='(i3.3)') + ', z='+string(rm_info[rmid].zfinal,format='(f6.4)')
  plot,[0],[0],xrange=xrange,yrange=yrange, $
    title=tag, pos=pos11, charsize=charsize,xticklen=ticklen, $
    xtickname=replicate(' ', 6L), yticks=3, yminor=2
  oplot, [6650,6850], [median(f_conti), median(f_conti) ],line=1
  oploterror,t_conti,f_conti,e_conti,psym=symcat(9), symsize=symsize
  if n_elements(t_conti_bad) ne 0 then oploterror,[t_conti_bad],[f_conti_bad],[e_conti_bad],psym=symcat(9), $
    color=cgcolor('red'),errcolor=cgcolor('red'), symsize=symsize
  ;xyouts, pos11[2] - 0.11, pos11[1] + 0.02, conti, /norm,charsize=1.5
  xyouts, pos11[2] + 0.03, pos11[3] - 0.5*(pos11[3]-pos11[1]), conti, /norm,charsize=1.5, $
    orien=90, align=0.5
  yrange=[min(f_line),max(f_line)]
  plot,[0],[0],xrange=xrange,yrange=yrange,xtitle='MJD-50000 [days]', $
    pos=pos12,charsize=charsize,xticklen=ticklen,xtickformat='(i0)',/noerase,yticks=3, yminor=2
  oplot, [6650,6850], [median(f_line),median(f_line)],line=1
  oploterror,t_line,f_line,e_line,psym=symcat(9),color=cgcolor('blue'),symsize=symsize
  if n_elements(t_line_bad) ne 0 then oploterror,[t_line_bad],[f_line_bad],[e_line_bad],psym=symcat(9), $
    color=cgcolor('red'),errcolor=cgcolor('red'), symsize=symsize
  ;xyouts, pos12[2] - 0.11, pos12[1] + 0.02, line2[indx_line], /norm,charsize=1.5,color=cgcolor('blue')
  xyouts, pos12[2] + 0.03, pos12[3] - 0.5*(pos12[3]-pos12[1]), line2[indx_line], /norm,charsize=1.5,color=cgcolor('blue'), orien=90, align=0.5
  xyouts, 0.04, 0.8325, 'Relative Flux [arbitrary units]',/norm,orientation=90,align=0.5,charsize=1.5
  ; now plot the CCF
  ;tag='RMID='+string(rmid,format='(i3.3)') + ', z='+string(rm_info[rmid].zfinal,format='(f6.4)')
  if ~keyword_set(xrange_ccf) then xrange_ccf=[-20,90]
  if keyword_set(wide) then xrange_ccf=[-90,90]
  if ~keyword_set(yrange_ccf) then yrange_ccf=[-0.5,1]
  plot,[0],[0],xrange=xrange_ccf,yrange=yrange_ccf,xtitle='Observed Time Delay [days]',ytitle='CCF', $
    pos=pos2,charsize=charsize,xticklen=ticklen*0.5,/noerase,/xsty
  oplot, xrange_ccf, [0,0],line=1
  oplot, [0,0], yrange_ccf, line=1
  if ~keyword_set(wide) then oplot, ccf_tau, ccf0 else begin
     oplot, ccf_tau_wide, ccf0_wide
     oplot, ccf_tau_wide, acf0_wide, line=1, color=cgcolor('red')
     xyouts,0.45,pos2[3]-0.03,/norm, 'ACF',color=cgcolor('red')
  endelse
  if ~keyword_set(tau_input) then begin
     tau_mea = [cent_tau[0], quantile_1d(0.16, cent_tau), quantile_1d(0.84, cent_tau) ]
  endif else begin
     ; make sure the lower error is before the upper error
     if tau_input[1] gt tau_input[2] then tau_mea=[tau_input[0], tau_input[0]+tau_input[2], tau_input[0]+tau_input[1]] else tau_mea=[tau_input[0], tau_input[0]+tau_input[1], tau_input[0]+tau_input[2]]
  endelse
  oplot, [tau_mea[0],tau_mea[0]],yrange_ccf,color=cgcolor('blue')
  oplot, [tau_mea[1],tau_mea[1]],yrange_ccf,color=cgcolor('blue'),line=2
  oplot, [tau_mea[2],tau_mea[2]],yrange_ccf,color=cgcolor('blue'),line=2
  xyouts, pos2[0] + 0.02, pos2[3] - 0.04, /norm, 'peak sig='+string(peak_sig,format='(f0.3)'),charsize=1.5

  print, tau_mea
  fmt='(i3.3, " ",a8," ", f6.4, " ", 3(f5.1, " "), f0.3, " ", 8(i5, " "), 3(" ", f6.3) )'
  tau_rest=tau_mea/(1. + rm_info[rmid].zfinal)
  tau_rest[1]=tau_rest[1]-tau_rest[0]
  tau_rest[2]=tau_rest[2]-tau_rest[0]
  print, tau_rest, format='(f0.1, " ", f0.1, " ", f0.1)'
  ; compute a virial product (VP)
  logvp=alog10(3d10*tau_rest[0]*24.*3600.*(sig_rms*1d5)^2/6.6742d-8/1.989d33)
  ; compute error bars on logvp, using 2-sided err for tau and 1-sided err for sig_rms
  elogvp_lo=sqrt( (tau_rest[1]/alog(10.D)/tau_rest[0])^2 + (2.*esig_rms/alog(10.D)/sig_rms)^2 )
  elogvp_hi=sqrt( (tau_rest[2]/alog(10.D)/tau_rest[0])^2 + (2.*esig_rms/alog(10.D)/sig_rms)^2 )
  if keyword_set(lun_tau) then printf, lun_tau,rmid, line, rm_info[rmid].zfinal, $
     tau_rest, peak_sig, round(sig_rms), round(esig_rms), round(sig_mean), round(esig_mean), round(fwhm_rms), round(efwhm_rms), round(fwhm_mean), round(efwhm_mean), logvp,-elogvp_lo,elogvp_hi, format=fmt

  ; now plot the mean and rms broad line spectrum
  rmsspecfile=topdir+rmtag+'/'+rmtag+'_'+line+'_w.dat'
  readcol,rmsspecfile,format='d,d,d', wave, flux, err
  ; force the rms profile errors outside the window to be zero
  ind_nonzero = where(abs(flux) gt 1d-6 )
  if ind_nonzero[0] ne -1 then begin
    err[0:(min(ind_nonzero)-1)] = -1.
    err[max(ind_nonzero)+1:*] = -1.
  endif
  if line eq 'hb' then wrange=[4750, 4950]
  if line eq 'mg2' then wrange=[2750,2850]
  if line eq 'ha' then wrange=[6500,6610]
  if line eq 'he2_4686' then wrange=[4620,4740]
  if ~keyword_set(vrange) then vrange=[-8000,8000] ; velocity
  ang=string(197B)
  cspeed=2.9979246d5
  velo=(alog(wave)-alog(wave_s))*cspeed
  plot, velo, flux, xrange=vrange, xticklen=ticklen*0.5,/noerase, /xsty, $
    xtitle=textoidl('Velocity [km s^{-1}]'), $  ; xtitle='Rest Wavelength ['+ang+']', $
    ytitle='RMS Spectrum',pos=pos3, charsize=charsize ;, yticks=4, yminor=2
  oplot, velo, err, color=cgcolor('red'),line=2  

  xyouts, pos3[2] - 0.11, pos3[3] - 0.04, line2[indx_line], /norm,charsize=1.5,color=cgcolor('blue')

  endplot

end
