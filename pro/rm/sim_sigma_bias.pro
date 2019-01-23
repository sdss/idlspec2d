; simulate quasar samples to demonstrate the selection biases
; start from the velocity dispersion function, assign BH mass,
; then quasar luminosity, then impose threshold luminosity, then compute
; SE mass

pro sim_sigma_bias,plot_true_bh=plot_true_BH, no_msigma=no_msigma

; the z=0 VDF from Bernardi et al. (2010), which is similar to the z=0.5
; VDF as shown in Bezanson et al. (2011). VDF may be different at z=1, but
; the data is sparse there. 

phi0=2.099 & sig0=113.78 & alpha=0.94 & beta1=1.85

; generate random points
ntrial = 10000000L
sig_trial = randomu(seed, ntrial)*480. + 20.
; compute probability
; 1) the function in Bernardi (2010)
prob = phi0*(sig_trial/sig0)^alpha*exp(-(sig_trial/sig0)^beta1)/gamma(alpha/beta1)*beta1/sig_trial
; 2) power-law with index alpha
 alpha_pl=-2. 
; prob=(sig_trial/200.)^alpha_pl
prob=prob/max(prob)
; now generate a random number between 0 and 1
ran=randomu(seed2, ntrial)
ind=where(prob ge ran, ngood)
sig_keep=sig_trial[ind]

; assigne BH using KH13, plus 0.3 dex scatter
sig_msig=0.44  ; 0.44 for all gal and 0.3 for elliptical in Gultekin09
logmbh = alog10(0.309*(sig_keep/200.)^4.38*1d9) + randomn(seed4, ngood)*sig_msig
; this is no intrinsic m-sigma relation
if keyword_Set(no_msigma) then $
  logmbh=8.6 + randomn(seed4, ngood)*sig_msig
mbh = 10.D^logmbh

; assign luminosity using an Eddington ratio distribution
sig_L=0.3 ; scatter in Eddington ratio distribution
mean_loglambda=-1.
lambda=10.D^(mean_loglambda + randomn(seed3, ngood)*sig_L )
Lbol = 1.26d38*mbh*lambda
Lbol_mean = 1.26d38*mbh*10.D^mean_loglambda ; the mean luminosity at fixed true BH mass
L5100=Lbol/10.
L5100_mean = Lbol_mean/10.
sig_FWHM=0.18
logFWHM = 0.5*(logmbh - 0.91 - 0.5*alog10(L5100_mean/1d44)) + sig_FWHM*randomn(seed5,ngood)
logmbh_vp06=0.91 + 0.5*alog10(L5100/1d44) + 2.*logFWHM

output={sigma:sig_keep, logmbh:logmbh, lambda:lambda, Lbol_mean:Lbol_mean, Lbol:Lbol, $
    L5100:L5100, logFWHM:logFWHM, logmbh_vp06:logmbh_vp06, phi0:phi0, sig0:sig0, $
    alpha:alpha, beta:beta1, alpha_pl:alpha_pl,sig_msig:sig_msig, sig_L:sig_L, sig_FWHM:sig_FWHM }
outfile='/data3/quasar/yshen/work/agn_host/sim/sim_sample.fits'
if keyword_set(no_msigma) then outfile='/data3/quasar/yshen/work/agn_host/sim/sim_sample_nomsigma.fits'
mwrfits, output, outfile, /create

sim=output
sigma=sim.sigma & l5100=sim.l5100 & logmbh=sim.logmbh & logmbh_vp06=sim.logmbh_vp06
; shuffle sigma by 0.1 dex
sigma=10.D^(alog10(sigma) + randomn(seed, ngood)*0.1 )
simsize=0.2
plot, sim.sigma, 10.D^sim.logmbh,xrange=[30,700],yrange=[1d6,1d11],/xlog,/ylog,/xsty,/ysty,psym=3
oplot, [10, 1000], 0.309*([10,1000]/200.)^4.38*1d9, line=2,thick=2,color=cgcolor('red')
iend=100
ind=where(alog10(L5100) gt 43.5)
oplot,sigma[ind[0:iend]],10.D^logmbh_vp06[ind[0:iend]],psym=symcat(16),color=cgcolor('green'),symsize=symsize 
ind=where(alog10(L5100) gt 44.)
oplot,sigma[ind[0:iend]],10.D^logmbh_vp06[ind[0:iend]],psym=symcat(16),color=cgcolor('blue'),symsize=symsize
ind=where(alog10(L5100) gt 44.5)
oplot,sigma[ind[0:iend]],10.D^logmbh_vp06[ind[0:iend]],psym=symcat(16),color=cgcolor('red'),symsize=symsize

if keyword_set(plot_true_BH) then begin
ind=where(alog10(L5100) gt 43.5)
oplot,sigma[ind[0:iend]],10.D^logmbh[ind[0:iend]],psym=symcat(9),color=cgcolor('green'),symsize=symsize
ind=where(alog10(L5100) gt 44.)
oplot,sigma[ind[0:iend]],10.D^logmbh[ind[0:iend]],psym=symcat(9),color=cgcolor('blue'),symsize=symsize
ind=where(alog10(L5100) gt 44.5)
oplot,sigma[ind[0:iend]],10.D^logmbh[ind[0:iend]],psym=symcat(9),color=cgcolor('red'),symsize=symsize
endif

xyouts, 0.15, 0.92, textoidl('logL_{5100}>43.5'),/norm,color=cgcolor('green'),charsize=1.5
xyouts, 0.15, 0.87, textoidl('logL_{5100}>44.0'),/norm,color=cgcolor('blue'),charsize=1.5
xyouts, 0.15, 0.82, textoidl('logL_{5100}>44.5'),/norm,color=cgcolor('red'),charsize=1.5

end

pro plot_sim_bias, no_msigma=no_msigma
font='helvetica'

file='/data3/quasar/yshen/work/agn_host/sim/sim_sample.fits'
if keyword_set(no_msigma) then file='/data3/quasar/yshen/work/agn_host/sim/sim_sample_nomsigma.fits'
sim=mrdfits(file,1)
sigma=sim.sigma & l5100=sim.l5100 & logmbh=sim.logmbh & logmbh_vp06=sim.logmbh_vp06
; shuffle sigma by 0.1 dex
ngood=n_elements(sigma)
sigma_mea=10.D^(alog10(sigma) + randomn(666, ngood)*0.1 )

lum_thre=[43.5,44., 44.5]

figfile='/data3/quasar/yshen/work/agn_host/figs/sim_bias.eps'
if keyword_set(no_msigma) then figfile='/data3/quasar/yshen/work/agn_host/figs/sim_bias_nomsigma.eps'
begplot,name=figfile, xsize=10,ysize=6.5,font=font
charsize=2
xtitle=textoidl('\sigma_* [km s^{-1}]')
;ytitle=textoidl('M_{BH,vir} (H\beta, VP06) [M')+sunsymbol()+']'
ytitle=textoidl('BH Mass [M')+sunsymbol()+']'
xrange=[30,700] & yrange=[1d6,3d10]

pos0=[0.09, 0.54, 0.39, 0.99]
ticklen=0.03

nuse=501L
; plot actual BH masses
tag=['A','B','C']
for i=0L, 2L do begin

  dx=pos0[2]-pos0[0] & dy=pos0[3]-pos0[1]
  pos=[ pos0[0] + dx*i, pos0[1], $
        pos0[2] + dx*i, pos0[3] ]

  if i eq 0 then plot,[0],[0],xrange=xrange,yrange=yrange,pos=pos,/nodata, $
     /noerase,xtickname=replicate(' ', 6L), /xsty, /ysty, /xlog,/ylog,charsize=1.5, $
    xticklen=ticklen, yticklen=ticklen else $
   plot,[0],[0],xrange=xrange,yrange=yrange,pos=pos,/nodata,charsize=1.5, $
     /noerase,xtickname=replicate(' ', 6L),ytickname=replicate(' ', 6L),/xsty, /ysty, /xlog,/ylog, $
    xticklen=ticklen, yticklen=ticklen
   oplot, [10, 1000], 0.309*([10,1000]/200.)^4.38*1d9, line=2,thick=5,color=cgcolor('black')


  ;djs_contourpts, sigma, 10.D^logmbh, /overplot, /nopoints

  ind=where(alog10(L5100) ge lum_thre[i] )
  xarr=alog10(sigma_mea[ind[0:nuse-1]]) & xsig=replicate(0.1, nuse)
  yarr=logmbh[ind[0:nuse-1]] & ysig=replicate(0.1,nuse)
  linmix_err, xarr,yarr, post, xsig=xsig, ysig=ysig
  oplot, sigma_mea[ind[0:nuse-1]], 10.D^logmbh[ind[0:nuse-1]],psym=symcat(9),color=cgcolor('blue')
  alpha1=median(post.alpha) & beta1=median(post.beta)
  oplot, [10., 1000.], 10.D^alpha1*[10.,1000.]^beta1, thick=5,color=cgcolor('blue')


  xyouts, 40, 1d10, textoidl('logL_{5100}>')+string(lum_thre[i],format='(f4.1)'),charsize=1.5
  if i eq 0 then xyouts, 40, 4d9, textoidl('M_{BH,true}'),charsize=1.5,color=cgcolor('blue')

  xyouts, 450, 2d6, tag[i], charsize=charsize
endfor
pos0=[0.09, 0.09, 0.39, 0.54]
; plot VP06 BH masses
tag=['D','E','F']
for i=0L, 2L do begin

  dx=pos0[2]-pos0[0] & dy=pos0[3]-pos0[1]
  pos=[ pos0[0] + dx*i, pos0[1], $
        pos0[2] + dx*i, pos0[3] ]

  if i eq 0 then plot,[0],[0],xrange=xrange,yrange=yrange,pos=pos,/nodata, $
     /noerase, /xsty, /ysty, /xlog,/ylog,charsize=1.5, $ 
    xticklen=ticklen, yticklen=ticklen else $
   plot,[0],[0],xrange=xrange,yrange=yrange,pos=pos,/nodata,charsize=1.5, $
     /noerase,ytickname=replicate(' ', 6L),/xsty, /ysty, /xlog,/ylog, $
    xticklen=ticklen, yticklen=ticklen
   oplot, [10, 1000], 0.309*([10,1000]/200.)^4.38*1d9, line=2,thick=5,color=cgcolor('black')
  ;djs_contourpts, sigma, 10.D^logmbh, /overplot, /nopoints

  ind=where(alog10(L5100) ge lum_thre[i] )
  xarr=alog10(sigma_mea[ind[0:nuse-1]]) & xsig=replicate(0.1, nuse)
  yarr=logmbh_vp06[ind[0:nuse-1]] & ysig=replicate(0.1,nuse)
  linmix_err, xarr,yarr, post, xsig=xsig, ysig=ysig

  oplot, sigma_mea[ind[0:nuse-1]], 10.D^logmbh_vp06[ind[0:nuse-1]],psym=symcat(9),color=cgcolor('red')
  alpha1=median(post.alpha) & beta1=median(post.beta)
  oplot, [10., 1000.], 10.D^alpha1*[10.,1000.]^beta1, thick=5,color=cgcolor('red')
  if i eq 0 then xyouts, 40, 1d10, textoidl('M_{BH,VP06}'),charsize=1.5,color=cgcolor('red')  

  xyouts, 450, 2d6, tag[i], charsize=charsize

endfor
xyouts, 0.5, 0.01, align=0.5, xtitle, /norm, charsize=charsize
xyouts, 0.04, 0.5, align=0.5, ytitle, /norm, charsize=charsize, orien=90


endplot


end
