; plot the S/N/pix as functions of mags

pro plot_sn_rm


  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)

  imag=(target.psfmag)[3,0:848]
  imag=reform(imag)
  

  ; get per epoch median S/N
  med_sn_ep = (target.med_SN)[*,0:848]
  ; get coadded median S/N
  fiber=indgen(849)+1
  rm_readspec, 0, fiber,mjd=56837,flux=flux,invvar=ivar
  sn = flux*sqrt(ivar)
  med_sn_coadd = median(sn, dim=1)

  fig=getenv('IDLRM_DIR')+'/misc/sn_plot.eps'
  begplot,name=fig,/color,xsize=6,ysize=5,/encap
  imag_more=fltarr(32, 849)
  for i=0L, 31 do imag_more[i,*]=imag
  
  plot, imag_more, med_sn_ep, psym=3,xtitle='SDSS i', ytitle=textoidl('Median S/N per pixel'), $
    pos=[0.15,0.12,0.95,0.95],yrange=[1,300],/ysty,/ylog
  oplot, imag, med_sn_coadd, psym=symcat(9),symsize=0.5,color=cgcolor('red')

  endplot
end
