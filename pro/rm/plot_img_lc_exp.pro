; plot the light curve for an example object
; this is for proposal purposes

pro plot_img_lc_exp

file='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/imaging/cfht_15A/rm028i_cfht.txt'
readcol,file,format='d,d', mjd_cfht_i, mag_cfht_i
file='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/imaging/cfht_15A/rm028g_cfht.txt'
readcol,file,format='d,d', mjd_cfht_g, mag_cfht_g

figfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/imaging/cfht_15A/lc_rm028.eps'
begplot,name=figfile,/color,xsize=6,ysize=4

plot, mjd_cfht_g, mag_cfht_g, psym=symcat(16),xtitle='MJD', ytitle='AB Magnitude', $
   yrange=[19.5, 18.5],/ysty,pos=[0.15, 0.15, 0.92, 0.98],xtickname=['56650', '56700','56750','56800','56850'], xrange=[56650, 56850],/xsty
oplot, mjd_cfht_i, mag_cfht_i, psym=symcat(16), color=cgcolor('red')

file='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/imaging/bok_photometry/bokrm_g.fits'
bok_g=mrdfits(file,1)
file='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/imaging/bok_photometry/bokrm_i.fits'
bok_i=mrdfits(file,1)

rmid=28
oplot, bok_g[rmid].mjd, bok_g[rmid].APERMAG7+0.05,psym=symcat(9)
oplot, bok_i[rmid].mjd, bok_i[rmid].APERMAG7+0.1,psym=symcat(9),color=cgcolor('red')


items=[' CFHT',' Bok']
legend, items[0], box=0,pos=[0.7,0.95],/norm, psym=symcat(16)
legend, items[1], box=0,pos=[0.7,0.89],/norm, psym=symcat(9)

xyouts, 0.2, 0.9, 'g band',/norm
xyouts, 0.4, 0.9, 'i band',/norm,color=cgcolor('red')

endplot

end
