; perform some tests of the flux vectors in different scenario

pro tmp_testmratio

   ; this is where the spFrame* files are
   datadir = '/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/'

   ; this is the exposure to be tested
   objname = 'spFrame-b1-00171812.fits.gz'

   spframe_read, datadir + objname, loglam=loglam1

   ; read in the sset for different cases
   sset_file1 = datadir + 'spFluxcalib-b1-00171812.fits.gz'
   sset1 = mrdfits(sset_file1,1)
   kindx1 = mrdfits(sset_file1,2)

   sset_file2 = datadir + 'recalib/test1/old/spFluxcalib-b1-00171812.fits.gz'
   sset2 = mrdfits(sset_file2,1)
   kindx2 = mrdfits(sset_file2,2)

   sset_file3 = datadir + 'recalib/test1/spFluxcalib-b1-00171812.fits.gz'
   sset3 = mrdfits(sset_file3,1)

   sset_file4 = datadir + 'recalib/test2/old/spFluxcalib-b1-00171812.fits.gz'
   sset4 = mrdfits(sset_file4,1)

   xrange=[3500,6300] & yrange=[0,50]
   plot, [0],[0],/nodata,xtitle='Wavelength',xrange=xrange,yrange=yrange,/xsty,/ysty

   ind = where(kindx1.qgood eq 1,ngood)
   ; for i=0L, ngood - 1L do oplot, 10.^kindx1[ind[i]].loglam,kindx1[ind[i]].mratio $
   ;    , psym=3
   
   ind = where(kindx2.qgood eq 1, ngood)
   for i=0L, ngood - 1L do oplot, 10.^kindx2[ind[i]].loglam,kindx2[ind[i]].mratio $
       , color=fsc_color('red'), psym=3

   oplot, 10.D^loglam1[*,0], bspline_valu(loglam1[*,0], sset1)
   oplot, 10.D^loglam1[*,0], bspline_valu(loglam1[*,0], sset2) $
       , color=fsc_color('cyan')


   oplot, 10.D^loglam1[*,0], bspline_valu(loglam1[*,0], sset3) $
       , color=fsc_color('green')

   oplot, 10.D^loglam1[*,0], bspline_valu(loglam1[*,0], sset4) $
       , color=fsc_color('magenta')

end

pro tmp_testmratio1

   ; this is where the spFrame* files are
   datadir = '/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/'

   ; this is the exposure to be tested
   objname = 'spFrame-b1-00171812.fits.gz'

   spframe_read, datadir + objname, loglam=loglam1

   sset_file = datadir + 'recalib/test2/spFluxcalib-b1-00171812.fits.gz'
   sset = mrdfits(sset_file,1)
   kindx = mrdfits(sset_file,2)
 
   colors=fsc_color(['opposite','green', 'red', 'cyan', 'magenta'])

   ind=where(kindx.qgood eq 1, ngood)

   xrange=[3500,6300] & yrange=[0,50]
   plot, [0],[0],/nodata,xtitle='Wavelength',xrange=xrange,yrange=yrange,/xsty,/ysty
   for i=0L,ngood - 1L do oplot, 10.^kindx[ind[i]].loglam,kindx[ind[i]].mratio $
       , color=colors[i mod 5], psym=3


   ; 
   loglam = kindx[ind].loglam
   mratio = kindx[ind].mratio
   mrativar = kindx[ind].mrativar
   ;iprox = indgen(ngood) ; [13,32,2,10]
   iprox = [15,11,12,20]

   nblue = 1L & nprox = n_elements(iprox)
   everyn = nblue * nprox * 10
   thisset = spflux_bspline(loglam[*,iprox], $
    mratio[*,iprox], mrativar[*,iprox], $
    everyn=everyn, outmask=mask_b_new)

   oplot, 10.D^loglam[*,0], bspline_valu(loglam[*,0], thisset, x2=x2),color=fsc_color('red') $
    , thick = 4

   oplot, 10.D^loglam[*,0], bspline_valu(loglam[*,0], sset, x2=x2),color=fsc_color('cyan')

end
