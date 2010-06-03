; Make refraction-color plots for QSOs with real data.
pro qsorplot

   columns = ['Z','ZWARNING','CLASS','PRIMTARGET','PLUG_RA','PLUG_DEC', $
    'OFFSETRA','OFFSETDEC','MODELFLUX','MODELFLUX_IVAR', $
    'SPECTROSYNFLUX','SPECTROSYNFLUX_IVAR']
   spall = hogg_mrdfits(filepath('spAll.fits', $
    root_dir=getenv('BOSS_SPECTRO_REDUX')), 1, columns=columns, nchunk=20)
   indx = where(strmatch(spall.class,'QSO*') $
    AND spall.zwarning EQ 0 $
    AND (spall.primtarget AND 2L^0 + 2L^1 + 2L^2 +2L^3 + 2L^4) NE 0 $
    AND spall.plug_dec GT -7.5 AND spall.plug_dec LT 7.5)
   spall = spall[indx]

   qsorefract, ztab=ztab, urefract=urefracttab, grefract=grefracttab, $
    ugcolor=ugcolortab, grcolor=grcolortab, ricolor=ricolortab, $
    izcolor=izcolortab
   ; Normalize refraction to 30 deg from zenith, e.g. on the equator
   urefracttab = tan(30./!radeg) * urefracttab
   grefracttab = tan(30./!radeg) * grefracttab

   modelmag = 22.5 - 2.5*alog10(spall.modelflux>0.1)
   ugcolor = modelmag[0,*] - modelmag[1,*]
   grcolor = modelmag[1,*] - modelmag[2,*]
   ricolor = modelmag[2,*] - modelmag[3,*]
   izcolor = modelmag[3,*] - modelmag[4,*]

   synthmag = 22.5 - 2.5*alog10(spall.spectrosynflux>0.1)
   grspectro = synthmag[1,*] - synthmag[2,*]
   rispectro = synthmag[2,*] - synthmag[3,*]
   gispectro = synthmag[1,*] - synthmag[3,*]

   rref = 0.5 * (spall.offsetdec[2] + spall.offsetdec[3])
   urefract = spall.offsetdec[0] - rref
   grefract = spall.offsetdec[1] - rref

   dz = 0.1
   nbin = 50
   zlo = 0.0 + findgen(nbin) * dz
   zhi = zlo + dz
   zmid = 0.5 * (zlo + zhi)
   csize = 2.0

   dfpsplot, 'qsorplot.ps', /color, /square

   urefract_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    urefract_hist[ibin] = median(urefract[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   djs_plot, spall.z, urefract, ps=3, xrange=[0,4], yrange=[-0.3,0.3], $
    xtitle='Redshift', ytitle='Refract-u', charsize=csize
   djs_oplot, zmid, urefract_hist, psym=-4, color='red'
   djs_oplot, ztab, urefracttab-1.1, color='cyan'

fitval = interpol(urefracttab,ztab,spall.z)
ii=where(spall.z lt 3)
print,djsig(urefract[ii])
print,djsig((urefract-fitval)[ii])

   grefract_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    grefract_hist[ibin] = median(grefract[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   djs_plot, spall.z, grefract, ps=3, xrange=[0,4], yrange=[-0.3,0.3], $
    xtitle='Redshift', ytitle='Refract-g', charsize=csize
   djs_oplot, zmid, grefract_hist, psym=-4, color='red'
   djs_oplot, ztab, grefracttab-0.42, color='cyan'

   ugcolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    ugcolor_hist[ibin] = median(ugcolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   djs_plot, spall.z, ugcolor, ps=3, xrange=[0,4], yrange=[-0.5,1.5], $
    xtitle='Redshift', ytitle='(u-g)', charsize=csize
   djs_oplot, zmid, ugcolor_hist, psym=-4, color='red'
   djs_oplot, ztab, ugcolortab+0.0, color='cyan'
djs_oplot, ztab, ugcolortab+0.15*(2.-ztab), color='green'

   grcolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    grcolor_hist[ibin] = median(grcolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   djs_plot, spall.z, grcolor, ps=3, xrange=[0,4], yrange=[-0.5,1.5], $
    xtitle='Redshift', ytitle='(g-r)', charsize=csize
   djs_oplot, zmid, grcolor_hist, psym=-4, color='red'
   djs_oplot, ztab, grcolortab+0.0, color='cyan'

   ricolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    ricolor_hist[ibin] = median(ricolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   djs_plot, spall.z, ricolor, ps=3, xrange=[0,4], yrange=[-0.5,1.5], $
    xtitle='Redshift', ytitle='(r-i)', charsize=csize
   djs_oplot, zmid, ricolor_hist, psym=-4, color='red'
   djs_oplot, ztab, ricolortab+0.0, color='cyan'

   izcolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    izcolor_hist[ibin] = median(izcolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   djs_plot, spall.z, izcolor, ps=3, xrange=[0,4], yrange=[-0.5,1.5], $
    xtitle='Redshift', ytitle='(i-z)', charsize=csize
   djs_oplot, zmid, izcolor_hist, psym=-4, color='red'
   djs_oplot, ztab, izcolortab+0.20, color='cyan'

   grdiff = grspectro - grcolor
   grsphist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    grsphist[ibin] = median(grdiff[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   djs_plot, spall.z, grdiff, psym=3, xrange=[0,5], yrange=[-1,1], charsize=csize
   djs_oplot, zmid, grsphist, psym=-4, color='red'
   djs_oplot, !x.crange, [0,0], color='cyan'

   ridiff = rispectro - ricolor
   risphist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    risphist[ibin] = median(ridiff[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   djs_plot, spall.z, ridiff, psym=3, xrange=[0,5], yrange=[-1,1], charsize=csize
   djs_oplot, zmid, risphist, psym=-4, color='red'
   djs_oplot, !x.crange, [0,0], color='cyan'

   gidiff = gispectro - (grcolor + ricolor)
   gisphist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    gisphist[ibin] = median(gidiff[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   djs_plot, spall.z, gidiff, psym=3, xrange=[0,5], yrange=[-1,1], charsize=csize
   djs_oplot, zmid, gisphist, psym=-4, color='red'
   djs_oplot, !x.crange, [0,0], color='cyan'

   dfpsclose
stop

end
