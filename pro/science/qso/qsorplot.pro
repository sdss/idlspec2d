; Make refraction-color plots for QSOs with real data.
pro qsorplot

   spall = mrdfits(filepath('spAll.fits', root_dir=getenv('SPECTRO_DATA')),1)
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

   ugcolor = spall.counts_model[0] - spall.counts_model[1]
   grcolor = spall.counts_model[1] - spall.counts_model[2]
   ricolor = spall.counts_model[2] - spall.counts_model[3]
   izcolor = spall.counts_model[3] - spall.counts_model[4]
   rref = 0.5 * (spall.offsetdec[2] + spall.offsetdec[3])
   urefract = spall.offsetdec[0] - rref
   grefract = spall.offsetdec[1] - rref

   dz = 0.1
   nbin = 50
   zlo = 0.0 + findgen(nbin) * dz
   zhi = zlo + dz
   zmid = 0.5 * (zlo + zhi)
   csize = 2.0

   urefract_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    urefract_hist[ibin] = median(urefract[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, urefract, ps=3, xrange=[0,4], yrange=[-0.3,0.3], $
    xtitle='Redshift', ytitle='Refract-u', charsize=csize
   soplot, zmid, urefract_hist, psym=-4, color='red'
   soplot, ztab, urefracttab-1.1, color='cyan'
stop
fitval = interpol(urefracttab,ztab,spall.z)
ii=where(spall.z lt 3)
print,djsig(urefract[ii])
print,djsig((urefract-fitval)[ii])

   grefract_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    grefract_hist[ibin] = median(grefract[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, grefract, ps=3, xrange=[0,4], yrange=[-0.3,0.3], $
    xtitle='Redshift', ytitle='Refract-g', charsize=csize
   soplot, zmid, grefract_hist, psym=-4, color='red'
   soplot, ztab, grefracttab-0.42, color='cyan'
stop

   ugcolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    ugcolor_hist[ibin] = median(ugcolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, ugcolor, ps=3, xrange=[0,4], yrange=[-0.5,1.5], $
    xtitle='Redshift', ytitle='(u-g)', charsize=csize
   soplot, zmid, ugcolor_hist, psym=-4, color='red'
   soplot, ztab, ugcolortab+0.0, color='cyan'
soplot, ztab, ugcolortab+0.15*(2.-ztab), color='green'
stop

   grcolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    grcolor_hist[ibin] = median(grcolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, grcolor, ps=3, xrange=[0,4], yrange=[-0.5,1.5], $
    xtitle='Redshift', ytitle='(g-r)', charsize=csize
   soplot, zmid, grcolor_hist, psym=-4, color='red'
   soplot, ztab, grcolortab+0.0, color='cyan'
stop

   ricolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    ricolor_hist[ibin] = median(ricolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, ricolor, ps=3, xrange=[0,4], yrange=[-0.5,1.5], $
    xtitle='Redshift', ytitle='(r-i)', charsize=csize
   soplot, zmid, ricolor_hist, psym=-4, color='red'
   soplot, ztab, ricolortab+0.0, color='cyan'
stop

   izcolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    izcolor_hist[ibin] = median(izcolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, izcolor, ps=3, xrange=[0,4], yrange=[-0.5,1.5], $
    xtitle='Redshift', ytitle='(i-z)', charsize=csize
   soplot, zmid, izcolor_hist, psym=-4, color='red'
   soplot, ztab, izcolortab+0.20, color='cyan'
stop


   grspectro = -2.5*alog10(spall.counts_synth[1] / spall.counts_synth[2])
   grdiff = grspectro - grcolor
   grsphist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    grsphist[ibin] = median(grdiff[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, grdiff, psym=3, xrange=[0,5], yrange=[-1,1], charsize=csize
   soplot, zmid, grsphist, psym=-4, color='red'
   soplot, !x.crange, [0,0], color='cyan'
stop

   rispectro = -2.5*alog10(spall.counts_synth[2] / spall.counts_synth[3])
   ridiff = rispectro - ricolor
   risphist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    risphist[ibin] = median(ridiff[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, ridiff, psym=3, xrange=[0,5], yrange=[-1,1], charsize=csize
   soplot, zmid, risphist, psym=-4, color='red'
   soplot, !x.crange, [0,0], color='cyan'
stop

   gispectro = -2.5*alog10(spall.counts_synth[1] / spall.counts_synth[3])
   gidiff = gispectro - (grcolor + ricolor)
   gisphist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    gisphist[ibin] = median(gidiff[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, gidiff, psym=3, xrange=[0,5], yrange=[-1,1], charsize=csize
   soplot, zmid, gisphist, psym=-4, color='red'
   soplot, !x.crange, [0,0], color='cyan'
stop

end
