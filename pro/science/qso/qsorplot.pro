; Make refraction-color plots for QSOs with real data.
pro qsorplot

   spall = mrdfits(filepath('spAll.fits', root_dir=getenv('SPECTRO_DATA')),1)
   indx = where(strmatch(spall.class,'QSO*') $
    AND spall.zwarning EQ 0 $
    AND spall.plug_dec GT -7.5 AND spall.plug_dec LT 7.5)
   spall = spall[indx]

   qsorefract, ztab=ztab, urefract=urefracttab, grefract=grefracttab, $
    ugcolor=ugcolortab, grcolor=grcolortab
   ; Normalize refraction to 30 deg from zenith, e.g. on the equator
   urefracttab = tan(30./!radeg) * urefracttab
   grefracttab = tan(30./!radeg) * grefracttab

   ugcolor = spall.counts_model[0] - spall.counts_model[1]
   grcolor = spall.counts_model[1] - spall.counts_model[2]
   rref = 0.5 * (spall.offsetdec[2] + spall.offsetdec[3])
   urefract = spall.offsetdec[0] - rref
   grefract = spall.offsetdec[1] - rref

   dz = 0.1
   nbin = 50
   zlo = 0.0 + findgen(nbin) * dz
   zhi = zlo + dz
   zmid = 0.5 * (zlo + zhi)
   urefract_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    urefract_hist[ibin] = median(urefract[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, urefract, ps=3, xrange=[0,4], yrange=[-0.3,0.3]
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
   splot, spall.z, grefract, ps=3, xrange=[0,4], yrange=[-0.3,0.3]
   soplot, zmid, grefract_hist, psym=-4, color='red'
   soplot, ztab, grefracttab-0.42, color='cyan'
stop

   ugcolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    ugcolor_hist[ibin] = median(ugcolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, ugcolor, ps=3, xrange=[0,4], yrange=[-0.5,2.5]
   soplot, zmid, ugcolor_hist, psym=-4, color='red'
   soplot, ztab, ugcolortab+0.40, color='cyan'
stop

   grcolor_hist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    grcolor_hist[ibin] = median(grcolor[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, grcolor, ps=3, xrange=[0,4], yrange=[-0.5,2.5]
   soplot, zmid, grcolor_hist, psym=-4, color='red'
   soplot, ztab, grcolortab+0.35, color='cyan'
stop

   grspectro = -2.5*alog10(spall.counts_synth[1] / spall.counts_synth[2])
   grdiff = grspectro - grcolor
   grsphist = fltarr(nbin)
   for ibin=0, nbin-1 do $
    grsphist[ibin] = median(grdiff[where(spall.z GE zlo[ibin] $
     AND spall.z LT zhi[ibin])])
   splot, spall.z, grdiff, psym=3, xrange=[0,5], yrange=[-1,1]
   soplot, zmid, grsphist, psym=-4, color='red'
   soplot, !x.crange, [0,0], color='cyan'


end
