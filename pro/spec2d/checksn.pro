pro checksn, flux, err, plug, wave, expres=expres, noplot=noplot, title=title

    if (keyword_set(expres)) then $
      readidlout, flux, sig=err, plug=plug, wave=wave, expres=expres

    w = 10^wave

    npix = (size(flux,/dim))[0]
    nspec = (size(flux,/dim))[1]

    snall = flux*0.0 
    nonzero = where(err GT 0)
    if (nonzero[0] NE -1) then $
      snall[nonzero] = flux[nonzero]/err[nonzero]
   
;
;	Let's first calculate median in each band
; 

    onespec = { g : 0.0, r : 0.0, i : 0.0, $           ; S/N entries
                gn1 : 0.0, rn1 : 0.0, in1 : 0.0, $     ; 1st noise estimate (median of err)
                gn2 : 0.0, rn2 : 0.0, in2 : 0.0, $     ; 2nd noise estimate (meanabsdev)
                gpix: 0L, rpix: 0L, ipix: 0L}          ; Number of pixels used

    allspec = replicate(onespec, nspec)

    for j=0, nspec - 1 do begin

	g = where(w[*,j] GT 4000 AND w[*,j] LT 5500 AND err[*,j] GT 0, ng)
        r = where(w[*,j] GT 5600 AND w[*,j] LT 6900 AND err[*,j] GT 0, nr)
        i = where(w[*,j] GT 6910 AND w[*,j] LT 8500 AND err[*,j] GT 0, ni)

        if (ng GT 0) then begin 
            level = median(flux[g,j],11)

            allspec[j].g = median(snall[g,j])
            allspec[j].gn1 = meanabsdev(flux[g,j]-level)
            allspec[j].gn2 = median(err[g,j])
            allspec[j].gpix = ng
        endif
        if (nr GT 0) then begin 
            level = median(flux[r,j],11)
            allspec[j].r = median(snall[r,j])
            allspec[j].rn1 = meanabsdev(flux[r,j]-level)
            allspec[j].rn2 = median(err[r,j])
            allspec[j].rpix = nr
        endif
        if (ni GT 0) then begin 
            level = median(flux[i,j],11)
            allspec[j].i = median(snall[i,j])
            allspec[j].in1 = meanabsdev(flux[i,j]-level)
            allspec[j].in2 = median(err[i,j])
            allspec[j].ipix = ni
        endif
     endfor

    nplug = n_elements(plug)
    nonsky = where(strtrim(plug.objtype,2) NE 'SKY' AND $
                   strtrim(plug.objtype,2) NE 'NA', nnonsky)
    goodg = where(allspec[nonsky].g GT 0.0, ngoodg)
    goodr = where(allspec[nonsky].r GT 0.0, ngoodr)
    goodi = where(allspec[nonsky].i GT 0.0, ngoodi)

    gfit = ladfit(plug[nonsky[goodg]].mag[1], alog10(allspec[nonsky[goodg]].g), absdev=gabsdev)
    gdiff = alog10(allspec[nonsky[goodg]].g) - poly(plug[nonsky[goodg]].mag[1], gfit)
    gstddev = stddev(gdiff)

    rfit = ladfit(plug[nonsky[goodr]].mag[2], alog10(allspec[nonsky[goodr]].r), absdev=rabsdev)
    rdiff = alog10(allspec[nonsky[goodr]].r) - poly(plug[nonsky[goodr]].mag[2], rfit)
    rstddev = stddev(rdiff)

    ifit = ladfit(plug[nonsky[goodi]].mag[3], alog10(allspec[nonsky[goodi]].i), absdev=iabsdev)
    idiff = alog10(allspec[nonsky[goodi]].i) - poly(plug[nonsky[goodi]].mag[3], ifit)
    istddev = stddev(idiff)
    

      xmin = 15.0
      xmax = 21.0
      xfit = [xmin, xmax]

      gfiducialfit = [6.12, -0.28]
      glimit = 10^poly(plug[nonsky].mag[1], gfiducialfit)
      glow = where(allspec[nonsky].g LT glimit $
         AND plug[nonsky].mag[1] LT xmax, nglow)

      rfiducialfit = [6.22, -0.28]
      rlimit = 10^poly(plug[nonsky].mag[2], rfiducialfit)
      rlow = where(allspec[nonsky].r LT rlimit $
         AND plug[nonsky].mag[2] LT xmax, nrlow)

      ifiducialfit = [6.05, -0.28]
      ilimit = 10^poly(plug[nonsky].mag[3], ifiducialfit)
      ilow = where(allspec[nonsky].i LT ilimit $
         AND plug[nonsky].mag[3] LT xmax, nilow)

      bright = where(gdiff GT 0.2 AND rdiff GT 0.2 AND idiff GT 0.2, nbright)

    if (NOT keyword_set(noplot)) then begin
      oldmulti = !p.multi  
      !p.multi = [0,2,3,0,1]
      plot, plug[nonsky[goodg]].mag[1], allspec[nonsky[goodg]].g, ps=1, /ylog, $
           ychars=2.0, xcharsize = 0.001, xr=[xmin,xmax], ymargin = [1,3], /xs, yr=[1,100], /ys, $
           title=title
      djs_oplot, xfit, 10^poly(xfit, gfit), color='red'
      djs_oplot, xfit, 10^poly(xfit, gfiducialfit), color='blue', thick=3
      xyouts, xmin+0.5, 3.0, string(format='(a,f5.2,f6.2,a)','log S/N = ', gfit, " * g'")
      xyouts, xmin+0.5, 2.0, string(format='(a,f5.3)','Deviation: ', gabsdev)


      plot, plug[nonsky[goodr]].mag[2], allspec[nonsky[goodr]].r, ps=1, /ylog, $
           ychars=2.0, xcharsize = 0.001, xr=[xmin,xmax], ymargin = [1,0], /xs, yr=[1,100], /ys
      djs_oplot, xfit, 10^poly(xfit, rfit), color='red'
      djs_oplot, xfit, 10^poly(xfit, rfiducialfit), color='blue', thick=3
      xyouts, xmin+0.5, 3.0, string(format='(a,f5.2,f6.2,a)','log S/N = ', rfit, " * r'")
      xyouts, xmin+0.5, 2.0, string(format='(a,f5.3)','Deviation: ', rabsdev)

      plot, plug[nonsky[goodi]].mag[3], allspec[nonsky[goodi]].i, ps=1, /ylog, $
           xchars=2.0, ychars=2.0, xr=[xmin,xmax], ymargin = [4,0], /xs, yr=[1,100], /ys
      djs_oplot, xfit, 10^poly(xfit, ifit), color='red'
      djs_oplot, xfit, 10^poly(xfit, ifiducialfit), color='blue', thick=3
      xyouts, xmin+0.5, 3.0, string(format='(a,f5.2,f6.2,a)','log S/N = ', ifit, " * i'")
      xyouts, xmin+0.5, 2.0, string(format='(a,f5.3)','Deviation: ', iabsdev)

      s_plotdev, plug[nonsky[goodg]].xfocal, plug[nonsky[goodg]].yfocal, gdiff, 3.0, $
             ychars=2.0, xcharsize = 0.001, ymargin=[1,3]
      s_plotdev, plug[nonsky[goodr]].xfocal, plug[nonsky[goodr]].yfocal, rdiff, 3.0, $
             ychars=2.0, xcharsize = 0.001, ymargin=[1,0]
      s_plotdev, plug[nonsky[goodi]].xfocal, plug[nonsky[goodi]].yfocal, idiff, 3.0, $
             ychars=2.0, xcharsize = 2.0, ymargin=[4,0]
      !p.multi  = oldmulti

    endif
    magsmin = [16.00, 16.5, 17.0, 17.5, 18., 18.5, 19., 19.5, 20.0, 20.5]

    magsmax = magsmin + 0.5

    nmag = n_elements(magsmin)
    bin = { mags: 0.5*(magsmin + magsmax), gsn: fltarr(nmag,2), $
             rsn: fltarr(nmag,2), isn: fltarr(nmag,2)}

    spectroid = (plug.fiberid - 1) / 320

  for imag=0, nmag -1 do begin

    for j=0,1 do begin

	g = where(spectroid EQ j AND $
                  plug.mag[1] GE magsmin[imag] AND plug.mag[1] LE magsmax[imag],ng)
        if (ng GT 0) then bin.gsn[imag,j] = median([allspec[g].g])

	r = where(spectroid EQ j AND $
                  plug.mag[2] GE magsmin[imag] AND plug.mag[2] LE magsmax[imag],nr)
        if (nr GT 0) then bin.rsn[imag,j] = median([allspec[r].r])

	i = where(spectroid EQ j AND $
                  plug.mag[3] GE magsmin[imag] AND plug.mag[3] LE magsmax[imag],ni)
        if (ni GT 0) then bin.isn[imag,j] = median([allspec[i].i])

      endfor
    endfor

   print, "               S/N g'             S/N r'              S/N i'"
   print, " MAG          1      2           1      2            1      2"

   print,transpose([[bin.mags], [bin.gsn], [bin.rsn], [bin.isn]]), $
     format='(f5.2,3x,f8.3,f8.3,3x,f8.3,f8.3,3x,f8.3,f8.3)' 

   nomap  = where(plug.fiberid EQ -1, nnomap)
   print, 'Total number of Fibers read in: ', nplug
   print, 'Number of Fibers NOT mapped: ', nnomap
   print, 'Number of NON-SKY or NON-NA Fibers: ', nnonsky
   print, 'Number of good fibers: ', ngoodg, ngoodr, ngoodi
   print, 'Number of low  fibers: ', nglow, nrlow, nilow


   print, ' Possible bright fibers: '
   for i=0, nbright - 1 do begin
     nearby = -1
     if (plug[bright[i]].fiberid LE 320) then $
       nearby = where(abs(plug[bright[i]].fiberid - plug.fiberid) LT 6 AND $
                      plug.fiberid LE 320) $
     else $
       nearby = where(abs(plug[bright[i]].fiberid - plug.fiberid) LT 6 AND $
                      plug.fiberid GT 320) 

     whop = where(plug[nearby].mag[2] LE 16, nwhop)
     if (nwhop EQ 0) then $
       print, plug[bright[i]]  $
     else print, 'Fiber ', plug[bright[i]].fiberid, ' is near bright fiber ', plug[whop].fiberid 

    endfor
return
end

