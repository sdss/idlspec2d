;+
; NAME:
;   checksn
;
; PURPOSE:
;   Calculate S/N per fiber in 3 bangs (g,r,i) and output plot and text 
;
; CALLING SEQUENCE:
;   checksn, flux, err, plug, wave, expres=expres, noplot=noplot, title=title
;
; OPTIONAL INPUTS:
;   expres     - Expression to glob in findfile, i.e. 'spMerge*fits'
;
; OPTIONAL KEYWORDS:
;   noplot     - just calculate table, and skip plot
;   title      - Add title to top of plot
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   flux       - spectra of input files
;   err        - corresponding error 
;   plug       - plugmap structure
;   wave       - corresponding log10 wavelengths
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   readidlout
;   sdsshead()
;   splog
;   sxpar()
;   s_plotdev
;
; REVISION HISTORY:
;   15-Apr-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------

pro checksn, flux, err, plug, wave, expres=expres, noplot=noplot, title=title

   if (keyword_set(expres)) then $
    readidlout, flux, sig=err, plug=plug, wave=wave, expres=expres, files=files

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

    gfit = fltarr(2)
    rfit = fltarr(2)
    ifit = fltarr(2)

    nonsky = where(strtrim(plug.objtype,2) NE 'SKY' AND $
                   strtrim(plug.objtype,2) NE 'NA', nnonsky)
    spectroid = (plug.fiberid - 1) / 320
    gdiff = fltarr(nnonsky)
    rdiff = fltarr(nnonsky)
    idiff = fltarr(nnonsky)
 

    badg = where(allspec[nonsky].g LE 0.0)
    goodg = where(allspec[nonsky].g GT 0.0, ngoodg)
    s1 = where(spectroid[nonsky] EQ 0)
    s2 = where(spectroid[nonsky] EQ 1)

    badi = where(allspec[nonsky].r LE 0.0)
    goodr = where(allspec[nonsky].r GT 0.0, ngoodr)

    badr = where(allspec[nonsky].i LE 0.0)
    goodi = where(allspec[nonsky].i GT 0.0, ngoodi)

    if (ngoodg GT 0) then begin
      gfit = ladfit(plug[nonsky[goodg]].mag[1], $
               alog10(allspec[nonsky[goodg]].g), absdev=gabsdev)
      gdiff[goodg] = alog10(allspec[nonsky[goodg]].g) - $
                     poly(plug[nonsky[goodg]].mag[1], gfit)
      gstddev = stddev(gdiff[goodg])
    endif

    if (ngoodr GT 0) then begin
      rfit = ladfit(plug[nonsky[goodr]].mag[2], $
              alog10(allspec[nonsky[goodr]].r), absdev=rabsdev)
      rdiff[goodr] = alog10(allspec[nonsky[goodr]].r) - $
              poly(plug[nonsky[goodr]].mag[2], rfit)
      rstddev = stddev(rdiff[goodr])
    endif

    if (ngoodi GT 0) then begin
      ifit = ladfit(plug[nonsky[goodi]].mag[3], $
           alog10(allspec[nonsky[goodi]].i), absdev=iabsdev)
      idiff[goodi] = alog10(allspec[nonsky[goodi]].i) - $
           poly(plug[nonsky[goodi]].mag[3], ifit)
      istddev = stddev(idiff[goodi])
    endif
    

      xmin = 15.0
      xmax = 21.0
      xfit = [xmin, xmax]

      gfiducialfit = [6.12, -0.28]
      glimit = 10^poly(plug[nonsky].mag[1], gfiducialfit)
      glow = where(allspec[nonsky].g LT glimit $
         AND plug[nonsky].mag[1] LT xmax, nglow)

      rfiducialfit = [6.18, -0.28]
      rlimit = 10^poly(plug[nonsky].mag[2], rfiducialfit)
      rlow = where(allspec[nonsky].r LT rlimit $
         AND plug[nonsky].mag[2] LT xmax, nrlow)

      ifiducialfit = [6.05, -0.28]
      ilimit = 10^poly(plug[nonsky].mag[3], ifiducialfit)
      ilow = where(allspec[nonsky].i LT ilimit $
         AND plug[nonsky].mag[3] LT xmax, nilow)

      bright = where(gdiff GT 0.2 AND rdiff GT 0.2 AND idiff GT 0.2, nbright)

      magsmin = [16.00, 16.5, 17.0, 17.5, 18., 18.5, 19., 19.5, 20.0, 20.5]

      magsmax = magsmin + 0.5

      nmag = n_elements(magsmin)
      bin = { mags: 0.5*(magsmin + magsmax), gsn: fltarr(nmag,2), $
             rsn: fltarr(nmag,2), isn: fltarr(nmag,2)}


      for imag=0, nmag -1 do begin

        for j=0,1 do begin

          g = where(spectroid EQ j AND $
             plug.mag[1] GE magsmin[imag] AND plug.mag[1] LE magsmax[imag],ng)
          if (ng GT 0) then bin.gsn[imag,j] = median([allspec[g].g $
                 * 10^(-gfit[1] * (plug[g].mag[1]-bin.mags[imag]))])

  	  r = where(spectroid EQ j AND $
            plug.mag[2] GE magsmin[imag] AND plug.mag[2] LE magsmax[imag],nr)
          if (nr GT 0) then bin.rsn[imag,j] = median([allspec[r].r $
                 * 10^(-rfit[1] * (plug[r].mag[2]-bin.mags[imag]))])

	  i = where(spectroid EQ j AND $
            plug.mag[3] GE magsmin[imag] AND plug.mag[3] LE magsmax[imag],ni)
          if (ni GT 0) then bin.isn[imag,j] = median([allspec[i].i $
                 * 10^(-ifit[1] * (plug[i].mag[3]-bin.mags[imag]))])
        endfor
      endfor

    if (NOT keyword_set(noplot)) then begin
      oldmulti = !p.multi  
      !p.multi = [0,2,3,0,1]
      plot, plug.mag[1], allspec.g, /nodata, /ylog, $
           ychars=2.0, xcharsize = 0.001, xr=[xmin,xmax], $
           ymargin = [1,3], /xs, yr=[0.8,100], /ys

      djs_oplot, xfit, 10^poly(xfit, gfit)
      djs_oplot, xfit, 10^poly(xfit, gfiducialfit), color='blue', thick=5
      if (keyword_set(title)) then xyouts, xmin+0.5, 110.0, title, $
            charsize=2.0

      if (s1[0] NE -1) then djs_oplot, plug[nonsky[s1]].mag[1], $
           allspec[nonsky[s1]].g > 0.9, ps=1, syms=0.7, color='red'
      if (s2[0] NE -1) then djs_oplot, plug[nonsky[s2]].mag[1], $
           allspec[nonsky[s2]].g > 0.9, ps=4, syms=0.7

      xyouts, xmin+0.5, 3.0, string(format='(a,f5.2,f6.2,a)','log S/N = ', $
           gfit, " * g'")

      if (size(gabsdev, /tname) NE 'UNDEFINED') then $
      xyouts, xmin+0.5, 2.0, string(format='(a,f5.3)','Deviation: ', gabsdev)

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.6], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.91], ps=1, $
         syms=0.7, color='red'
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.62], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.9], $
         'Spec 1', color=djs_icolor('red')

      djs_oplot, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.6], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.85], ps=4, syms=0.7
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.62], $
         10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.84], 'Spec 2'

      gminutes = (bin.gsn[6,*]/10.0)^2 * 60.0
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.9], $
           string(format='(f5.1)', gminutes[0]), color=djs_icolor('red')
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.84], $
           string(format='(f5.1)', gminutes[1])


      plot, plug.mag[2], allspec.r, /nodata, /ylog, ytitle = 'S/N', $
           ychars=2.0, xcharsize = 0.001, xr=[xmin,xmax], ymargin = [1,0], $
           /xs, yr=[0.8,100], /ys

      djs_oplot, xfit, 10^poly(xfit, rfit)
      djs_oplot, xfit, 10^poly(xfit, rfiducialfit), color='blue', thick=5

      if (s1[0] NE -1) then djs_oplot, plug[nonsky[s1]].mag[2], $
           allspec[nonsky[s1]].r > 0.9, ps=1, syms=0.7, color='red'
      if (s2[0] NE -1) then djs_oplot, plug[nonsky[s2]].mag[2], $
           allspec[nonsky[s2]].r > 0.9, ps=4, syms=0.7

      xyouts, xmin+0.5, 3.0, string(format='(a,f5.2,f6.2,a)','log S/N = ', $
        rfit, " * r'")
      if (size(rabsdev, /tname) NE 'UNDEFINED') then $
      xyouts, xmin+0.5, 2.0, string(format='(a,f5.3)','Deviation: ', rabsdev)
      rminutes = (bin.rsn[6,*]/10.0)^2 * 60.0
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.9], $
           string(format='(f5.1)', rminutes[0]), color=djs_icolor('red')
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.84], $
           string(format='(f5.1)', rminutes[1]) 

      plot, plug.mag[3], allspec.i, /nodata, /ylog, $
           xtitle = 'Fiber Magnitude', ychars=2.0, xchars = 2.0, $
           xr=[xmin,xmax], ymargin = [4,0], /xs, yr=[0.8,100], /ys

      djs_oplot, xfit, 10^poly(xfit, ifit)
      djs_oplot, xfit, 10^poly(xfit, ifiducialfit), color='blue', thick=5

      if (s1[0] NE -1) then djs_oplot, plug[nonsky[s1]].mag[3], $
           allspec[nonsky[s1]].i > 0.9, ps=1, syms=0.7, color='red'
      if (s2[0] NE -1) then djs_oplot, plug[nonsky[s2]].mag[3], $
           allspec[nonsky[s2]].i > 0.9, ps=4, syms=0.7

      xyouts, xmin+0.5, 3.0, string(format='(a,f5.2,f6.2,a)','log S/N = ', $
        ifit, " * i'")

      if (size(iabsdev, /tname) NE 'UNDEFINED') then $
        xyouts, xmin+0.5, 2.0, string(format='(a,f5.3)','Deviation: ', iabsdev)
      iminutes = (bin.isn[5,*]/10.0)^2 * 52.0
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.9], $
           string(format='(f5.1)', iminutes[0]), color=djs_icolor('red')
      xyouts, [!x.crange[0] + (!x.crange[1] - !x.crange[0])*0.8], $
           10^[!y.crange[0] + (!y.crange[1] - !y.crange[0])*0.84], $
           string(format='(f5.1)', iminutes[1])

    if (ngoodg GT 0) then $
        s_plotdev, plug[nonsky[goodg]].xfocal, plug[nonsky[goodg]].yfocal, $
           gdiff[goodg], 3.0, ychars=2.0, xcharsize = 0.001, $
           ymargin=[1,3], cap=1.0 $
    else plot, plug.xfocal, plug.yfocal, ychars=2.0, xchars=0.001, $
            ymargin=[1,3], /nodata
    if (nglow GT 0) then djs_oplot, plug[nonsky[glow]].xfocal, $
           plug[nonsky[glow]].yfocal, ps=1

    if (ngoodr GT 0) then $
      s_plotdev, plug[nonsky[goodr]].xfocal, plug[nonsky[goodr]].yfocal, $
           rdiff[goodr], 3.0, ychars=2.0, xcharsize = 0.001, $
           ymargin=[1,0], cap=1.0 $
    else plot, plug.xfocal, plug.yfocal, ychars=2.0, xchars=0.001, $
            ymargin=[1,0], /nodata
    if (nrlow GT 0) then djs_oplot, plug[nonsky[rlow]].xfocal, $
            plug[nonsky[rlow]].yfocal, ps=1

    if (ngoodi GT 0) then $
      s_plotdev, plug[nonsky[goodi]].xfocal, plug[nonsky[goodi]].yfocal, $
           idiff[goodi], 3.0, ychars=2.0, xcharsize = 2.0, $
           ymargin=[4,0], cap=1.0 $
    else plot, plug.xfocal, plug.yfocal, ychars=2.0, xchars=2.0, $
            ymargin=[4,0], /nodata
    if (nilow GT 0) then djs_oplot, plug[nonsky[ilow]].xfocal, $
           plug[nonsky[ilow]].yfocal, ps=1

      !p.multi  = oldmulti

    endif

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

   splog, "               S/N g'             S/N r'              S/N i'",/noname
   splog, " MAG          1      2           1      2            1      2",/noname

   splog,transpose([[bin.mags], [bin.gsn], [bin.rsn], [bin.isn]]), $
     format='(f5.2,3x,f8.3,f8.3,3x,f8.3,f8.3,3x,f8.3,f8.3)' ,/noname

   splog, ' Estimated Minutes in G: ', transpose(gminutes)
   splog, ' Estimated Minutes in R: ', transpose(rminutes)
   splog, ' Estimated Minutes in I: ', transpose(iminutes)

   nomap  = where(plug.fiberid EQ -1, nnomap)
   splog, 'Total number of Fibers read in: ', nplug
   splog, 'Number of Fibers NOT mapped: ', nnomap
   splog, 'Number of NON-SKY or NON-NA Fibers: ', nnonsky
   splog, 'Number of good fibers: ', ngoodg, ngoodr, ngoodi
   splog, 'Number of low  fibers: ', nglow, nrlow, nilow


   splog, ' Possible bright fibers: '
   for i=0, nbright - 1 do begin
     nearby = -1
     if (plug[bright[i]].fiberid LE 320) then $
       nearby = where(abs(plug[nonsky[bright[i]]].fiberid - $
                   plug.fiberid) LT 6 AND plug.fiberid LE 320) $
     else $
       nearby = where(abs(plug[nonsky[bright[i]]].fiberid - $
                   plug.fiberid) LT 6 AND plug.fiberid GT 320) 

     whop = where(plug[nearby].mag[2] LE 16, nwhop)
     if (nwhop EQ 0) then $
       splog, plug[nonsky[bright[i]]],/noname $
     else splog, 'Fiber ', plug[nonsky[bright[i]]].fiberid, $
          ' is near bright fiber ', plug[whop].fiberid

    endfor


    if (n_elements(files) GT 0) then begin
      hdr = sdsshead(files[0])
     splog, sxpar(hdr,'PLATEID'), sxpar(hdr,'DATE-OBS'), $
            sxpar(hdr,'MJD'), sxpar(hdr, 'TILEID'), $
            sxpar(hdr,'RA'), sxpar(hdr,'DEC'), sxpar(hdr, 'EXPTIME'), $
            format='(i4,a11, i6, i7, f10.4, f11.4, i8)', /noname
    endif


return
end
