;------------------------------------------------------------------------------
pro zplot_circle, radius, label=label, ltheta=ltheta, _EXTRA=KeywordsForPlot

   if (n_elements(radius) GT 1) then begin
      for i=0, n_elements(radius)-1 do begin
         if (keyword_set(label)) then $
          zplot_circle, radius[i], label=label[i], ltheta=ltheta, $
           _EXTRA=KeywordsForPlot $
         else $
          zplot_circle, radius[i], ltheta=ltheta, _EXTRA=KeywordsForPlot
      endfor
      return
   endif

   nsamp = 100 ; Number of samples for half the circle
   theta = 180. * findgen(nsamp) / (nsamp * !radeg)
   xplot = radius[0] * cos(theta)
   yplot = radius[0] * sin(theta)
   xplot = [xplot, reverse(xplot), xplot[0]]
   yplot = [yplot, -yplot, yplot[0]]
   djs_oplot, xplot, yplot, _EXTRA=KeywordsForPlot

   if (keyword_set(label) AND NOT keyword_set(ltheta)) then ltheta = 90.
   if (n_elements(ltheta) EQ 1) then begin
      if (NOT keyword_set(label)) then begin
         label = strtrim(string(radius[0]),2)
         ; Trim trailing zeros from this string if they're after a decimal point
         if (strpos(label,'.') NE -1) then begin
            ipos = strlen(label) - 1
            while (ipos GT 0 AND strmid(label,ipos,1) EQ '0') do begin
               label = strmid(label,0,ipos)
               ipos = ipos - 1
            endwhile
         endif
      endif

      xplot = cos(ltheta / !radeg) $
       * (radius[0] + 0.01*(!x.crange[1]-!x.crange[0]))
      yplot = sin(ltheta / !radeg) $
       * (radius[0] + 0.01*(!y.crange[1]-!y.crange[0]))
      align = sin(0.5 * ltheta / !radeg)
      djs_xyouts, xplot, yplot, 'z=' + label, align=align
   endif

   return
end
;------------------------------------------------------------------------------
pro zplot_exclude, theta, ramid=ramid

   if (NOT keyword_set(ramid)) then ramid = 0
   nplot = 40
   maxrad = sqrt( (max(abs(!x.crange)))^2 + (max(abs(!y.crange)))^2 )

   xplot = cos((theta - ramid) / !radeg)
   yplot = sin((theta - ramid) / !radeg)

   for i=0, nplot-1 do begin
      rfac = i * maxrad / nplot
      djs_oplot, rfac*xplot, rfac*yplot, linestyle=2, _EXTRA=KeywordsForPlot
   endfor
   djs_oplot, [0, rfac*xplot[0]], [0, rfac*yplot[0]]
   djs_oplot, [0, rfac*xplot[1]], [0, rfac*yplot[1]]

   return
end
;------------------------------------------------------------------------------
pro zplot_exclude_galaxy, _EXTRA=extra

   ; The Galactic plane crosses the equatorial plane at RA= 102.86, 282.86
   zplot_exclude, 102.86 + [-15,15], _EXTRA=extra
   zplot_exclude, 282.86 + [-15,15], _EXTRA=extra

   return
end
;------------------------------------------------------------------------------
pro zplot_label, ramid=ramid

   xplot = [!x.crange[1], 0, !x.crange[0], 0]
   yplot = [0, !y.crange[1], 0, !y.crange[0]]
   for i=0, n_elements(xplot)-1 do begin
      theta = (atan(yplot[i], xplot[i]) * !radeg - ramid) / 15. ; in hours
      if (theta LT 0) then theta = theta + 24.
      if (xplot[i] EQ 0) then align = 0.5 $
       else if (xplot[i] GT 0) then align = 0.0 $
       else align = 1.0
      if (yplot[i] GT 0) then thisy = 1.01 * yplot[i] $
       else if (yplot[i] LT 0) then thisy = 1.04 * yplot[i] $
       else thisy = yplot[i]
      djs_xyouts, xplot[i], thisy, $
       ' '+strtrim(string(round(theta)),2)+'^h ', align=align
   endfor

   return
end
;------------------------------------------------------------------------------
pro zplot

   ramid = 0.0

   ; Exclude plates 202,260,326
;   plate = [265,266,267,268,269,270,271,274,275,276,277,278,279, $
;    280,281,282,283,284,285,291,292,293,297,298,299, $
;    300,301,302,303,304,305,306,307,308,309, $
;    310,311,312,313,314,315, $
;    338,339,340,341,342,343,344,345,346,348, $
;    373,374,379,380,381,382,383,384,385,386,387, $
;    388,389,390,391,392,393,394,395,396,397,398,399, $
;    400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416]
   platelist, plist=plist
   ii = where(abs(plist.dec) LT 0.25 AND plist.qsurvey $
    AND strtrim(plist.status1d,2) EQ 'Done')

   readspec, plist[ii].plate, mjd=plist[ii].mjd, plug=plug, zans=zans

   xplot = zans.z * cos((plug.ra - ramid) / !radeg)
   yplot = zans.z * sin((plug.ra - ramid) / !radeg)

   spec_gal = strtrim(zans.class,2) EQ 'GALAXY' AND zans.zwarning EQ 0
   spec_qso = strtrim(zans.class,2) EQ 'QSO' AND zans.zwarning EQ 0
   target_main = (plug.primtarget AND 2^6) NE 0 $
    OR (plug.primtarget AND 2^7) NE 0 $
    OR (plug.primtarget AND 2^8) NE 0
   target_brg = (plug.primtarget AND 2^5) NE 0 $
    OR (plug.primtarget AND 2^26) NE 0
   target_brg = target_brg AND (target_main EQ 0) ; Don't let an object be
                                                  ; both in MAIN and BRG.
   target_qso = (plug.primtarget AND 2^0) NE 0 $
    OR (plug.primtarget AND 2^1) NE 0 $
    OR (plug.primtarget AND 2^2) NE 0 $
    OR (plug.primtarget AND 2^3) NE 0 $
    OR (plug.primtarget AND 2^4) NE 0

   imain = where(spec_gal AND target_brg EQ 0, nmain)
   ibrg = where(spec_gal AND target_brg, nbrg)
   iqso = where(spec_qso, nqso)

   !x.margin = [4,4]
   !y.margin = [2,2]

   ;----------
   ; Plot to z=0.15

   dfpsplot, 'zplot-main.ps', /color, /square
   zmax = 0.151
;   !p.region = [-zmax,zmax,-zmax,zmax]
   plot, [0], [0], /nodata, xrange=[-zmax,zmax], yrange=[-zmax,zmax], $
    xstyle=1, ystyle=1, xticks=1, yticks=1, $
    xtickname=[' ',' '], ytickname=[' ',' '], $
    title=string(nmain, format='("EQUATORIAL STRIPE    (", i6, " galaxies)")')
   djs_oplot, xplot[imain], yplot[imain], ps=3
   djs_oplot, xplot[ibrg], yplot[ibrg], ps=3, color='red'
   djs_oplot, xplot[iqso], yplot[iqso], ps=3, color='blue'
   zplot_circle, [0.05, 0.10], ltheta=70
   zplot_circle, [0.15]
   zplot_exclude_galaxy, ramid=ramid
   zplot_label, ramid=ramid
   dfpsclose

   ;----------
   ; Plot to z=0.60

   dfpsplot, 'zplot-brg.ps', /color, /square
   zmax = 0.601
   plot, [0], [0], /nodata, xrange=[-zmax,zmax], yrange=[-zmax,zmax], $
    xstyle=1, ystyle=1, xticks=1, yticks=1, $
    xtickname=[' ',' '], ytickname=[' ',' '], $
    title=string(nbrg, format='("EQUATORIAL STRIPE    (", i6, " LRGs)")')
   djs_oplot, xplot[imain], yplot[imain], ps=3
   djs_oplot, xplot[ibrg], yplot[ibrg], ps=1, symsize=0.25, color='red'
   djs_oplot, xplot[iqso], yplot[iqso], ps=3, symsize=0.25, color='blue'
   zplot_circle, [0.20, 0.40], ltheta=70
   zplot_circle, [0.60]
   zplot_exclude_galaxy, ramid=ramid
   zplot_label, ramid=ramid
   dfpsclose

   ;----------
   ; Plot to z=5

   dfpsplot, 'zplot-qso.ps', /color, /square
   zmax = 5.01
   plot, [0], [0], /nodata, xrange=[-zmax,zmax], yrange=[-zmax,zmax], $
    xstyle=1, ystyle=1, xticks=1, yticks=1, $
    xtickname=[' ',' '], ytickname=[' ',' '], $
    title=string(nqso, format='("EQUATORIAL STRIPE    (", i6, " QSOs)")')
   djs_oplot, xplot[imain], yplot[imain], ps=3
   djs_oplot, xplot[ibrg], yplot[ibrg], ps=3, color='red'
   djs_oplot, xplot[iqso], yplot[iqso], ps=1, $
    symsize=zans[iqso].z/8., color='blue'
   zplot_circle, [1,2,3,4], ltheta=70
   zplot_circle, [5]
   zplot_exclude_galaxy, ramid=ramid
   zplot_label, ramid=ramid
   dfpsclose

   ;----------
   ; Plot to z=5

   zmin = 0.01
   logzmin = alog10(zmin)
   xplot = (alog10(zans.z>zmin)-logzmin) * cos((plug.ra - ramid) / !radeg)
   yplot = (alog10(zans.z>zmin)-logzmin) * sin((plug.ra - ramid) / !radeg)

   dfpsplot, 'zplot-log.ps', /color, /square
   logzmax = alog10(5.01)
   lzrange = logzmax - logzmin
   plot, [0], [0], /nodata, xrange=[-lzrange,lzrange], $
    yrange=[-lzrange,lzrange], $
    xstyle=1, ystyle=1, xticks=1, yticks=1, $
    xtickname=[' ',' '], ytickname=[' ',' '], $
    title=string(nmain+nbrg+nqso, format='("EQUATORIAL STRIPE    (", i6, " objects)")')
   djs_oplot, xplot[imain], yplot[imain], ps=3
   djs_oplot, xplot[ibrg], yplot[ibrg], ps=3, color='red'
   djs_oplot, xplot[iqso], yplot[iqso], ps=1, $
    symsize=zans[iqso].z/8., color='blue'
   zplot_circle, alog10([0.03,0.1,0.5,2])-logzmin, $
    label=['0.03', '0.1','0.5','2'], ltheta=45
   zplot_circle, alog10(5)-logzmin, label='5', ltheta=45
   zplot_exclude_galaxy, ramid=ramid
   zplot_label, ramid=ramid
   dfpsclose

   return
end
