; Tests for looking at the velocity dispersion fits for plate 406.
; Loop from the highest-S/N galaxies to the lowest.
pro vdisptest, plate, fiberid, mjd=mjd, doplot=doplot, debug=debug, $
 eigenfile=eigenfile, columns=columns, $
 brightest=brightest, slowplot=slowplot, plotfile=plotfile1

   if (NOT keyword_set(plate)) then plate = 406
   if (NOT keyword_set(fiberid)) then fiberid = 1 + lindgen(640)
   if (NOT keyword_set(eigenfile)) then eigenfile='spEigenElodie.fits'
   if (NOT keyword_set(columns)) then columns=lindgen(24)

;   if (n_elements(slowplot) EQ 0) then slowplot = 1
;debug = 1
;doplot = 1
;brightsort = 1
;plate = 406
;mjd = 52238
;slowplot = 1

   stime0 = systime(1)

   readspec, plate, fiberid, mjd=mjd, flux=objflux, invvar=objivar, $
    synflux=synflux, zans=zans
   readspec, plate, mjd=mjd, fiberid[0], loglam=loglam ; same for every object

   if (keyword_set(plotfile1)) then begin
      if (size(plotfile1,/tname) EQ 'STRING') then plotfile = plotfile1 $
       else plotfile = string(plate, mjd, $
        format='("vdisp-",i4.4,"-",i5.5,".ps")')
      dfpsplot, plotfile, /color
   endif
   if (keyword_set(debug) OR keyword_set(plotfile)) then doplot = 1
   if (keyword_set(plotfile)) then debug = 0

   if (keyword_set(brightest)) then begin
      igal = where(strtrim(zans.class,2) EQ 'GALAXY' AND zans.zwarning EQ 0, $
       ngal)
      igal = igal[ reverse(sort(zans[igal].sn_median)) ] ; Sort by S/N
   endif else begin
      ngal = (n_elements(zans))
      igal = lindgen(ngal)
   endelse

   for jj=0, ngal-1 do begin
      print, 'Working on object ', jj+1, ' of ', ngal

      plottitle = string(plate, zans[igal[jj]].fiberid, $
       format='("Vel. Disp. Plate ",i4," Fiber ",i3)')

      restwave = 10^loglam / (1 + zans[igal[jj]].z)
      vdans1 = vdispfit(objflux[*,igal[jj]], objivar[*,igal[jj]], $
       loglam, zobj=zans[igal[jj]].z, $
       eigenfile=eigenfile, eigendir=eigendir, columns=columns, $
       yfit=yfit1, plottitle=plottitle, doplot=doplot, debug=debug)
      if (jj EQ 0) then begin
         vdans = vdans1
         yfit = yfit1
      endif else begin
         vdans = [[vdans],[vdans1]]
         yfit = [[yfit],[yfit1]]
      endelse

      print,'Fiber = ', zans[igal[jj]].fiberid, ' sigma=', vdans1.vdisp, $
       ' +/- ', vdans1.vdisp_err
      print, 'Chi^2/DOF = ', vdans1.vdispchi2 / (vdans1.vdispdof > 1)

      ipix = where(objivar[*,igal[jj]] GT 0 AND yfit1 NE 0, npix)
      if (keyword_set(doplot) AND npix GT 1) then begin
         !p.multi = 0
         csize = 2
         ymax = 1.25 * max(djs_median(objflux[ipix,igal[jj]],width=101))
         ymin = -0.2 * ymax
         djs_plot, [restwave[ipix]], [objflux[ipix,igal[jj]]], $
          yrange=[ymin,ymax], $
          color='default', xtitle='Rest-frame Wavelength', ytitle='Flux', $
          title=plottitle, charsize=csize
         djs_oplot, [restwave[ipix]], [yfit1[ipix]], color='red'
         djs_oplot, !x.crange, [0,0]
         djs_oplot, restwave[ipix], $
          objflux[ipix,igal[jj]]-synflux[ipix,igal[jj]], $
          color='blue'
         djs_oplot, restwave[ipix], objflux[ipix,igal[jj]]-yfit1[ipix], $
          color='red'

         xplot = 0.9 * !x.crange[0] + 0.1 * !x.crange[1]
         yplot = 0.1 * !y.crange[0] + 0.9 * !y.crange[1]
         djs_xyouts, xplot, yplot, string(vdans1.vdisp, vdans1.vdisp_err, $
          format='("\sigma = ", f4.0, " +/- ", f4.0)'), charsize=csize
         yplot = 0.2 * !y.crange[0] + 0.8 * !y.crange[1]
         djs_xyouts, xplot, yplot, $
          string(vdans1.vdispchi2 / (vdans1.vdispdof > 1), $
          format='("\chi^2/DOF = ", f4.2)'), charsize=csize

         if (keyword_set(debug)) then begin
            print, 'Press any key...'
            cc = strupcase(get_kbrd(1))
         endif

         chivec1 = (objflux[ipix,igal[jj]]-yfit1[ipix]) $
          * sqrt(objivar[ipix,igal[jj]])
         chivec2 = (objflux[ipix,igal[jj]]-synflux[ipix,igal[jj]]) $
          * sqrt(objivar[ipix,igal[jj]])
         binsz = 0.1
         !p.multi = 0
         plothist, chivec1, bin=binsz, xrange=[-10,10], $
          xtitle=textoidl('\chi'), ytitle=textoidl('\chi Distribution'), $
          title=plottitle, charsize=csize
         plothist, chivec2, bin=binsz, /overplot, color=djs_icolor('blue')
         xplot = !x.crange[0] + findgen(101)*(!x.crange[1]-!x.crange[0])/100
         yplot = exp(-0.5*xplot^2) * npix * binsz / sqrt(2.*!pi)
         djs_oplot, xplot, yplot, color='red'

         if (keyword_set(debug)) then begin
            print, 'Press any key...'
            cc = strupcase(get_kbrd(1))
         endif
      endif

      ; Now fit for a different number of stellar eigentemplates
      if (keyword_set(slowplot)) then begin
         neigen = 50
         vdmany = 0
         for ieigen=1, neigen do begin
            ; Burles counter...
            print, format='("Num eigen ",i5," of ",i5,a1,$)', $
             ieigen, neigen, string(13b)

            vdans1 = vdispfit(objflux[*,igal[jj]], objivar[*,igal[jj]], $
             loglam, zobj=zans[igal[jj]].z, $
             eigenfile=eigenfile, eigendir=eigendir, columns=lindgen(ieigen))
            if (NOT keyword_set(vdmany)) then vdmany = vdans1 $
             else vdmany = [vdmany,vdans1]
         endfor
         if (NOT keyword_set(vdall)) then vdall = vdmany $
          else vdall = [[vdall],[vdmany]]

         !p.multi = [0,1,3]
         !x.range = 0
         !y.range = 0
         csize = 2
         djs_plot, lindgen(neigen)+1, vdmany.vdisp, /ynozero, $
          xtitle='Number of eigenspectra', ytitle='Vel-Disp. [km/s]', $
          charsize=csize, title=plottitle
         oploterr, lindgen(neigen)+1, vdmany.vdisp, vdmany.vdisp_err
         djs_plot, lindgen(neigen)+1, vdmany.vdispchi2, /ynozero, psym=-4, $
          xtitle='Number of eigenspectra', ytitle='\chi^2 of vel-disp.', $
          charsize=csize
         chi2diff = vdmany[0:neigen-2].vdispchi2 - vdmany[1:neigen-1].vdispchi2
         ymax = max(chi2diff[2:neigen-2])
         djs_plot, lindgen(neigen)+1, chi2diff, /ynozero, psym=-4, $
          yrange=[-0.1*ymax,1.1*ymax], $
          xtitle='Number of eigenspectra', ytitle='\Delta(\chi^2)', $
          charsize=csize
         djs_oplot, !x.crange, [0,0]
         !p.multi = 0

         if (keyword_set(debug)) then begin
            print, 'Press any key...'
            cc = strupcase(get_kbrd(1))
         endif
      endif
   endfor

   splog, 'Total time = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'

   if (keyword_set(plotfile)) then dfpsclose
stop

   return
end
