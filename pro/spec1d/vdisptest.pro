; Tests for looking at the velocity dispersion fits for plate 406.
; Loop from the highest-S/N galaxies to the lowest.
pro vdisptest, plate, mjd=mjd, doplot=doplot, debug=debug

eigenfile='spEigenElodie.fits'
columns=lindgen(24)
debug = 1
doplot = 1
brightsort = 1
plate = 406
mjd = 52238
slowplot = 1

   if (keyword_set(debug)) then doplot = 1

   stime0 = systime(1)

   readspec, plate, mjd=mjd, flux=objflux, invvar=objivar, $
    synflux=synflux, zans=zans
   readspec, plate, mjd=mjd, 1, loglam=loglam ; same for every object
   if (keyword_set(brightsort)) then begin
      igal = where(strtrim(zans.class,2) EQ 'GALAXY' AND zans.zwarning EQ 0, $
       ngal)
      igal = igal[ reverse(sort(zans[igal].sn_median)) ] ; Sort by S/N
   endif else begin
      ngal = (n_elements(zans))
      igal = lindgen(ngal)
   endelse

   for jj=0, ngal-1 do begin
      print, 'Working on object ', jj+1, ' of ', ngal

      wave = 10^loglam / (1 + zans[igal[jj]].z)
      vdans1 = vdispfit(objflux[*,igal[jj]], objivar[*,igal[jj]], $
       loglam, zobj=zans[igal[jj]].z, $
       eigenfile=eigenfile, eigendir=eigendir, columns=columns, $
       yfit=yfit1, doplot=doplot, debug=debug)
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
         ymax = 1.25 * max(djs_median(objflux[ipix,igal[jj]],width=101))
         ymin = -0.2 * ymax
         djs_plot, [wave[ipix]], [objflux[ipix,igal[jj]]], yrange=[ymin,ymax], $
          color='default', xtitle='Rest-frame Wavelength', ytitle='Flux'
         djs_oplot, [wave[ipix]], [yfit1[ipix]], color='red'
         djs_oplot, !x.crange, [0,0]
         djs_oplot, wave[ipix], objflux[ipix,igal[jj]]-synflux[ipix,igal[jj]], $
          color='blue'
         djs_oplot, wave[ipix], objflux[ipix,igal[jj]]-yfit1[ipix], $
          color='red'

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
          xtitle='\chi', ytitle=textoidl('\chi Distribution')
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
         !p.multi = [0,1,3]
         !x.range = 0
         !y.range = 0
         csize = 2
         djs_plot, lindgen(neigen)+1, vdmany.vdisp, /ynozero, $
          xtitle='Number of eigenspectra', ytitle='Vel-Disp. [km/s]', $
          charsize=csize
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

stop
   return
end
