; This code analyzes the focal plane mapping data taken for plate 174
; on MJD 51779.  See sdss-commish/265 for the posted results.
;------------------------------------------------------------------------------
pro focal_binplot, xdata, ydata, _EXTRA=EXTRA

   binsz = 50
   nbin = 650 / binsz
   xmin = -325

   xmed = fltarr(nbin)
   ymed = fltarr(nbin)
   ysig = fltarr(nbin)

   djs_plot, xdata, ydata, _EXTRA=EXTRA

   ;----------
   ; Fit the median Y in bins of 50 mm, and overplot in green.

   for ibin=0, nbin-1 do begin
      indx = where(xdata GE (xmin + ibin*binsz) $
       AND xdata LT (xmin + (ibin+1)*binsz))
      ymed[ibin] = median(ydata[indx])
      ysig[ibin] = djsig(ydata[indx])
      xmed[ibin] = xmin + (ibin+0.5)*binsz
      djs_oplot, [xmed[ibin]], [ymed[ibin]], color='green', ps=4
      djs_oploterr, xmed[ibin], ymed[ibin], yerr=ysig[ibin], /cap, color='green'
print, xmed[ibin], ymed[ibin], ysig[ibin]
   endfor

   ;----------
   ; Fit a sinusoidal curve to the median Y values, and overplot in red.

   bestchisq = 9e9
   for xoffset=0,350,10 do begin
   for scale=0.1,3,0.02 do begin
      ymodel = sin(scale*(xmed + xoffset) / !radeg)
      res = linfit(ymodel, ymed, chisq=chisq)
      if (chisq LT bestchisq) then begin
         bestchisq = chisq
         bestymodel = ymodel
         bestscale = scale
         bestxoffset = xoffset
         bestres = res
      endif
   endfor
   endfor

   print, 'Y = SIN(',bestscale,' * (X + ', bestxoffset,'))'
   djs_oplot, xmed, bestres[0] + bestres[1]*bestymodel, color='red'

   return
end

;------------------------------------------------------------------------------
pro focal_fit, xfocal, yfocal, xmeasure, ymeasure, weights, xfit, yfit

   dx = 0.03 + 0.02 * (findgen(21) - 10)
   dy = -0.02 + 0.02 * (findgen(21) - 10)
   scale = 0.9996 + 0.0001 * (findgen(11) - 5)
   rot = 0.020 + 0.002 * (findgen(11) - 5) ; In degrees

   bestdist = 1.e9

   for ix=0, n_elements(dx)-1 do begin
   for iy=0, n_elements(dy)-1 do begin
   for is=0, n_elements(scale)-1 do begin
   for ir=0, n_elements(rot)-1 do begin
      cosrot = cos(rot[ir] / !radeg)
      sinrot = sin(rot[ir] / !radeg)
      xmeasure2 = scale[is] * (xmeasure * cosrot + ymeasure * sinrot) + dx[ix]
      ymeasure2 = scale[is] * (ymeasure * cosrot - xmeasure * sinrot) + dy[iy]
      dist2 = total(weights * ( (xmeasure2-xfocal)^2 + (ymeasure2-yfocal)^2 ))
      if (dist2 LT bestdist) then begin
print,dx[ix],dy[iy],scale[is],rot[ir],dist2
         bestdist = dist2
         xfit = xmeasure2
         yfit = ymeasure2
      endif
   endfor
   endfor
   endfor
   endfor

   return
end
;------------------------------------------------------------------------------
pro focal1

   !p.thick = 2
   !x.thick = 2
   !y.thick = 2
   dfpsplot, 'focalplane.ps', /color

   filtsz = 21 ; Filter size for median filtering
   wrange = [4600., 5600.] ; Wavelength range to average flux
   minflux = 200 ; Minimum mean flux to use a fiber

   plugfile = 'plPlugMapM-0174-51775-01.par'

   junk = { filelist, expnum:  0L, dra   :  0.0, ddec  :  0.0 }
   dlist = { filelist, 6124, 0.0, 0.0 }
   dlist = [dlist, { filelist, 6125, -1.0,  0.0 }]
   dlist = [dlist, { filelist, 6128,  1.0,  0.0 }]
   dlist = [dlist, { filelist, 6129,  0.0, -1.0 }]
   dlist = [dlist, { filelist, 6130,  0.0,  1.0 }]
   dlist = [dlist, { filelist, 6131,  0.0,  0.0 }]

;   dlist = { filelist, 6131,  0.0,  0.0 }
   dlist = [dlist, { filelist, 6132, -2.0,  0.0 }]
   dlist = [dlist, { filelist, 6133,  2.0,  0.0 }]
   dlist = [dlist, { filelist, 6134,  0.0, -2.0 }]
   dlist = [dlist, { filelist, 6135,  0.0,  2.0 }]
   dlist = [dlist, { filelist, 6136,  0.0,  0.0 }]

;   dlist = { filelist, 6136,  0.0,  0.0 }
   dlist = [dlist, { filelist, 6137, -3.0,  0.0 }]
   dlist = [dlist, { filelist, 6138,  3.0,  0.0 }]
   dlist = [dlist, { filelist, 6139,  0.0, -3.0 }]
   dlist = [dlist, { filelist, 6140,  0.0,  3.0 }]
   dlist = [dlist, { filelist, 6141,  0.0,  0.0 }]

   ;----------
   ; Construct the data array for the mean fluxes

   nexp = n_elements(dlist)
   nfiber = 640
   fluxarr = fltarr(nfiber, nexp)
   xyoffset = fltarr(nfiber, 2)
   xyofferr = fltarr(nfiber, 2)

   ;----------
   ; Read in the plug map

   yanny_read, plugfile, pp
   plug = *pp[0]
   plug = sortplugmap(plug)

   ;----------
   ; Read the wavelength mapping

   wsetfile1 = 'wset-51779-0174-00006114-b1.fits'
   wset1 = mrdfits(wsetfile1,1)
   xx = 0
   traceset2xy, wset1, xx, loglam1
   wsetfile2 = 'wset-51779-0174-00006114-b2.fits'
   wset2 = mrdfits(wsetfile2,1)
   traceset2xy, wset2, xx, loglam2
   waves = 10.^[[loglam1], [loglam2]]
   wmask = waves GE wrange[0] AND waves LE wrange[1]

   ;----------
   ; Read the fluxes and compute the mean in the specified wavelength range

   for iexp=0, nexp-1 do begin
      fluxfile1 = 'sci-0174-b1-' + string(dlist[iexp].expnum,format='(i8.8)') $
       + '.fits'
      flux1 = mrdfits(fluxfile1)
      ivar1 = mrdfits(fluxfile1, 1)

      fluxfile2 = 'sci-0174-b2-' + string(dlist[iexp].expnum,format='(i8.8)') $
       + '.fits'
      flux2 = mrdfits(fluxfile2)
      ivar2 = mrdfits(fluxfile2, 1)

      flux = [[flux1], [flux2]]
      ivar = [[ivar1], [ivar2]]

      ; Interpolate over bad pixels, then median filter...
      flux = djs_maskinterp(flux, ivar EQ 0, iaxis=0)
      for i=0, nfiber-1 do $
       flux[*,i] = djs_median(flux[*,i], width=filtsz)

      ; Compute the mean flux within the specified wavelength range
      fluxarr[*,iexp] = total(flux * wmask, 1) / total(wmask, 1)

      ; Insist that at least 50% of the pixels are good for a fiber,
      ; or set the flux equal to zero.
      fracgood = total((ivar GT 0) * wmask, 1) / total(wmask ,1)
      fluxarr[*,iexp] = fluxarr[*,iexp] * (fracgood GT 0.5)
   endfor

   ;----------
   ; Loop over each dimension, finding centroids

   for iaxis=0, 1 do begin

      if (iaxis EQ 0) then begin
         iuse = where(dlist.ddec EQ 0)
         dx = dlist[iuse].dra
      endif else begin
         iuse = where(dlist.dra EQ 0)
         dx = dlist[iuse].ddec
      endelse

      for ifiber=0, nfiber-1 do begin

         if (!version.release LT '5.4') then begin
            coeff = poly_fit(dx, transpose(fluxarr[ifiber,iuse]), 2, $
             yfit, junk, junk, corrm)
         endif else begin
            coeff = poly_fit(dx, transpose(fluxarr[ifiber,iuse]), 2, $
             yfit=yfit, covar=corrm, /double)
         endelse

         xyoffset[ifiber,iaxis] = -0.5 * coeff[1] / coeff[2]
         xyofferr[ifiber,iaxis] = - 0.5 * sqrt(corrm[1,1]) / coeff[2] $
          + 0.5 * coeff[1] * sqrt(corrm[2,2]) / (coeff[2])^2

      endfor
   endfor

   ;----------
   ; Select good points to fit

   qgood = bytarr(nfiber) + 1
   qgood = qgood AND (fluxarr[*,0] GT minflux)
;   for iexp=0, nexp-1 do $
;    qgood = qgood AND (fluxarr[*,iexp] GT minflux)
   qgood = qgood AND abs(xyoffset[*,0]) LT 1.5 AND abs(xyoffset[*,1] LT 1.5) $
    AND xyofferr[*,0] LT 0.30 AND xyofferr[*,1] LT 0.30 $
    AND finite(xyoffset[*,0]) AND finite(xyoffset[*,1])
   ibad = where(finite(xyoffset[*,0]) EQ 0 OR finite(xyoffset[*,1]) EQ 0)
   if (ibad[0] NE -1) then begin
      xyoffset[ibad,0] = 0
      xyoffset[ibad,1] = 0
   endif


   ;----------
   ; Allow shifts in X, Y, scale and rotation to the measurements

   xfocal = plug.xfocal
   yfocal = plug.yfocal
   xmeasure = plug.xfocal + xyoffset[*,0]
   ymeasure = plug.yfocal + xyoffset[*,1]
   focal_fit, xfocal, yfocal, xmeasure, ymeasure, qgood, xfit, yfit

   ;----------
   ; Reject the worst points and refit

   dvec = sqrt( (xfit - xfocal)^2 + (yfit - yfocal)^2 )
   djs_iterstat, dvec, sigrej=3.5, mask=outmask
   qgood = qgood * outmask
   focal_fit, xfocal, yfocal, xmeasure, ymeasure, qgood, xfit, yfit

   !p.multi = 0
   igood = where(qgood)
   pscale = 50
   hsize = !d.x_size / 100.
   csize = 1.6
   plot, [0], [0], xrange=[-350,350], yrange=[-350,350], /xstyle, /ystyle, $
    xtitle='X focal [mm]', ytitle='Y focal [mm]', $
    title='Focal Plane Distortions at 5000 Ang', charsize=csize
   arrow, xfocal[igood], yfocal[igood], $
    (xfocal + pscale*(xfit-xfocal))[igood], $
    (yfocal + pscale*(yfit-yfocal))[igood], $
    hsize=hsize, /data
   arrow, -300, -300, -300+pscale, -300, hsize=hsize, /data
   xyouts, -300+pscale+20, -300, '1 arcsec', charsize=csize

print,djsig((xfit-xfocal)[igood]),djsig((yfit-yfocal)[igood])
print,djsig( sqrt( ((xfit-xfocal)[igood])^2 + ((yfit-yfocal)[igood])^2 ) )

wait, 5
   !p.multi = [0,2,2]
   yrange = [-1,1]
   psym = 1
   symsize = 0.5
   focal_binplot, xfocal[igood], (xfit-xfocal)[igood], $
    psym=psym, symsize=symsize, yrange=yrange, $
    xtitle='X focal [mm]', ytitle='\Delta X [arcsec]', $
    title='Focal Plane Mapping at 5000 Ang'
   focal_binplot, xfocal[igood], (yfit-yfocal)[igood], $
    psym=psym, symsize=symsize, yrange=yrange, $
    xtitle='X focal [mm]', ytitle='\Delta Y [arcsec]'
   focal_binplot, yfocal[igood], (xfit-xfocal)[igood], $
    psym=psym, symsize=symsize, yrange=yrange, $
    xtitle='Y focal [mm]', ytitle='\Delta X [arcsec]'
   focal_binplot, yfocal[igood], (yfit-yfocal)[igood], $
    psym=psym, symsize=symsize, yrange=yrange, $
    xtitle='Y focal [mm]', ytitle='\Delta Y [arcsec]'

  dfpsclose

end
