;+
; NAME:
;   lrgmodel_plots
;
; PURPOSE:
;   Generate goodness-of-fit plots from the LRG model photo-z fits.
;
; CALLING SEQUENCE:
;   lrgmodel_plots, [ suffix, public=, /recalibrate, /colorcut, $
;    subsamp=, _EXTRA= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   suffix         - Suffix name for plot names, which is 'photoz'+SUFFIX+'.ps'
;   public         - Select public data; this keyword is passed
;                    to LRGMODEL_READ_SPALL()
;   recalibrate    - Recalibrate the photometry using SDSS_RECALIBRATE
;   colorcut       - If set, then make Nikhil's color-cuts on the 2dF data
;   subsamp        - Subsample the spAll galaxies by this factor; default to 4
;   _EXTRA         - Additional keywords to pass to LRGMODEL_PHOTOZ(),
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine reads both the SDSS spectra from spAll, and the SDSS-2dF
;   redshifts.
;
;   Some of the plots are replicated with the SDSS-2dF points in color.
;   The SDSS-2dF 2003A objects are in green, and the 2003B objects are in blue.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_oplot
;   djs_plot
;   hogg_scatterplot
;   lrgmodel_append_twodf()
;   lrgmodel_read_spall()
;   plothist
;   sdss_recalibrate
;
; REVISION HISTORY:
;   18-Dec-2003  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro lrgmodel_plots, suffix, public=public, recalibrate=recalibrate, $
 colorcut=colorcut, subsamp=subsamp, _EXTRA=KeywordsForPhotoz

   if (NOT keyword_set(suffix)) then suffix = ''
   if (NOT keyword_set(subsamp)) then subsamp = 4
   csize = 2.5
   thick = 2

   ;----------
   ; Read the data files

   spall = lrgmodel_read_spall(public=public)
   if (subsamp GT 1) then begin
      indx = lindgen(floor(float(n_elements(spall))/subsamp)) * subsamp
      spall = spall[indx]
   endif

;spall = spall[0:10000] ; Test ???
   nsdss = n_elements(spall)
   spall = lrgmodel_append_twodf(spall,'2003A',colorcut=colorcut, nadd=nadd1)
   spall = lrgmodel_append_twodf(spall,'2003B',colorcut=colorcut, nadd=nadd2)

   ;----------
   ; Optionally recalibrate the photometry

   if (keyword_set(recalibrate)) then begin
      splog,'Recalibrating...'
      sdss_recalibrate,spall
      splog,'...done'
   endif

   ;----------
   ; Trim to Nikhil's LRG target selection

   lrgbits = sdf_select_target(spall, /all)
   qtrim = lrgbits NE 0
   itrim = where(qtrim, ntrim)
   spall = spall[itrim]
   if (nsdss GT 0) then nsdss = total(qtrim[0:nsdss-1])
   if (nadd1 GT 0) then nadd1 = total(qtrim[nsdss:nsdss+nadd1-1])
   if (nadd2 GT 0) then nadd2 = total(qtrim[nsdss+nadd1:nsdss+nadd1+nadd2-1])

   ;----------
   ; Compute the colors

   grcolor = - 2.5*alog10((spall.modelflux[1]>0.1) / (spall.modelflux[2]>0.1))
   ricolor = - 2.5*alog10((spall.modelflux[2]>0.1) / (spall.modelflux[3]>0.1))

   ;----------
   ; Fit the photo-z's

   splog,'Fitting photo-zs'

   zfit = lrgmodel_photoz(spall.modelflux, spall.modelflux_ivar, $
    z_err=zfit_err, extinction=spall.extinction, /abcorrect, $
    _EXTRA=KeywordsForPhotoz)

   zdiff = (zfit - spall.z)
   azdiff = abs(zdiff)

   ;----------
   ; Print the standard deviations in redshift bins

   print, " zMin zMax   Nobj  <z-zFit>    68%cf    95%cf    99%cf"
   for i=0, 17 do begin
      if (i LT 17) then begin
         zmin = 0.05 * i
         zmax = 0.05 * (i+1)
      endif else begin
         zmin = 0.10
         zmax = 0.80
      endelse
      indx = where(spall.z GT zmin AND spall.z LE zmax, ct)
      if (ct GT 0) then begin
         isort = indx[sort(azdiff[indx])]
         i68 = isort[0.68*ct]
         i95 = isort[0.95*ct]
         i99 = isort[0.99*ct]
         print, zmin, zmax, ct, median([zdiff[indx]]), $
          azdiff[i68], azdiff[i95], azdiff[i99], $
          format='(2f5.2,i7,f10.4,3f9.4)'
      endif
   endfor

   ;----------
   ; PLOT: Photo-z vs. spectro-z (points)

   dfpsplot, 'photoz'+suffix+'.ps', /square, /color

   djs_plot, spall[0:nsdss-1].z, zfit[0:nsdss-1], $
    psym=3, xr=[0,0.8], yr=[0,0.8], $
    xtitle='Spectroscopic Z', ytitle='Photometric Z', charsize=csize, $
    thick=thick, charthick=thick
   djs_oplot, !x.crange, !x.crange, color='red', thick=thick

   ; Re-plot with the SDSS-2dF points in color
   djs_plot, spall[0:nsdss-1].z, zfit[0:nsdss-1], $
    psym=3, xr=[0,0.8], yr=[0,0.8], $
    xtitle='Spectroscopic Z', ytitle='Photometric Z', charsize=csize, $
    thick=thick, charthick=thick
   indx1 = nsdss + lindgen(nadd1)
   indx2 = nsdss + nadd1 + lindgen(nadd2)
   djs_oplot, spall[indx1].z, zfit[indx1], psym=1, symsize=0.25, $
    color='green', thick=thick
   djs_oplot, spall[indx2].z, zfit[indx2], psym=1, symsize=0.25, $
    color='blue', thick=thick
   djs_oplot, !x.crange, !x.crange, color='red', thick=thick

   ;----------
   ; PLOT: Photo-z vs. spectro-z (greyscale)

   ii = where(zfit GT 0.01)
   hogg_scatterplot, spall[ii].z, zfit[ii], xr=[0,0.8], yr=[0,0.8], $
    xtitle='Spectroscopic Z', ytitle='Photometric Z', /conditional, $
    quantile=[0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975], charsize=csize
   djs_oplot, !x.crange, !x.crange, color='red'
   djs_oplot, [0.1,0.1], !y.crange, color='red'

   ;----------
   ; PLOT: Delta-z vs. spectro-z (greyscale)

   ii = where(zfit GT 0.01)
   hogg_scatterplot, spall[ii].z, zdiff[ii], xr=[0,0.8], yr=[-0.3,0.3], $
    xtitle='Spectroscopic Z', ytitle='z(photo) - z(spectro)', /conditional, $
    quantile=[0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975], charsize=csize
   djs_oplot, !x.crange, [0,0], color='red'
   djs_oplot, [0.1,0.1], !y.crange, color='red'

   ;----------
   ; PLOT: Delta-z vs. photo-z (greyscale)

   ii = where(zfit GT 0.01)
   hogg_scatterplot, zfit[ii], zdiff[ii], xr=[0,0.8], yr=[-0.3,0.3], $
    xtitle='Photometric Z', ytitle='z(photo) - z(spectro)', /conditional, $
    quantile=[0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975], charsize=csize
   djs_oplot, !x.crange, [0,0], color='red'
   djs_oplot, [0.1,0.1], !y.crange, color='red'

   ;----------
   ; PLOT: Delta-z histogram

   jj = where(spall.z GT 0.10)
   binsz = 0.005
   xrange = [-0.3,0.3]
   plothist, zdiff[jj], xrange=xrange, bin=binsz, /xstyle, $
    xtitle='z(photo) - z(spectro)', charsize=csize
   xvec = findgen((xrange[1]-xrange[0])/binsz) * binsz + xrange[0]
   yvec1 = 6000*exp(-(xvec/0.023)^2)
   yvec2 = 700*exp(-(xvec/0.023/3)^2)
   djs_oplot, xvec, yvec1, color='red'
   djs_oplot, xvec, yvec2, color='red'
   djs_oplot, xvec, yvec1+yvec2, color='green'

   ;----------
   ; PLOT: Delta-z vs. magnitude (greyscale)

   rmag = 22.5 - 2.5*alog10(spall.modelflux[2] > 1)
   ii = where(zfit GT 0.01)
   hogg_scatterplot, rmag[ii], zdiff[ii], xr=[14,23], yr=[-0.3,0.3], $
    xtitle='r-mag', ytitle='z(photo) - z(spectro)', /conditional, $
    quantile=[0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975], charsize=csize
   djs_oplot, !x.crange, [0,0], color='red'

   ; Re-plot with the SDSS-2dF points in color
   hogg_scatterplot, rmag[ii], zdiff[ii], xr=[14,23], yr=[-0.3,0.3], $
    xtitle='r-mag', ytitle='z(photo) - z(spectro)', /conditional, $
    quantile=[0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975], charsize=csize
   djs_oplot, !x.crange, [0,0], color='red'
   djs_oplot, rmag[indx1], zdiff[indx1], psym=1, symsize=0.25, color='green'
   djs_oplot, rmag[indx2], zdiff[indx2], psym=1, symsize=0.25, color='blue'

   ;----------
   ; PLOT: Delta-z vs. (g-r) (greyscale)

   ii = where(zfit GT 0.01)
   hogg_scatterplot, grcolor[ii], zdiff[ii], xr=[0.5,3.0], yr=[-0.3,0.3], $
    xtitle='(g-r)', ytitle='z(photo) - z(spectro)', /conditional, $
    quantile=[0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975], charsize=csize
   djs_oplot, !x.crange, [0,0], color='red'

   ; Re-plot with the SDSS-2dF points in color
   hogg_scatterplot, grcolor[ii], zdiff[ii], xr=[0.5,3.0], yr=[-0.3,0.3], $
    xtitle='(g-r)', ytitle='z(photo) - z(spectro)', /conditional, $
    quantile=[0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975], charsize=csize
   djs_oplot, !x.crange, [0,0], color='red'
   djs_oplot, grcolor[indx1], zdiff[indx1], psym=1, symsize=0.25, color='green'
   djs_oplot, grcolor[indx2], zdiff[indx2], psym=1, symsize=0.25, color='blue'

   ;----------
   ; PLOT: Delta-z vs. [(r-i)-(g-r)/8.] (greyscale)

   xaxis = ricolor - grcolor/8.
   ii = where(zfit GT 0.01)
   hogg_scatterplot, xaxis[ii], zdiff[ii], xr=[-0.5,2.0], yr=[-0.3,0.3], $
    xtitle='(r-i)-(g-r)/8', ytitle='z(photo) - z(spectro)', /conditional, $
    quantile=[0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975], charsize=csize
   djs_oplot, !x.crange, [0,0], color='red'

   ; Re-plot with the SDSS-2dF points in color
   hogg_scatterplot, xaxis[ii], zdiff[ii], xr=[-0.5,2.0], yr=[-0.3,0.3], $
    xtitle='(r-i)-(g-r)/8', ytitle='z(photo) - z(spectro)', /conditional, $
    quantile=[0.025, 0.05, 0.16, 0.5, 0.84, 0.95, 0.975], charsize=csize
   djs_oplot, !x.crange, [0,0], color='red'
   djs_oplot, xaxis[indx1], zdiff[indx1], psym=1, symsize=0.25, color='green'
   djs_oplot, xaxis[indx2], zdiff[indx2], psym=1, symsize=0.25, color='blue'

   dfpsclose

   return
end
;------------------------------------------------------------------------------
