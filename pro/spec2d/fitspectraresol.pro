;+
; NAME:
;   fitspectraresol
;
; PURPOSE:
;   Fit polynomials to the emission line width for each arc fiber bundle
;
; CALLING SEQUENCE:
;   reslset = fitspectraresol(arc_flux, arc_fluxivar, xcen, wset $
;    [ ncoeff=, xmin=, xmax=, medresol=, numbundles= ] )
;
; INPUTS:
;   arc_flux     - arc image extracted flux
;   arc_fluxivar - corresponding inverse variance
;   xcen         - xpeaks of arc lines [ntrace,nlines]
;   wset         - Trace set for initial wavelength solution in row number ROW.
;
; OPTIONAL KEYWORDS:
;   ncoeff     - Order of legendre polynomial to apply to width vs. row;
;                default to 4.
;   xmin       - Lowest row number for trace set; default to 0.
;   xmax       - Highest row number for trace set; default to 2047.
;   numbundles - The number of fiber bundles
;   quick      - Flag used during quick reduction on the mountain in
;                apo routines using only middle region of ccd
;
; OUTPUTS:
;   reslset   - Traceset structure containing fit coefficients
;
; OPTIONAL OUTPUTS:
;   medresol  - Median resolution in each of the 4 quadrants
;               of the CCD, ordered LL,LR,UL,UR.
;
; COMMENTS:
;   Used to fill arcstruct.dispset, which can then be applied
;   to PSF-corrected sky subtraction.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   splog
;   xy2traceset
;
; REVISION HISTORY:
;   Basede on fitdispersion.pro
;   25-Nov-2020  Written by Hector Ibarra-Medel, UIUC
;-
;------------------------------------------------------------------------------
function fitspectraresol, arc_flux, arc_fluxivar, xcen_inp, wset, $
                        ncoeff=ncoeff, xmin=xmin, xmax=xmax, medresol=medresol, $
                        numBundles = numBundles, quick=quick, resol_final=resol_final, $
                        bundlefibers=bundlefibers, arc_test_file=arc_test_file, waves=waves

   if (NOT keyword_set(ncoeff)) then ncoeff = 4
   if (NOT keyword_set(xmin)) then xmin = 0.0
   if (NOT keyword_set(xmax)) then xmax = 2047.0

   nline = (size(xcen_inp,/dimen))[1]
   ntrace = (size(xcen_inp,/dimen))[0]
   npix = (size(arc_flux,/dimen))[0]
   
   if (nline LT 3) then begin
      splog, 'WARNING: Too few lines for spectral resolution traceset (', nline, ')'
      return, -1
   endif

   ;----------
   ; Sort XCEN, which is necessary for the call to EXTRACT_IMAGE.

   isort = sort(xcen_inp[0,*])
   xcen = xcen_inp[*,isort]

   ;----------
   ; Construct a mask that is nonzero only for those pixels within
   ; +/- 12 pixels from each arc line center on each fiber.

   arcmask = make_array(size=size(arc_flux), /byte)
   for itrace=0, ntrace-1 do begin
      for iline=0, nline-1 do begin
         i1 = floor(xcen[itrace,iline] - 12) > 0
         i2 = ceil(xcen[itrace,iline] + 12) < (npix - 1)
         if (i1 LT npix-1 AND i2 GT 0) then arcmask[i1:i2,itrace] = 1
      endfor
   endfor

   ;----------
   ; Extract arc lines from the [NPIX,NTRACE] image, measuring the
   ; wavelength sigmas during extraction for every fiber.

   ;extract_image, arc_flux, arc_fluxivar*arcmask, xcen, sigma, $
   ; arclineflux, arclineivar, ansimage=ansimage, wfixed=[1,1], $
   ; highrej=10, lowrej=10, relative=1, npoly=5
   
   traceset2xy,wset,dum,loglam
   lamb=10^loglam
   dpix=5; pixel window to perform the gaussian fit
   spc_r= fltarr(nline,ntrace)
   if keyword_set(arc_test_file) then begin
        ;mydevice = !D.NAME
        ;print, REPSTR(arc_test_file,'.fits','.ps')
        ;dfpsplot,REPSTR(arc_test_file,'.fits','_499.ps'),/square,/color
        ;loadct,39
        subtitle = 'fitspectraresol: '+strjoin((strsplit(arc_test_file,'-',/extract))[1:2],' ')+'fiber '
   endif
   splog,'Calculating the spectra resolution'
   for il=0, ntrace-1 do begin
     indx=transpose(uint(xcen[il,*]))
     for it=0, nline-1 do begin
       x_in=lamb[indx[it]-dpix:indx[it]+dpix,il]
       y_in=arc_flux[indx[it]-dpix:indx[it]+dpix,il]
       yfit = mpfitpeak(x_in, y_in, a, error=sy, nterms=3)
       ;Rs=lamb[indx[it]]/(a[2]*2.0*sqrt(2.0*alog(2.0)))
       Rs=a[2]*2.0*sqrt(2.0*alog(2.0))
       spc_r[it,il]=Rs
       if il eq 499 then begin
         splog,'Line '+strtrim(string(it))+' '+string(Rs)+' '+string(lamb[indx[it]])
         if keyword_set(arc_test_file) then begin
            plot,x_in,y_in, title = subtitle+strtrim(il+1,2)
            oplot,x_in,yfit,color=djs_icolor('red')
         endif
       endif
     endfor
   endfor
;   if keyword_set(arc_test_file) then begin
;        device, /close
;        SET_PLOT, mydevice
;   endif

   if keyword_set(arc_test_file) then begin
        foreach ifib, [1,100,200,250,300,400,499] do begin
            ;mydevice = !D.NAME
            ;splog, REPSTR(arc_test_file,'.fits','.ps')
            ;dfpsplot,REPSTR(arc_test_file,'.fits','_'+strtrim(ifib,2)+'.ps'),/square,/color
            ;loadct,39
            for il=0, ntrace-1 do begin
                indx=transpose(uint(xcen[il,*]))
                for it=0, nline-1 do begin
                    x_in=lamb[indx[it]-dpix:indx[it]+dpix,il]
                    y_in=arc_flux[indx[it]-dpix:indx[it]+dpix,il]
                    yfit = mpfitpeak(x_in, y_in, a, error=sy, nterms=3)
                    if il eq ifib then begin
                        plot,x_in,y_in, title = subtitle+strtrim(il+1,2)
                        oplot,x_in,yfit,color=djs_icolor('red')
                    endif
                endfor
            endfor
            ;device, /close
            ;SET_PLOT, mydevice
        endforeach
   endif

   ;----------
   ; Mask values
   gmask = (spc_r GT 0) * (spc_r LT 4.5)
   igood = where(gmask)
   resol = make_array(size=size(spc_r), /float)
   if (igood[0] NE -1) then $
     resol[igood] = spc_r[igood]
   
   ;----------
   ; Perform median across bundles on good arclines only
   ; somewhat tedious, but it works

;   resol = reform(resol,nline,20,numbundles)
;   gmask = reform(gmask,nline,20,numbundles)
   resol_bundle = fltarr(nline,numbundles)

   bundlenum = intarr(total(bundlefibers))
   resol_final = fltarr(nline, ntrace)
   
   start = 0
   for i= 0, numbundles-1 do begin
       endf = start + bundlefibers[i]
       bundlenum[start: endf-1] = i
       start = endf
   endfor

   t_lo = intarr(numbundles)
   for ibun=1, numbundles-1 do t_lo[ibun]=total(bundlefibers[0:ibun-1])
   t_hi = t_lo + bundlefibers -1

   for iline=0, nline-1 do begin
     for j=0, numbundles-1 do begin
        bunfib = where(bundlenum eq j)
        ss = where(gmask[iline,bunfib] AND resol[iline,bunfib] GT 0, ct)
        if (ct GE 0.5*bundlefibers[j]) then $
         resol_bundle[iline,j] = djs_median(resol[iline,bunfib[ss]])
        resol_final[iline,[t_lo[j]:t_hi[j]]] = resol_bundle[iline,j]
     endfor
   endfor

;   resol_final = congrid(resol_bundle, nline, ntrace)
   ;----------
   ; Turn the widths back into a traceset.

   ; ASB: adding maxdev=0.2
   ;ncoeff=4
   xy2traceset, transpose(xcen), resol_final, reslset, inmask=(resol_final GT 0), $
    ncoeff=ncoeff, xmin=xmin, xmax=xmax, maxdev=0.2

   if keyword_set(arc_test_file) then begin
;    debug_fits = 'fitspectraresol.fits'
     mwrfits_named, resol,         arc_test_file, name='RESOL', /create
     mwrfits_named, gmask,         arc_test_file, name='GMASK'
     mwrfits_named, resol_bundle,  arc_test_file, name='RESOL_BUNDLE'
     mwrfits_named, resol_final,   arc_test_file, name='RESOL_FINAL'
     mwrfits_named, reslset,       arc_test_file, name='RESLSET'
     mwrfits_named, waves,         arc_test_file, name='WAVES'
     mwrfits_named, spc_r,         arc_test_file, name='SPC_R'
   endif

   ;----------
   ; Compute the widths in each of 4 quandrants on the CCD
   ; as the median of the unmasked pixels around the lines being fit

   ; matt modify for 4x4 grid, use quadrupole terms and only half ccd
   ; for quick extract 
   ; otherwise do old quadrant fit
   traceset2xy, reslset, xx, resol_fit
 
   ;plot,lamb[transpose(uint(xcen[250,*]))],resol_final[*,250]
   ;oplot,lamb[*,250],resol_fit[*,250],color=djs_icolor('red')
   ;dfpsclose
 
   if keyword_set(quick) then begin
       x1 = [npix/4,3*npix/8,3*npix/8,5*npix/8]
       x2 = [3*npix/8-1,5*npix/8-1,5*npix/8-1,3*npix/4-1]
       y1 = [3*ntrace/8,ntrace/4,5*ntrace/8,3*ntrace/8]
       y2 = [5*ntrace/8-1,3*ntrace/8-1,3*ntrace/4-1,5*ntrace/8-1]
   endif else begin
       x1 = [0,0,npix/2,npix/2]
       x2 = [npix/2-1,npix/2-1,npix-1,npix-1]
       y1 = [0,ntrace/2,0,ntrace/2]
       y2 = [ntrace/2-1,ntrace-1,ntrace/2-1,ntrace-1]
   endelse

   medresol = fltarr(4)
   for i=0,3 do begin
      indx = where(arcmask[x1[i]:x2[i],y1[i]:y2[i]],ct)
      if (ct GT 0) then $
       medresol[i] = $
        median([ (resol_fit[x1[i]:x2[i],y1[i]:y2[i]])[indx] ])
   endfor
   
   splog, 'Median wavelength R = ' $
    + string(medresol,format='(4f5.2)') + ' pix (LL LR UL UR)' ;left bottom top right

   return, reslset
end
;------------------------------------------------------------------------------
