;+
; NAME:
;   smearcorr
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:  
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;   smearcorr is used to calculate and write to file a low order
;   polynomial function which registers the flux in the science exposures
;   to the flux in a fiducial exposure.  
;
; EXAMPLES:
;
; BUGS:
;
;  Blue wavelength region is hardwired: b1 = findgen(60)*4.0e-3 + 3.568
;  Red wavelength region is hardwired : r1 = findgen(54)*4.0e-3 + 3.756
;  Order of polynomial is hardwired:  3
;
; PROCEDURES CALLED:
;   djs_iterstat
;   mrdfits()
;   mwrfits
;   pixelmask_bits()
;   traceset2xy
;   xy2traceset
;   spmedian_rebin
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   17-Oct-2000  Formerly fluxcorr_new -- written by S. Burles
;   09-Oct-2002  Revised by C. Tremonti to calibrate point sources to the 
;                smear and galaxies to the calib images. Added rebin_exposure()
;                and smear_coeff() to streamline. 
;-
;------------------------------------------------------------------------------
function median_rebin, loglam, flux, ivar, color, mask=mask, sigrej=sigrej, $
         outwave = outwave, sn = sn, quality = quality

   if NOT keyword_set(sigrej) then sigrej = 20.0

   ;---------
   ; Define new wavelength sampling

   if color eq 'r' then begin
     w1 = findgen(54)*4.0e-3 + 3.756
     w2 = findgen(54)*4.0e-3 + w1[1]
   endif 
   if color eq 'b' then begin
     w1 = findgen(60)*4.0e-3 + 3.568
     w2 = findgen(60)*4.0e-3 + w1[1]
   endif
   if color eq 'full' then begin
     w1 = findgen(96)*4.0e-3 + 3.568
     w2 = findgen(96)*4.0e-3 + w1[1]
   endif
   range = transpose([[w1],[w2]])
   outwave = djs_median(range,1)

   ;--------------------------
   ; Set up vectors for output
   nr = (size(range))[2]
   ntrace = (size(flux))[2]
   fit = fltarr(nr, ntrace)
   mask = fltarr(nr, ntrace)

   ;--------------------------
   ; Do iterative rejection on each bin of each fiber
   for itrace=0,ntrace-1 do begin
     for irange = 0, nr -1 do begin
        inside = where(loglam[*,itrace] GE range[0,irange] $
                 AND loglam[*,itrace] LT range[1,irange], ninside)
        if ninside GT 0 then begin
           good = where(ivar[inside,itrace] GT 0, ngood)
       
           if ngood GT 1 then begin 
              djs_iterstat, flux[inside[good],itrace], median=md, sigma=sig
              fit[irange, itrace] = md
              if ngood GT 0.5 * ninside then $
                   mask[irange, itrace]  = total(ivar[inside[good],itrace])
              sn = sig * sqrt(mask[irange, itrace])
              if sn GT sigrej OR sn LE 0 then mask[irange, itrace] = 0.0
           endif
        endif
     endfor
   endfor

   sn = djs_median(flux * sqrt(ivar),1)

   badval = ( fibermask_bits('NOPLUG')   OR fibermask_bits('BADTRACE') $
       OR fibermask_bits('BADFLAT')  OR fibermask_bits('BADARC')   $
       OR fibermask_bits('NEARWHOPPER') )

   quality = (mask[0,*] AND badval)

   return, fit
end

;------------------------------------------------------------------------------
pro frame_flux_tweak, bloglam, rloglam, bflux, rflux, bivar, rivar, $
    best_exp, plugtag, corrfile

   expid = plugtag[uniq(plugtag.expid)].expid
   splog, 'Best Exposure =', best_exp
   indx = where(plugtag.expid eq best_exp)

   ;----------
   ; Create red & blue smoothed best images

   bfit = spmedian_rebin(bloglam[*,indx], bflux[*,indx], bivar[*,indx], $
          'b', outwave = bwave, mask = bmask, sn = bsn, quality = bquality)

   rfit = spmedian_rebin(rloglam[*,indx], rflux[*,indx], rivar[*,indx], $
          'r', outwave = rwave, mask = rmask, sn = rsn, quality = rquality)

   ;---------------------------------------------------------
   ; Now put blue and red together
    
   bestflux = [bfit,rfit]
   bestivar = [bmask, rmask]
   bestsnmed = transpose([[bsn],[rsn]])
   qbest =  (bquality OR rquality)

   npix = (size(bestflux,/dimens))[0]
   nfiber = (size(bestflux,/dimens))[1]
   wave = [bwave,rwave] # replicate(1,nfiber)

   ; --------------------------------------------------------
   ;  Create basic mask
    
   bestmask = lonarr(nfiber)

   ;----------
   ; Loop through the science images

   nfiles = n_elements(corrfile)

   for ifile=0, nfiles-1 do begin

      indx = where(plugtag.expid eq expid[ifile])
      splog, 'Fluxing science image =', expid[ifile]

      ;------------
      ; Create coeffients of fit and set them to a zero-order polynomaial
      ; with an amplitude of 1.0
 
      fitimg = bestflux*0.0 + 1.0
      xy2traceset, wave, fitimg, corrset, ncoeff=4, /silent
      corrset.coeff[0,*] = 1.0
      corrset.coeff[1:*,*] = 0.0
      thismask  = bestmask 

      ;--------------
      ; If the science image and best image are the same,
      ; force their ratio to be unity.
      if (expid[ifile] EQ best_exp) then begin
        mwrfits, corrset, corrfile[ifile], /create
      ; mwrfits, thismask, corrfile[ifile]
        continue
      endif 

      ;-------------
      ; Create rebinned blue and red science images

; Need to figure out the best way to set mask bits!!!

      bscifit = spmedian_rebin(bloglam[*,indx], bflux[*,indx], bivar[*,indx], $
                'b', mask = bscimask, sn = bscisn, quality = bsciquality)

      rscifit = spmedian_rebin(rloglam[*,indx], rflux[*,indx], rivar[*,indx], $
                'r', mask = rscimask, sn = rscisn, quality = rsciquality)
       
      ;----------------
      ; Combine red and blue

      sciflux = [bscifit,rscifit]
      sciivar = [bscimask,rscimask]
      scisnmed = transpose([[bscisn],[rscisn]])
      qsci = lonarr(nfiber)
      qsci = qsci OR bsciquality
      qsci = qsci OR rsciquality

      ;-------------------------------------------------------------
      ; Determine the S/N of the expsure and act accrodingly
      ; poly 3:  High S/N => 4 order legendre fit (in loglam)
      ; poly 2:  Medium S/N => 3 order legendre fit
      ; poly 1:  Low S/N => find average of nearest high S/N 

      fiber_coeff = lonarr(nfiber) + 1

      poly3 = where(scisnmed[0,*] GT 2.5 AND scisnmed[1,*] GT 5.0  $
               AND  bestsnmed[0,*] GT 2.5 AND  bestsnmed[1,*] GT 5.0 $
               AND  qsci EQ 0 AND qbest EQ 0, npoly3) 
      if poly3[0] NE -1 then fiber_coeff[poly3] = 3 

      poly2 = where(scisnmed[0,*] GT 1.0 AND scisnmed[1,*] GT 2.0  $
               AND  bestsnmed[0,*] GT 1.0 AND bestsnmed[1,*] GT 2.0 $
               AND  qsci EQ 0 AND qbest EQ 0 AND fiber_coeff NE 3, npoly2) 
      if poly2[0] NE -1 then fiber_coeff[poly2] = 2

      poly1 = where(fiber_coeff eq 1 and $
                    scisnmed[0,*] ne 0 and scisnmed[1,*] ne 0 and $
                    bestsnmed[0,*] ne 0 and bestsnmed[1,*] ne 0 and $
                    strmatch(plugtag[indx].objtype, '*SKY*') ne 1, npoly1)

      ;-----------------
      ; Calculate polynomial coeffs for the 3 cases
   
      if npoly3 gt 0 then begin
        xy2traceset, wave[*,poly3], bestflux[*,poly3], polyset, $
            invvar=bestivar[*,poly3], ncoeff=4, inputfunc=sciflux[*,poly3], $
            lower = 3, upper = 3
        corrset.coeff[*,poly3] = polyset.coeff
      endif

      if npoly2 gt 0 then begin
        xy2traceset, wave[*,poly2], bestflux[*,poly2], polyset, $
            invvar=bestivar[*,poly2], ncoeff=3, inputfunc=sciflux[*,poly2], $
            lower = 3, upper = 3
        corrset.coeff[0,poly2] = polyset.coeff[0,*]
        corrset.coeff[1,poly2] = polyset.coeff[1,*]
        corrset.coeff[2,poly2] = polyset.coeff[2,*]
      endif

      ;---------------------
      ; For the lowest S/N case use the median of the nearest 5 high or 
      ; moderate S/N fibers, then adjust the zeropoint

      if npoly1 gt 0 and (npoly2 gt 0 or npoly3 gt 0) then begin
        if npoly3 gt 0 then hisn = poly3 
        if npoly2 gt 0 then hisn = poly2
        if npoly3 gt 0 and npoly2 gt 0 then hisn = [poly3, poly2] 

        for ifib = 0, npoly1 - 1 do begin
          dist = (plugtag[poly1[ifib]].xfocal - plugtag[hisn].xfocal)^2 + $
                 (plugtag[poly1[ifib]].yfocal - plugtag[hisn].yfocal)^2
          nearindx = sort(dist)
          near5 = hisn[nearindx[0:4]]
          corrset.coeff[*,poly1[ifib]] = djs_median(corrset.coeff[*,near5], 2)
        endfor
           
        traceset2xy, corrset, wave, corrtemp
        xy2traceset, wave[*,poly1], bestflux[*,poly1], polyset, $
          invvar=bestivar[*,poly1], ncoeff=1, $
          inputfunc=sciflux[*,poly1] * corrtemp[*,poly1], lower = 3, upper = 3

        zptcor = rebin(polyset.coeff, 4, npoly1)
        corrset.coeff[*,poly1] = corrset.coeff[*,poly1] * zptcor
      endif

      ;----------
      ; Identify sky fibers from the plug map   
      sky = where(strmatch(plugtag[indx].objtype, '*SKY*'))
      corrset.coeff[0,sky] = 1.0
      corrset.coeff[1:*,sky] = 0.0

      ;---------------
      ; Reject any bad corrections and replace with low S/N solution
      ; solution.  Bad vectors are those with fit coefficients outside
      ; of pre-set boundaries.  The boundaries on coeffs 1 & 2 increase 
      ; with increasing spread in the coeff0 values -- this should help
      ; to accomodate really bad plates.

      ; Normalize coefficients 
      coef0 = corrset.coeff[0,*]  
      coef1 = corrset.coeff[1,*] / coef0 
      coef2 = corrset.coeff[2,*] / coef0 
      coef3 = corrset.coeff[3,*] / coef0 
 
      meanclip, coef0, coef0mean, coef0sig, clipsig = 5
      coef0 = coef0 / coef0mean
      coef0sig = (coef0sig / coef0mean) > 0.10
    
      ; Adjust boundary values for plate quality
      coef1bound = 3.5 * coef0sig  
      coef2bound = 2.5 * coef0sig 
      coef3bound = 2.5 * coef0sig 

      isbad = coef0 LT  (1 - 5 * coef0sig) OR coef0 GT (1 + 5 * coef0sig) OR $
              coef1 LT -coef1bound OR coef1 GT coef1bound OR $
              coef2 LT -coef2bound OR coef2 GT coef2bound OR $
              coef2 LT -coef3bound OR coef2 GT coef3bound

      bad = where(isbad, nbad)
      if bad[0] NE -1 then begin
        splog, 'Warning: Large deviations in flux correction '
        splog, 'Warning: Replacing with zero point shift:', $
        string(bad + 1)

      !P.MULTI = [0, 1, 2]
      for ii = 0, nbad - 1 do begin
        ifib = bad[ii]
        plot, 10.0^wave[*,ifib], corrtemp[*,ifib], yr=[0.0, 2.0], /nodata, $
                ytitle='Best Frame Flux / Science Frame Flux', title='Bad Vector'
        oplot, 10.0^wave[*,ifib], bestflux[*,ifib]/sciflux[*,ifib], psym=6, $
               syms=0.5, thick=3
        djs_oplot, 10.0^wave[*,ifib], corrtemp[*,ifib], color='red', thick=3
       endfor      
 
        xy2traceset, wave[*,bad], bestflux[*,bad], polyset, $
            invvar=bestivar[*,bad], ncoeff=1, inputfunc=sciflux[*,bad], $
            lower = 3, upper = 3
        corrset.coeff[0,bad] = polyset.coeff
        corrset.coeff[1:*,bad] = 0.0
      endif

      ;---------------
      ; Set maskbits and append to end of corrfile
         
      ;thismask = thismask OR (fibersn EQ 2) * pixelmask_bits('SMEARMEDSN')
      ;thismask = thismask OR (fibersn EQ 3) * pixelmask_bits('SMEARHIGHSN')

      ;------------
      ; Write out as FITS

      mwrfits, corrset, corrfile[ifile], /create
      ;mwrfits, thismask, corrfile[ifile]

      ;------------
      ; Plot correction vectors

      if keyword_set(noplot) then continue

      traceset2xy, corrset, wave, corrimage

      ;---------------------
      ; Show individual fits to standards and examples of 3 S/N bins

      !P.MULTI = [0, 1, 2]
      std = where(strmatch(plugtag[indx].objtype, '*_STD*'), nstd)
      for ii = 0, nstd - 1 do begin
        istd = std[ii]
        plot, 10.0^wave[*,istd], corrimage[*,istd], yr=[0.5, 1.5], /nodata, $
                ytitle='Best Frame Flux / Science Frame Flux', $
                title='Standard Star'
        oplot, 10.0^wave[*,istd], bestflux[*,istd]/sciflux[*,istd], psym=6, $
               syms=0.5, thick=3
        djs_oplot, 10.0^wave[*,istd], corrimage[*,istd], color='red', thick=3
      endfor 
      
      for ii = 0, (10 < npoly2) - 1 do begin
        istd = poly2[ii]
        plot, 10.0^wave[*,istd], corrimage[*,istd], yr=[0.5, 1.5], /nodata, $
                ytitle='Best Frame Flux / Science Frame Flux', title='Poly 2'
        oplot, 10.0^wave[*,istd], bestflux[*,istd]/sciflux[*,istd], psym=6, $
               syms=0.5, thick=3
        djs_oplot, 10.0^wave[*,istd], corrimage[*,istd], color='red', thick=3
      endfor 

      for ii = 0, (10 < npoly1) - 1 do begin
        istd = poly1[ii]
        plot, 10.0^wave[*,istd], corrimage[*,istd], yr=[0.0, 2.0], /nodata, $
                ytitle='Best Frame Flux / Science Frame Flux', title='Poly 1'
        oplot, 10.0^wave[*,istd], bestflux[*,istd]/sciflux[*,istd], psym=6, $
               syms=0.5, thick=3
        djs_oplot, 10.0^wave[*,istd], corrimage[*,istd], color='red', thick=3
      endfor 

      ;-----------------
      ; Plot all low and medium S/N correction vectors

      !P.MULTI = [0, 1, 2]

      djs_plot, 10.0^wave, corrimage, /nodata, yr=[0.4, 1.6], $
                xr=[min(10.0^wave)-100,max(10.0^wave)+100], $
                /xstyle, /ystyle, xtitle='\lambda [\AA]', $
                ytitle='Best Frame Flux / Science Frame Flux'

      for iobj=0, npoly2 -1 do $
        djs_oplot, 10.0^wave[*,poly2[iobj]], corrimage[*,poly2[iobj]], $
                   color='blue', nsum=5
      for iobj=0, npoly1 -1 do $
        djs_oplot, 10.0^wave[*,poly1[iobj]], corrimage[*,poly1[iobj]], $
                   color='magenta', nsum=5

      legend, ['High S/N', 'Med S/N', 'Low S/N', 'Standard Star'], psym=0, $
              thick = 3, color=djs_icolor(['green', 'blue', 'magenta', 'red'])

      ;-----------------
      ; Plot all high S/N correction vectors + standards

      djs_plot, 10.0^wave, corrimage, /nodata, yr=[0.4, 1.6], $
                xr=[min(10.0^wave)-100,max(10.0^wave)+100], $
                /xstyle, /ystyle, xtitle='\lambda [\AA]', $
                ytitle='Best Frame Flux / Science Frame Flux'

      for iobj=0, npoly3 -1 do $
        djs_oplot, 10.0^wave[*,poly3[iobj]], corrimage[*,poly3[iobj]], $
                   color='green', nsum=5

      std = where(strmatch(plugtag[indx].objtype, '*_STD*'), nstd)
      for iobj=0, nstd -1 do $
        djs_oplot, 10.0^wave[*,std[iobj]], corrimage[*,std[iobj]], $
                   color='red', nsum=5, thick=2

      !P.MULTI = 0
   endfor
   return
end
