;+
; NAME:
;   frame_flux_tweak
;
; PURPOSE:
;   To genreate low order correctons to each red-blue exposure pair which 
;   match the flux (on a fiber-by-fiber basis) to that in a fiducial exposure
;   selected to have the best S/N and spectrophotometry. This is necessary 
;   so that multiple exposures can be combined with effective cosmic ray 
;   rejection.  The flux vectors of each exposure should be multiplied by the 
;   resultant 4th order legnedre polynomial.
;
; CALLING SEQUENCE:
;   frame_flux_tweak, bloglam, rloglam, bflux, rflux, bivar, rivar, $
;     best_exp, plugtag, corrfile, diag=, title=
;
; INPUTS:
;   bloglam  - blue wavelength array in (log10 Ang.) [nblue_pix, nfiber*nexp]
;   rloglam  - red wavelength array in (log10 Ang.) [nred_pix, nfiber*nexp] 
;   bflux    - blue flux array (uncalibrated) [nblue_pix, nfiber*nexp]
;   rflux    - red flux array (uncalibrated) [nred_pix, nfiber*nexp]
;   bivar    - blue inverse variance [nblue_pix, nfiber*nexp]
;   rivar    - red inverse variance [nred_pix, nfiber*nexp]
;   best_exp - string containing the exposure ID of the "best" exposure --
;              this is used as the fiducial exposure against which the 
;              others are ratioed
;   plugtag  - plugmap-like structure which also contains the exposure ID of 
;              each spectrum [nfiber*nexp]
;   corrfile - names of output FITS files (1/exposure) [nexp]. The traditional
;              names are spFluxcorr-eeeeeeee-s.fits where "eeeeeeee" is the 
;              exposure ID and "s" is the spectrograph ID (1/2)
;
; OPTIONAL KEYWORDS:
;   diag  - if set show plots of individual fibers at the 3 S/N levels
;
; OUTPUTS:  
;   FITS files named in "corrfile" are generated containing the 4th order
;   legendere coefficients of the ratio of the best exposure and each of 
;   the others (on a fiber-to-fiber basis).  Plots are also generated.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Objects with high S/N are fit with a 4th order legendre
;   Objects with medium S/N are fit with a 3rd order legendre
;   Objects with Low S/N are assigned the median correction of the five nearest
;       fibers.  The zeropoint is then fit.
;   Bad fits -- those wildly deviant from the others -- are replaced with a 
;       zeropoint shift.
;
; EXAMPLES:
;
; BUGS:
;  The median corrections assigned to low S/N fibers are probably imperfect
;  because of the mix of point and extended sources on the plate.  Also
;  the shape of the correction may not be well fit by a low order polynomial
;  -- particularly if there is a mis-match in the dichroic region.
; 
;  Mask bits not currently returned!! (no record of hi/med/low S/N correction)
;  
; PROCEDURES CALLED:
;   djs_iterstat
;   djs_oplot
;   djs_plot
;   legend
;   mwrfits
;   pixelmask_bits()
;   splog
;   spmedian_rebin()
;   traceset2xy
;   xy2traceset
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   17-Oct-2000  Formerly fluxcorr_new -- written by S. Burles
;   09-Oct-2002  Revised by C. Tremonti to calibrate point sources to the 
;                smear and galaxies to the calib images.
;   12-Aug-2003  Changed by C. Tremonti to work on pairs of science images.
;                Split off "spmedian_rebin" & altered handling of medium and 
;                low S/N cases
;-
;------------------------------------------------------------------------------

pro frame_flux_tweak, bloglam, rloglam, bflux, rflux, bivar, rivar, $
    best_exp, plugtag, corrfile, diag = diag

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
   medloglam = [bwave,rwave] # replicate(1,nfiber) 
   wave = 10.0^medloglam[*,0]

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
      xy2traceset, medloglam, fitimg, corrset, ncoeff=4, /silent
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
        xy2traceset, medloglam[*,poly3], bestflux[*,poly3], polyset, $
            invvar=bestivar[*,poly3], ncoeff=4, inputfunc=sciflux[*,poly3], $
            lower = 3, upper = 3
        corrset.coeff[*,poly3] = polyset.coeff
      endif

      if npoly2 gt 0 then begin
        xy2traceset, medloglam[*,poly2], bestflux[*,poly2], polyset, $
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
           
        traceset2xy, corrset, medloglam, corrtemp
        xy2traceset, medloglam[*,poly1], bestflux[*,poly1], polyset, $
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
      ; Reject any bad corrections and replace with zeropoint shift only.
      ; Bad vectors are those with fit coefficients deviant by > 6-sigma.   
      ; Sigma is computed of the zero order coeff since this quantity is 
      ; robustly measured and a good indicator of plate quality

      ; Normalize coefficients 
      coef0 = corrset.coeff[0,*]  
      coef1 = corrset.coeff[1,*] / coef0 
      coef2 = corrset.coeff[2,*] / coef0 
      coef3 = corrset.coeff[3,*] / coef0 

      djs_iterstat, coef0, sigrej=6, mean=coef0mean, sigma=coef0sig 
      coef0 = coef0 / coef0mean
      coef0sig = (coef0sig / coef0mean) > 0.10
    
      isbad = (coef0 LT  (1 - 6*coef0sig) OR coef0 GT (1 + 6*coef0sig)) OR $
              (coef1 LT -6*coef0sig OR coef1 GT 6*coef0sig) OR $
              (coef2 LT -6*coef0sig OR coef2 GT 6*coef0sig) OR $
              (coef3 LT -6*coef0sig OR coef3 GT 6*coef0sig) 

      bad = where(isbad, nbad)
      if bad[0] NE -1 then begin
        splog, 'Warning: Large deviations in flux correction '
        splog, 'Warning: Replacing with zero point shift:', $
        string(bad + 1)
    
        if keyword_set(diag) then begin 
          !P.MULTI = [0, 1, 2]
          for ii = 0, nbad - 1 do begin
            ifib = bad[ii]
            plot, wave, corrtemp[*,ifib], yr=[0.0, 2.0], /nodata, $
              ytitle='Best Frame Flux / Science Frame Flux', title='Bad Vector'
            oplot, wave, bestflux[*,ifib]/sciflux[*,ifib], psym=6, syms=0.5, $
              thick=3
            djs_oplot, wave, corrtemp[*,ifib], color='red', thick=3
          endfor      
        endif
 
        xy2traceset, medloglam[*,bad], bestflux[*,bad], polyset, $
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

      traceset2xy, corrset, medloglam, corrimage

      ;---------------------
      ; Show individual fits to standards and examples of 3 S/N bins

      if keyword_set(diag) then begin
        !P.MULTI = [0, 1, 2]
        std = where(strmatch(plugtag[indx].objtype, '*_STD*'), nstd)
        for ii = 0, nstd - 1 do begin
          istd = std[ii]
          plot, wave, corrimage[*,istd], yr=[0.5, 1.5], /nodata, $
            ytitle='Best Frame Flux / Science Frame Flux', title='Standard Star'
          oplot, wave, bestflux[*,istd]/sciflux[*,istd], psym=6, syms=0.5, $
            thick=3
          djs_oplot, wave, corrimage[*,istd], color='red', thick=3
        endfor 
      
        for ii = 0, (10 < npoly2) - 1 do begin
          istd = poly2[ii]
          plot, wave, corrimage[*,istd], yr=[0.5, 1.5], /nodata, $
            ytitle='Best Frame Flux / Science Frame Flux', title='Poly 2'
          oplot, wave, bestflux[*,istd]/sciflux[*,istd], psym=6, $
            syms=0.5, thick=3
          djs_oplot, wave, corrimage[*,istd], color='red', thick=3
        endfor 

        for ii = 0, (10 < npoly1) - 1 do begin
          istd = poly1[ii]
          plot, wave, corrimage[*,istd], yr=[0.0, 2.0], /nodata, $
            ytitle='Best Frame Flux / Science Frame Flux', title='Poly 1'
          oplot, wave, bestflux[*,istd]/sciflux[*,istd], psym=6, $
            syms=0.5, thick=3
          djs_oplot, wave, corrimage[*,istd], color='red', thick=3
        endfor 
      endif

      ;-----------------
      ; Plot all low and medium S/N correction vectors

      !P.MULTI = [0, 1, 2]

      djs_plot, wave, corrimage, /nodata, yr=[0.4, 1.6], $
        xr=[min(wave)-100, max(wave)+100], /xstyle, /ystyle, $
        xtitle='\lambda [\AA]', ytitle='Best Frame Flux / Science Frame Flux', $
        title = 'Science: ' + expid[ifile] + ' Best: ' + best_exp 


      for iobj=0, npoly2 -1 do $
        djs_oplot, wave, corrimage[*,poly2[iobj]], color='blue', nsum=5
      for iobj=0, npoly1 -1 do $
        djs_oplot, wave, corrimage[*,poly1[iobj]], color='magenta', nsum=5

      legend, ['High S/N', 'Med S/N', 'Low S/N', 'Standard Star'], psym=0, $
              thick = 3, color=djs_icolor(['green', 'blue', 'magenta', 'red'])

      ;-----------------
      ; Plot all high S/N correction vectors + standards

      djs_plot, wave, corrimage, /nodata, yr=[0.4, 1.6], $
        xr=[min(wave)-100,max(wave)+100], /xstyle, /ystyle, $
        xtitle='\lambda [\AA]', ytitle='Best Frame Flux / Science Frame Flux'

      for iobj=0, npoly3 -1 do $
        djs_oplot, wave, corrimage[*,poly3[iobj]], color='green', nsum=5

      std = where(strmatch(plugtag[indx].objtype, '*_STD*'), nstd)
      for iobj=0, nstd -1 do $
        djs_oplot, wave, corrimage[*,std[iobj]], color='red', nsum=5, thick=2

      !P.MULTI = 0
   endfor
   return
end

