function median_rebin, flux, ivar, loglam, range, mask=mask

   nr = (size(range))[2]
   ntrace = (size(flux))[2]

   fit = fltarr(nr, ntrace)
   mask = fit

   for itrace=0,ntrace-1 do begin
     for irange = 0, nr -1 do begin
        inside = where(loglam[*,itrace] GE range[0,irange] $
                 AND loglam[*,itrace] LT range[1,irange], ninside)
        if ninside GT 0 then begin
           good = where(ivar[inside,itrace] GT 0, ngood)
       
           if ngood GT 1 then begin 
              fit[irange, itrace] = median(flux[inside[good],itrace])
              if ngood GT 0.5 * ninside then $
                   mask[irange, itrace]  = total(ivar[inside[good],itrace])
           endif
        endif
     endfor
   endfor

return, fit
end
              
pro fluxcorr_new, bsmearfile, rsmearfile, bscifile, rscifile, corrfile

   ;----------
   ; Read the plug-map file (for identifying sky fibers)

   bsmearflux = mrdfits(bsmearfile,0)
   bsmearivar = mrdfits(bsmearfile,1)
   bsmearmask = mrdfits(bsmearfile,2)
   bsmearset  = mrdfits(bsmearfile,3)
   plugmap    = mrdfits(bsmearfile,5)
   traceset2xy, bsmearset, xx, bsmearloglam

   b1 = findgen(60)*4.0e-3 + 3.568
   b2 = findgen(60)*4.0e-3 + b1[1]

   brange = transpose([[b1],[b2]])
   bwave = djs_median(brange,1)

   bfit = median_rebin(bsmearflux, bsmearivar, bsmearloglam, brange, $
           mask=bmask)



   rsmearflux = mrdfits(rsmearfile,0)
   rsmearivar = mrdfits(rsmearfile,1)
   rsmearmask = mrdfits(rsmearfile,2)
   rsmearset  = mrdfits(rsmearfile,3)
   traceset2xy, rsmearset, xx, rsmearloglam

   r1 = findgen(54)*4.0e-3 + 3.756
   r2 = findgen(54)*4.0e-3 + r1[1]

   rrange = transpose([[r1],[r2]])
   rwave = djs_median(rrange,1)

   rfit = median_rebin(rsmearflux, rsmearivar, rsmearloglam, rrange, $
           mask=rmask)

   bsn = djs_median(bsmearflux * sqrt(bsmearivar),1)
   rsn = djs_median(rsmearflux * sqrt(rsmearivar),1)
   smearsnmed = transpose([[bsn],[rsn]])
 
   ;---------------------------------------------------------
   ;  now put blue and red together
   ;

   smearflux = [bfit,rfit]
   smearivar = [bmask, rmask]
   nfiber = (size(smearflux,/dimens))[1]
   wave = [bwave,rwave] # replicate(1,nfiber)

   ncoeff = 3 ; ???

   nfiles = n_elements(corrfile)

   for ifile = 0, nfiles - 1 do begin
 
   if (bsmearfile EQ bscifile[ifile] AND $
       rsmearfile EQ rscifile[ifile]) then begin
      ;----------
      ; Special case: If the science image and smear image is the same,
      ; then force their ratio to be unity.

      fitimg = smearflux*0.0 + 1.0
      xy2traceset, wave, fitimg, corrset, ncoeff=ncoeff

      corrset.coeff[0,*] = 1.0
      corrset.coeff[1:*,*] = 0.0

   endif else begin

     bsciflux = mrdfits(bscifile[ifile],0)
     bsciivar = mrdfits(bscifile[ifile],1)
     bscimask = mrdfits(bscifile[ifile],2)
     bsciset  = mrdfits(bscifile[ifile],3)
     traceset2xy, bsciset, xx, bsciloglam

     bfit = median_rebin(bsciflux, bsciivar, bsciloglam, brange, $
           mask=bmask)

     rsciflux = mrdfits(rscifile[ifile],0)
     rsciivar = mrdfits(rscifile[ifile],1)
     rscimask = mrdfits(rscifile[ifile],2)
     rsciset  = mrdfits(rscifile[ifile],3)
     traceset2xy, rsciset, xx, rsciloglam

     rfit = median_rebin(rsciflux, rsciivar, rsciloglam, rrange, $
           mask=rmask)
     bsn = djs_median(bsciflux * sqrt(bsciivar),1)
     rsn = djs_median(rsciflux * sqrt(rsciivar),1)
     scisnmed = transpose([[bsn],[rsn]])

     sciflux = [bfit,rfit]
     sciivar = [bmask,rmask]

     ;-------------------------------------------------------------
     ;   3 levels of S/N
     ;    level 1:  High S/N, independent coefficient solution
     ;    level 2:  Med  S/N, Scaled Spectrophoto solution
     ;    level 3:  Low  S/N, median spectrophoto solution

     highsn = where(scisnmed[0,*] GT 2.5 $
               AND  scisnmed[1,*] GT 5.0 $
               AND  smearsnmed[0,*] GT 1.0 $
               AND  smearsnmed[1,*] GT 1.0)

     if highsn[0] EQ -1 then return

     spectrophoto = where(strtrim(plugmap[highsn].objtype,2) EQ $
                           'SPECTROPHOTO_STD', nspectrophoto)

     if spectrophoto[0] EQ -1 then begin
        splog, "ABORT: No spectrophoto with high S/N"
        return
     endif

     spectrophoto = highsn[spectrophoto]

     fibersn = lonarr(nfiber) + 3
     fibersn[highsn] = 1

     medsn = where((scisnmed[0,*] GT 1.0 $
                OR scisnmed[1,*] GT 2.0 ) $
               AND  (smearsnmed[0,*] GT 0.2 $
                OR smearsnmed[1,*] GT 0.5) AND fibersn NE 1) 

     if medsn[0] NE -1 then fibersn[medsn] = 2

     xy2traceset, wave, smearflux, highsnset, $
        invvar=smearivar, ncoeff=ncoeff, inputfunc=sciflux, $
        lower = 3, upper = 3

     traceset2xy, highsnset, wave, highsnimage 

     medianset = highsnset
     finalset  = highsnset

     if nspectrophoto EQ 1 then begin
         mediancoeff = highsnset.coeff[*,spectrophoto] 
         splog, "WARNING: Only 1 spectrophoto with high S/N"
     endif else mediancoeff = djs_median(highsnset.coeff[*,spectrophoto],2) 

     medianset.coeff = mediancoeff # replicate(1,nfiber)

     traceset2xy, medianset, wave, medianimage

     smearnormflux = smearflux  
     smearnormivar = smearivar
     divideflat, smearnormflux, smearnormivar, medianimage 

     if medsn[0] NE -1 then begin     
 
       xy2traceset, wave, smearnormflux, medsnset, $
          invvar=smearnormivar, ncoeff=1, inputfunc=sciflux, $
          lower = 3, upper = 3

       finalset.coeff[*,medsn] = mediancoeff # medsnset.coeff[0,medsn]
     endif

     lowsn = where(fibersn EQ 3, nlowsn)
     if nlowsn GT 0 then $
         finalset.coeff[*,lowsn] = mediancoeff # replicate(1,nlowsn)

     corrset = finalset
     traceset2xy, corrset, wave, corrimage

   endelse

   mwrfits, corrset, corrfile[ifile], /create

   endfor

   return
end
