;+
; NAME:
;   fluxcorr_new
;
; PURPOSE:
;   Compute the frame-by-frame relative scaling for each object on a plate.
;
; CALLING SEQUENCE:
;   fluxcorr_new, bsmearfile, rsmearfile, bscifile, rscifile, corrfile
;
; INPUTS:
;   bsmearfile - spFrame fits file chosen as blue smear image
;   rsmearfile - spFrame fits file chosen as red smear image
;   bscifile   - spFrame fits file(s) blue science image
;   rscifile   - spFrame fits file(s) red science image
;   corrfile   - Fits file to output flux correction vectors.
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;   fluxcorr_new is used to calculate and write to file a low order
;    polynomial function which registers the flux in the science exposures
;    to the median flux levels in the smear exposures.
;
;   Based on S/N, the produced function per fiber can be given by one
;    of three methods.  For the highest S/N fibers, typically 10\% of the
;    fibers including spectrophoto stds, the full 3rd order fit is done.
;   For medium S/N fibers, only a one parameter scaling of the best 
;     spectrophoto correction is produced.
;   For the lowest S/N fibers, the spectrophoto correction is used without
;     any scaling.  Beware, this means a stellar spectrophoto correction 
;     is applied to all low S/N fibers, even if they are not
;     seeing-limited targets.
;
; EXAMPLES:
;
; BUGS:
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
;
; INTERNAL SUPPORT ROUTINES:
;   median_rebin():  Used to rebin spectra in large wavelength blocks
;                    passed in parameter range 
;
; REVISION HISTORY:
;   17-Oct-2000  Written by S. Burles
;     
;-
;------------------------------------------------------------------------------
function median_rebin, flux, ivar, loglam, range, mask=mask, sigrej=sigrej

   if NOT keyword_set(sigrej) then sigrej = 20.0

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

return, fit
end
              
;------------------------------------------------------------------------------
pro fluxcorr_new, bsmearfile, rsmearfile, bscifile, rscifile, corrfile

   splog, 'Smear image blue=', bsmearfile
   splog, 'Smear image red= ', rsmearfile

   ;----------
   ; Read the plug-map file (for identifying sky fibers)

   bsmearflux = mrdfits(bsmearfile,0,bhdr)
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



   rsmearflux = mrdfits(rsmearfile,0,rhdr)
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
   badval = ( fibermask_bits('NOPLUG')   OR fibermask_bits('BADTRACE') $
       OR fibermask_bits('BADFLAT')  OR fibermask_bits('BADARC')   $
       OR fibermask_bits('NEARWHOPPER') )
   qsmear =  (bsmearmask[0,*] OR rsmearmask[0,*]) AND badval

 
   ;---------------------------------------------------------
   ;  now put blue and red together
   ;

   smearflux = [bfit,rfit]
   smearivar = [bmask, rmask]
   nfiber = (size(smearflux,/dimens))[1]

   ; --------------------------------------------------------
   ;  Create basic smear mask
   ;
   smearmask = lonarr(nfiber) 
   if (strtrim(sxpar(bhdr, 'FLAVOR'),2) EQ 'smear' AND $
       strtrim(sxpar(rhdr, 'FLAVOR'),2) EQ 'smear') then $
         smearmask = smearmask OR pixelmask_bits('SMEARIMAGE')

   ;----------


   wave = [bwave,rwave] # replicate(1,nfiber)

   ncoeff = 3 ; ???

   nfiles = n_elements(corrfile)

   for ifile=0, nfiles-1 do begin

      splog, 'Fluxing science image blue=', bscifile[ifile]
      splog, 'Fluxing science image red=', rscifile[ifile]
 
      fitimg = smearflux*0.0 + 1.0
      xy2traceset, wave, fitimg, corrset, ncoeff=ncoeff

      corrset.coeff[0,*] = 1.0
      corrset.coeff[1:*,*] = 0.0

      ;  if they are the same, set both high and med SN flagbits

      thismask  = smearmask OR pixelmask_bits('SMEARHIGHSN')
      thismask  = thismask  OR pixelmask_bits('SMEARMEDSN')
      qsci = lonarr(nfiber)

   ; Special case: If the science image and smear image is the same,
   ; then force their ratio to be unity.

     if (NOT (bsmearfile EQ bscifile[ifile] AND $
          rsmearfile EQ rscifile[ifile])) then begin

       thismask  = smearmask 
       bfit = bwave # replicate(0,nfiber)
       bmask = bwave # replicate(0,nfiber)
       bsn = fltarr(nfiber) 
       checkblue = 0

       if bscifile[ifile] NE '' then begin
         bsciflux = mrdfits(bscifile[ifile],0)
         bsciivar = mrdfits(bscifile[ifile],1)
         bscimask = mrdfits(bscifile[ifile],2)
         bsciset  = mrdfits(bscifile[ifile],3)
         traceset2xy, bsciset, xx, bsciloglam

         bfit = median_rebin(bsciflux, bsciivar, bsciloglam, brange, $
             mask=bmask)
         bsn = djs_median(bsciflux * sqrt(bsciivar),1)
         checkblue = 1
         qsci = qsci OR ( bscimask[0,*] AND badval)

       endif else $
         splog, 'WARNING: no blue science frame for ', corrfile[ifile]
       
 
       rfit = rwave # replicate(0,nfiber)
       rmask = rwave # replicate(0,nfiber)
       rsn = fltarr(nfiber) 
       checkred = 0

       if rscifile[ifile] NE '' then begin

         rsciflux = mrdfits(rscifile[ifile],0)
         rsciivar = mrdfits(rscifile[ifile],1)
         rscimask = mrdfits(rscifile[ifile],2)
         rsciset  = mrdfits(rscifile[ifile],3)
         traceset2xy, rsciset, xx, rsciloglam

         rfit = median_rebin(rsciflux, rsciivar, rsciloglam, rrange, $
           mask=rmask)
         rsn = djs_median(rsciflux * sqrt(rsciivar),1)
         checkred = 1
         qsci = qsci OR ( rscimask[0,*] AND badval)

       endif else $
         splog, 'WARNING: no red science frame for ', corrfile[ifile]

       scisnmed = transpose([[bsn],[rsn]])

       sciflux = [bfit,rfit]
       sciivar = [bmask,rmask]

     ;-------------------------------------------------------------
     ;   3 levels of S/N
     ;    level 1:  High S/N, independent coefficient solution
     ;    level 2:  Med  S/N, Scaled Spectrophoto solution
     ;    level 3:  Low  S/N, median spectrophoto solution

       highsn = where((scisnmed[0,*] GT 2.5 OR checkblue EQ 0) $
               AND  (scisnmed[1,*] GT 5.0 OR checkred EQ 0) $
               AND  smearsnmed[0,*] GT 1.0 $
               AND  smearsnmed[1,*] GT 1.0 $
               AND  qsci EQ 0 AND qsmear EQ 0)

       spectrophoto = -1L
       fibersn = lonarr(nfiber) + 3

       if highsn[0] NE -1 then begin

         spectrophoto = where(strtrim(plugmap[highsn].objtype,2) EQ $
                           'SPECTROPHOTO_STD', nspectrophoto)

         if spectrophoto[0] EQ -1 then $
              spectrophoto = where(strtrim(plugmap[highsn].objtype,2) EQ $
                           'REDDEN_STD', nspectrophoto)

         fibersn[highsn] = 1
       endif

       if spectrophoto[0] EQ -1 then begin
             splog, "WARNING: No spectrophoto with high S/N"
             corrset.coeff = 0.0
       endif else begin

         spectrophoto = highsn[spectrophoto]


         medsn = where((scisnmed[0,*] GT 1.0 $
                OR scisnmed[1,*] GT 2.0 ) $
               AND  (smearsnmed[0,*] GT 0.2 $
                OR smearsnmed[1,*] GT 0.5) AND fibersn NE 1 ) 

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
         divideflat, smearnormflux, invvar=smearnormivar, medianimage 

         if medsn[0] NE -1 then begin     
 
           xy2traceset, wave, smearnormflux, medsnset, $
            invvar=smearnormivar, ncoeff=1, inputfunc=sciflux, $
            lower = 3, upper = 3, yfit=yy

           finalset.coeff[*,medsn] = mediancoeff # medsnset.coeff[0,medsn]
         endif

         lowsn = where(fibersn EQ 3, nlowsn)
         if nlowsn GT 0 then $
           finalset.coeff[*,lowsn] = mediancoeff # replicate(1,nlowsn)

         corrset = finalset
         traceset2xy, corrset, wave, corrimage

         totalsmearsn = total(smearsnmed,1)
         bad = where(djs_median(corrimage,1) GT $
              (4.0+totalsmearsn)*median(corrimage),nbad)
         if bad[0] NE -1 then begin
           splog, 'Warning: Large deviations in flux correction '
           splog, 'Warning: Replacing with median solution in fibers:', $
                string(bad + 1)
           corrset.coeff[*,bad] = mediancoeff # replicate(1,nbad)
           fibersn[bad] = 3
         endif
         ; Set maskbits and append to end of corrfile
         
         thismask = thismask OR (fibersn LE 2) * pixelmask_bits('SMEARMEDSN')
         thismask = thismask OR (fibersn EQ 1) * pixelmask_bits('SMEARHIGHSN')

       endelse
     endif

     mwrfits, corrset, corrfile[ifile], /create
     mwrfits, thismask, corrfile[ifile]

   endfor

   return
end
