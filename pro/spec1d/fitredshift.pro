;+
; NAME:
;   fitredshift
;
; PURPOSE:
;   Find the most significant correlation peak (soon to be peaks)
;     return optional keywords in units of pixels
;
; CALLING SEQUENCE:
;    fitredshift, fluxfft, starfft, $
;     [ nsearch=, zfit=, z_err=, veldispfit=, veldisp_err=, /doplot ]
;
; INPUTS:
;   fluxfft    - complex fft of prepared galaxy spectrum
;   fluxerr    - error vector for flux (only test if it is zero)
;   starfft    - complex fft of stellar template
;   starerr    - error vector for flux (only test if it is zero)
;
; OPTIONAL KEYWORDS:
;   nsearch    - number of peaks to search, almost always only 1 is searched
;   zmin       - Minimum z (in pixels) to allow (should be < 0)
;   doplot     - plot the correlation peak and fits in an xwindow
;
; OPTIONAL OUTPUTS:
;   zfit       - best fit z, this should evolve into an array of z's
;                  with accompanying z_errs and z_confidences
;   z_err      - centroid errors from gaussfit of correlation peak
;   veldisp    - sigma of cross-correlation peak
;   veldisp_err- error on sigma
;
; COMMENTS:
;
;   Use doplot keyword to see how well peak is being fit
;   Still need to work on exact selection criteria for MOST significant peak
;     or even better: measure all peaks with probability > 1%     
;
; EXAMPLES:
;
; BUGS:
;
; Hardwired exclusion of blueshifts greater than 100 pixels
;        this helps the noisiest cases
;
; PROCEDURES CALLED:
;   findmaxarea
;   curvefit
;   gauss_funct
;
; REVISION HISTORY:
;   25-Mar-2000  Written by S. Burles, FNAL
;   26-Jun-2000  D. Finkbeiner - modified to properly weight corr
;                                vector for variable overlap
;   26-Jun-2000 (v 1.5) altered algorithm to set peak search boundary
;       at the half height of the peak, rather than zero.  This fixes
;       cases where there is a close double peak. 
;-
;------------------------------------------------------------------------------
; This routine locates the 20 highest peaks, and measures
; the peak, center, and area of each.  xcen and peak
; contain the best velues.

pro findmaxarea, look, xcen, peak, maxarea, cen=cen, area=area, pks=pks

   ruin = look
   nx = n_elements(look)
   x = findgen(nx)
   area = fltarr(20)
   dev = fltarr(20)
   cen = fltarr(20)
   pks = fltarr(20)

   for i=0, 19 do begin

      pks[i] = max(ruin,velcen)

      if pks[i] GT 0.0 then begin

         ; find bounds of positive deviation with cap at +/- 3 pixels
         lowerbound = max([where(ruin LT 0 AND x LT velcen),velcen-4]) + 1
         upperbound = min([where(ruin LT 0 AND x GT velcen),velcen+4]) - 1
         area[i] = total(ruin[lowerbound:upperbound]) 

         dev[i] = stddev([look[velcen-100:lowerbound-1],$
                          look[upperbound+1:velcen+100]], /double)
         cen[i] = velcen
         ruin[lowerbound:upperbound] = 0.0
      endif else i=20
   endfor

   meandev = mean(dev)
   dev = dev/meandev

   maxarea = max(area,areaplace)
   maxpks = max(pks,pksplace)

   place = pksplace
   if (areaplace NE pksplace) then begin
      print, 'Max area position does not equal Max peak position', $
                cen[areaplace], cen[pksplace]
      print, 'Using Area + 1.5*peak to decide'
      maxcomb = max(area + 1.5*pks,place)
   endif
     
   xcen = cen[place]
   peak = pks[place]
   return
end

;------------------------------------------------------------------------------ 
pro fitredshift, fluxfft, fluxerr, starfft, starerr, $
 nsearch=nsearch, zmin=zmin, zfit=zfit, z_err=z_err, $
 veldispfit=veldisp, veldisp_err=veldisp_err, zconf=zconf, doplot=doplot

   if (NOT keyword_set(nsearch)) then nsearch = 5

; returned value defaults
   zfit = 0.
   z_err = 999.
   zconf = 0.
   veldisp = 0.
   veldisp_err = 999.

; check dimensions - ASSUME everything is already padded to 2^N
   IF stddev([n_elements(fluxfft), n_elements(fluxerr),  $
             n_elements(starfft), n_elements(starerr)]) NE 0 THEN BEGIN 
       help, fluxfft, fluxerr, starfft, starerr 
       message, 'dimensions do not match!'
   ENDIF 

   IF (total(abs(fluxfft)) EQ 0) OR (total(abs(starfft)) EQ 0) THEN BEGIN 
       print, 'FITREDSHIFT:  FAILED - array full of zeros'
       return
   ENDIF 


; The following compiles GAUSS_FUNCT, which is in the routine GAUSSFIT
   junk = gaussfit(findgen(7),findgen(7))

; compute correlation of object flux to stellar template
   corr = float(fft(fluxfft * conj(starfft),/inverse))
   nx = n_elements(corr)
   pad = nx / 2

; shift corr vector so there is no discontinuity at zero lag. 
   corr = shift(corr, pad)

; Need to fill an array of length fluxfilt which records the number
; of good pixels cross-correlated between star and galaxy as a function
; of shift.  We construct this from the error vectors. 
  
   denom = shift(double( $
         fft(fft(fluxerr NE 0)*conj(fft(starerr NE 0)), /inverse)), pad) *nx

   reweight = pad / (denom > 10.0) ; never use less than 10 pix


; This loop finds the redshift by searching the nsearch highest peaks

   good = 0
   x = lindgen(nx)

   for i=0, nsearch-1 do begin

      newcorr = corr * sqrt(reweight) ; A hack!!!???

      if (keyword_set(zmin)) then begin
; ???
         newcorr[0:pad+zmin-1] = 0.0
      endif

      findmaxarea, newcorr, velcen, peak, cen=cen, area=area, pks=pks

      ; Let xtemp be centered about velcen for all newcorr values above 0.0

      lowerbound = max(where(newcorr LT peak/2 AND x LT velcen))
      upperbound = min(where(newcorr LT peak/2 AND x GT velcen))
      xtemp = x[lowerbound:upperbound] - velcen
      if (!version.release LT '5.4') then begin
         parabola = poly_fit(xtemp, newcorr[xtemp+velcen], 2, yfit)
      endif else begin
         parabola = poly_fit(double(xtemp), double(newcorr[xtemp+velcen]), 2, $
          yfit=yfit)
      endelse

      if (parabola[2] GE 0.0) then begin
          print, 'peak is not well fit at ', velcen
          newcorr[xtemp+velcen] = 0.0
      endif else if (total(newcorr[xtemp+velcen]) LT 0.0) then begin
          print, 'total newcorr is less than zero at ', velcen
          newcorr[xtemp+velcen] = 0.0
      endif else if (total(yfit) LT 0.0) then begin
          print, 'total fit is less than zero at ', velcen
          newcorr[xtemp+velcen] = 0.0
      endif else begin
          i = nsearch
          good = 1
; let's attempt to fit a gaussian with no background terms
      endelse

   endfor

   if (NOT good) then begin
      print, 'No good peaks found'
      return
   endif

   xcen = (-0.5 * parabola[1]/parabola[2])
   height = (poly([xcen], parabola))[0]

   ytemp = yfit > 0.0
   guesssig = sqrt(total((xtemp-xcen)^2 * ytemp) / total(ytemp))

   left = long(velcen + xcen - 1) - lindgen(100)
   right = long(velcen + xcen + 1) + lindgen(100)
   asig = stddev(newcorr[left]-newcorr[right], /double) / sqrt(2.)

   ; Here's my attempt to fit a gaussian to the correlation peak
   ; The main problem here is to decide where the baseline of the
   ; gaussian falls.  FWHM and sigma depend on where the gaussian fit
   ; goes to zero.  Very troubling.

; new xtemp
   lowerbound = max(where(newcorr LT height/2 AND x LT velcen))-2
   upperbound = min(where(newcorr LT height/2 AND x GT velcen))+2
   xtemp2 = x[lowerbound:upperbound] - velcen



   if (asig GT 0) then weights = xtemp2*0.0 + 1.0/asig^2 $
    else weights = xtemp2*0.0 + 1.0

   base = -height/3. ; empirical guess of gaussian baseline

   a = double([height, xcen, guesssig, base])
;   gaussf = curvefit(xtemp,newcorr[xtemp+velcen] + height/3.0, weights, $
;    a+0D, gausserrors, function_name="GAUSS_FUNCT")
;   gaussf = gaussf - height/3.0

   gaussf = curvefit(xtemp2,newcorr[xtemp2+velcen], weights, $
    a, gausserrors, function_name="GAUSS_FUNCT", iter=iter)

   IF abs(a[1]) GT 3 THEN BEGIN 
       message, 'curvefit failed to converge!!!', /info
       print, 'Reverting to parabola solution'
       
   ENDIF 


;   if (keyword_set(doplot)) then begin
;      wset,0
      djs_plot, x-velcen, newcorr, ps =10, xr=[-40,40], $
       title='Best correlation peak w/fits (Green:gauss, Red: Parabola)'
      djs_oplot, xtemp, poly(xtemp, parabola), color='red'
      djs_oplot, xtemp2, gaussf, color='green'
;   endif

   zfit = velcen + a[1] - pad
   z_err = gausserrors[1]
   veldisp = a[2]
   veldisp_err = gausserrors[2]

;   zconf = a[0] ; height of gaussian

   twopiei = 2.0 * !dpi * complex(0.0,1.0)
   knums = fft_wavenums(n_elements(starfft))
   phase = exp( - twopiei * knums * zfit)
   model = double(fft(starfft*phase,/inv))
   gal = double(fft(fluxfft, /inv))
   
   starmask = where(starerr NE 0)
   chisq = total((model-gal)^2*fluxerr*starmask)/total(starmask)
   zconf = chisq

   return
end
