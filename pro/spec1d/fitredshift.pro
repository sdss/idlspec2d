;+
; NAME:
;   fitredshift
;
; PURPOSE:
;   Find the most significant correlation peak (soon to be peaks)
;     return optional keywords in units of pixels
;
; CALLING SEQUENCE:
;    fitredshift, fluxfft, starfft, pad=pad, ngoodpixels=ngoodpixels, $
;          nsearch=nsearch, zfit=z, z_err=z_err, $
;          veldispfit=veldisp, veldisp_err=veldisp_err, doplot=doplot
;
; INPUTS:
;   fluxfft    - complex fft of prepared galaxy spectrum
;   starfft    - complex fft of stellar template
;   pad        - number of padding pixels used in preparing ffts
;
; OPTIONAL KEYWORDS:
;   ngoodpixels- long array indicating the number of pixels overlap at each lag
;   nsearch    - number of peaks to search, almost always only 1 is searched
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
;
; EXAMPLES:
;
; BUGS:
;
;	Hardwired exclusion of blueshifts greater than 100 pixels
;        this helps the noisiest cases
;
;	Need to compile GAUSSFIT beforehand, which includes the 
;        fitting function GAUSS_FUNCT 
;
; PROCEDURES CALLED:
;   findmaxarea
;   curvefit
;   gauss_funct
;
; REVISION HISTORY:
;   25-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;
;
;;	This is used to locate the 20 highest peaks, and measure
;       the peak, center, and area of each.  xcen and peak
;       contain the best velues.
;
;
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

	 ; find bounds of positive deivation with cap at +/- 3 pixels
         lowerbound = max([where(ruin LT 0 AND x LT velcen),velcen-4]) + 1
         upperbound = min([where(ruin LT 0 AND x GT velcen),velcen+4]) - 1
         area[i] = total(ruin[lowerbound:upperbound]) 

	 dev[i] = stddev([look[velcen-100:lowerbound-1],$
                          look[upperbound+1:velcen+100]])
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
   
pro fitredshift, fluxfft, starfft, pad=pad, ngoodpixels=ngoodpixels, $
          nsearch=nsearch, zfit=z, z_err=z_err, $
          veldispfit=veldisp, veldisp_err=veldisp_err, doplot=doplot

        if (NOT keyword_set(nsearch)) then nsearch=5

        corr =  float(fft(fluxfft * conj(starfft),1))
        
        nx = n_elements(corr)
        if (NOT keyword_set(pad)) then pad = nx/2

	if (NOT keyword_set(ngoodpixels)) then begin
          denom = pad*(1.0 - abs(findgen(nx)-pad)/pad)
          compensate = pad/(denom > 1.0) 
        endif else $
          compensate = max(ngoodpixels)/(ngoodpixels > 1.0) 

 
        corr = shift(corr,pad)
      ;  compensate = shift(compensate,pad)
 

      ;  This loop finds the redshift by searching the nsearch highest peaks

	good = 0
        x = lindgen(nx)
        for i=0,nsearch-1 do begin

           look = corr*sqrt(compensate)	

    ;!!!!!!! danger here !!!!!!!
    ;  I'm not allowing for blueshifts larger than 100 pixels 
    ;  Do this by zeroing out look from 0 to pad - 101
    ; !!!!!!!!!!!!!!!!!!!!!!!!!!

	  look[0:pad-101] = 0.0

	   findmaxarea, look, velcen, peak, cen=cen, area=area, pks=pks
           ; peak = max(look,velcen)

           ; let xtemp be centered about velcen for all corr values above
           ; 0.0

	   lowerbound = max(where(corr LT 0 AND x LT velcen))
	   upperbound = min(where(corr LT 0 AND x GT velcen))
           xtemp = x[lowerbound:upperbound]-velcen
           parabola = poly_fit(xtemp, corr[xtemp+velcen], 2)

           yfit = poly(xtemp, parabola)


           if (parabola[2] GE 0.0) then begin
             print, 'peak is not well fit at ', velcen
             corr[xtemp+velcen] = 0.0
           endif else if (total(corr[xtemp+velcen]) LT 0.0) then begin
             print, 'total corr is less than zero at ', velcen
             corr[xtemp+velcen] = 0.0
           endif else if (total(yfit) LT 0.0) then begin
             print, 'total fit is less than zero at ', velcen
             corr[xtemp+velcen] = 0.0
           endif else begin
             i=5
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
        guesssig = sqrt(total((xtemp-xcen)^2 * ytemp)/total(ytemp))


        left = long(velcen + xcen - 1) - lindgen(100)
        right = long(velcen + xcen + 1) + lindgen(100)
        asig = stddev(corr[left]-corr[right])/sqrt(2)


	; Here's my attempt to fit a gaussian to the correlation peak
        ;  The main problem here is to decide where the baseline of the
        ;   gaussian falls.  FWHM and sigma depend on where the gaussian fit
        ;   goes to zero.  Very troubling.

        if (asig GT 0) then weights = xtemp*0.0 + 1.0/asig^2 

        a = [height, xcen, guesssig]
        gaussf = curvefit(xtemp,corr[xtemp+velcen] + height/3.0, weights, $
                         a, gausserrors, function_name = "GAUSS_FUNCT")

        gaussf = gaussf - height/3.0

	;  4 parameter gaussian fit below does not always converge

;        gaussf = gaussfit(xtemp,corr[xtemp+velcen], a, nterms=4, $
;                           estimates = [height, xcen, guesssig, -height/3.0])
;        print, a[3], height/3.0

	if (keyword_set(doplot)) then begin
             window,0
             plot, x-velcen, corr, ps =10, xr=[-20,20], $
              title='Best correlation peak w/fits (Green:gauss, Red: Parabola)'
             oplot, xtemp, poly(xtemp, parabola), color=500
             oplot, xtemp, gaussf, color=50000
        endif


        z = velcen + a[1] - pad
        z_err = gausserrors[1]
        veldisp = a[2]
        veldisp_err = gausserrors[2]

        return
end
