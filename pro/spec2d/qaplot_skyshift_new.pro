;+
; NAME:
;   qaplot_skyshift_new
;
; PURPOSE:
;   Generate the QA histogram plots of the sky line positions
;
;
; CALLING SEQUENCE:
;   qaplot_skyshift_new, wset, xsky, skywaves, skyshift, objhdr, [title= ]
;
; INPUTS:
;   wset       - Wavelength solution from arc line spectra
;   xsky       - Pixel position of sky lines [NFIBER,NLINE]
;   skywaves   - Wavelengths of sky lines used for computing shifts (Ang)
;   skyshift   - Shifts to apply to sky lines in pix [NROW,NTRACE]
;   objhdr     - Header of the spFrame file
;
; OPTIONAL KEYWORDS:
;   title      - TITLE of plot
;
; OUTPUTS:
;   Output plots only
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_oplot
;   djs_plot
;   traceset2pix
;   GaussFit
;
; REVISION HISTORY:
;   24-May-2021 Written by Hector Ibarra, UIUC
;-
;------------------------------------------------------------------------------
pro qaplot_skyshift_new, wset, xsky, skywaves, skyshift, objhdr, title=title
   dims = size(xsky, /dimens)
   nfiber = dims[0]
   nskyline = n_elements(skywaves)
   ysky = dindgen(nfiber) # (dblarr(nskyline)+1)
   xpredict = transpose( traceset2pix(wset, alog10(skywaves)) )
   scalefac = 1.0 ; Scale factor in pixel shift per tic mark on the X axis
   ;----------
   binsize=0.1
   xplot2 =  (xsky - xpredict) / scalefac
   xplot3 =  (xsky - xpredict - skyshift) / scalefac
   max1=max(xplot2)
   max2=max(xplot3)
   min1=min(xplot2)
   min2=min(xplot3)
   if max1 lt max2 then begin
     max_val=max2+0.2
   endif else begin
     max_val=max1+0.2
   endelse
   if min1 lt min2 then begin
     min_val=min1-0.2
   endif else begin
     min_val=min2-0.2
   endelse
   max_val=1.5
   min_val=-1.5
   tp=0
   for i=0, nskyline-1 do begin
    indx_val = i mod 4
    if indx_val eq 0 then begin
      !p.multi = [0,2,4]
      tp=0
    endif
    ;Plot sky line positions before shift
    lambd=skywaves[i]
    tit = string(lambd, format='("Shift line (",f6.1,")")')
    lin = string(lambd, format='(f6.1)')
    h = HISTOGRAM(xplot2[*,i], BINSIZE=binsize, LOCATIONS=loc)
    plot, loc + (binsize / 2.0), h, xrange=[min_val,max_val], $
    title=title+' '+tit+' before shift', xtitle='pixel shift (xsky-xpredict)', ytitle='N', $
    charsize=0.8, /xstyle, /ystyle
    binCenters = loc + (binsize / 2.0)
    if n_elements(h) gt 3 then begin
       yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
       djs_oplot, binCenters, yfit, COLOR='red', THICK=2 
       maxfit = String(coeff[0], FORMAT='(F0.3)')
       centerfit = String(coeff[1], FORMAT='(F0.3)')
       sigma = String(coeff[2], FORMAT='(F0.3)')
       avg=string(mean(xplot2[*,i]), FORMAT='(F0.3)')
       std=string(stdev(xplot2[*,i]), FORMAT='(F0.3)')
       ;xyouts, 0.35, 0.15+0.25*(3-indx_val), /NORMAL, 'Maximum: ' + maxfit, charsize=1.0;, COLOR='navy'
       xyouts, 0.32, 0.10+0.25*(3-indx_val)+0.05, /NORMAL, 'Center: ' + centerfit, charsize=1.0;, COLOR='navy'
       xyouts, 0.32, 0.05+0.25*(3-indx_val)+0.05, /NORMAL, 'Sigma: ' + sigma, charsize=1.0;, COLOR='navy'
       ; Save the shifts in the header file
       splog, 'Average and STD  position of '+lin+' before shift: '+avg+' , '+std+' (pixels)'
       sxaddpar, objhdr, 'AVGBSH'+strtrim(string(i),2), avg, $
        lin+' line average position before shift (pixels)'
       sxaddpar, objhdr, 'STDBSH'+strtrim(string(i),2), std, $
        lin+' line std position before shift (pixels)'
       sxaddpar, objhdr, 'SIGBSH'+strtrim(string(i),2), sigma, $
        lin+' line gaussian fit sigma position before shift (pixels)'
       sxaddpar, objhdr, 'CENBSH'+strtrim(string(i),2), centerfit, $
        lin+' line gaussian fit center position before shift (pixels)'
    endif
    ;----------
    ; Plot sky line positions after shift
    h = HISTOGRAM(xplot3[*,i], BINSIZE=binsize, LOCATIONS=loc)
    plot, loc + (binsize / 2.0), h, $
    xrange=[min_val,max_val], $
    title=title+' '+tit+' after shift', xtitle='pixel shift (xsky-xpredict-skyshift)', ytitle='N', $
    charsize=0.8, /xstyle, /ystyle;, /fill
    binCenters = loc + (binsize / 2.0)
    if n_elements(h) gt 3 then begin
       yfit = GaussFit(binCenters, h, coeff, NTERMS=3)
       djs_oplot, binCenters, yfit, COLOR='green', THICK=2
       maxfit = String(coeff[0], FORMAT='(F0.3)')
       centerfit = String(coeff[1], FORMAT='(F0.3)')
       sigma = String(coeff[2], FORMAT='(F0.3)')
       avg=string(mean(xplot3[*,i]), FORMAT='(F0.3)')
       std=string(stdev(xplot3[*,i]), FORMAT='(F0.3)')
       ;xyouts, 0.85, 0.15+0.25*(3-indx_val), /NORMAL, 'Maximum: ' + maxfit, charsize=1.0;, COLOR='navy'
       xyouts, 0.82, 0.10+0.25*(3-indx_val)+0.05, /NORMAL, 'Center: ' + centerfit, charsize=1.0;, COLOR='navy'
       xyouts, 0.82, 0.05+0.25*(3-indx_val)+0.05, /NORMAL, 'Sigma: ' + sigma, charsize=1.0;, COLOR='navy' 
       ; Save the shifts in the header file
       splog, 'Average and STD position of '+lin+' after shift: '+avg+' , '+std+' (pixels)'
       sxaddpar, objhdr, 'AVGASH'+strtrim(string(i),2), avg, $
        lin+' line average position after shift (pixels)'
       sxaddpar, objhdr, 'STDASH'+strtrim(string(i),2), std, $
        lin+' line std position after shift (pixels)'
       sxaddpar, objhdr, 'SIGASH'+strtrim(string(i),2), sigma, $
        lin+' line gaussian fit sigma position after shift (pixels)'
       sxaddpar, objhdr, 'CENASH'+strtrim(string(i),2), centerfit, $
        lin+' line gaussian fit center position after shift (pixels)'
    endif
    if indx_val eq 3 then begin
      !p.multi = 0
      tp=1
    endif
   endfor  
   if tp eq 0 then !p.multi = 0
   return
end
