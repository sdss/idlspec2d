;------------------------------------------------------------------------------
; This routine used to be used by VDISPFIT, but is no longer as of 10-Feb-02.
;------------------------------------------------------------------------------
pro findchi2min, x, chi2, minchi2, minsigma, errsigma, npts=npts, $
 deltachisq=deltachisq, nfine=nfine, $
 plottitle=plottitle, doplot=doplot, debug=debug

   if (NOT keyword_set(deltachisq)) then deltachisq = 1.0
   if (NOT keyword_set(nfine)) then nfine = 100

   ; spline chi2 curve and return desired values

   nx = n_elements(x)
   if (nx EQ 0) then return

   minx = min(x, max=maxx)
   y2 = spl_init(x, chi2)

   x2 = findgen(nx*nfine+1)/(nx*nfine) * (maxx - minx) + minx
   chifit = spl_interp(x, chi2, y2, x2)

   minchi2 = min(chifit,fitplace)

   minsigma = x2[fitplace] 

   ; attempt to scale chi^2 to get errors correct ???
  
   usethisdelta = deltachisq
   if (keyword_set(npts)) then begin
     scale =  minchi2/npts
     if (scale GT 1.0) then usethisdelta = deltachisq * scale
   endif

   range = where(chifit - minchi2 LT usethisdelta)

   if (range[0] EQ -1) then return

   uppersigma = x2[max(range)] - minsigma
   lowersigma = minsigma - x2[min(range)] 

   errsigma = 0.5*(uppersigma + lowersigma)

   if (keyword_set(doplot)) then begin
      plot, x, chi2, psym=1, $
       title=plottitle
      djs_oplot, x2, chifit
      djs_oplot, x2[range], chifit[range], color='red'

      ; Wait for a keystroke...
      if (keyword_set(debug)) then begin
         print, 'Press any key...'
         cc = strupcase(get_kbrd(1))
      endif

   endif

   return
end    
;------------------------------------------------------------------------------
