pro mfindchi2min, x, chi2, minchi2, minsigma, errsigma, npts=npts, $
          deltachisq = deltachisq, nfine = nfine, chi2scale = chi2scale,$
         doplot=doplot

   if (NOT keyword_set(deltachisq)) then deltachisq = 1.0
   if (NOT keyword_set(nfine)) then nfine = 100

   ; spline chi2 curve and return desired values

;   if (keyword_set(npts)) then chi2=chi2/npts

   nx = n_elements(x)
   if (nx EQ 0) then return

   minx = min(x, max=maxx)
   y2 = spl_init(x, chi2)

   x2 = findgen(nx*nfine+1)/(nx*nfine) * (maxx - minx) + minx
   chifit = spl_interp(x, chi2, y2, x2)

   minchi2 = min(chifit,fitplace)

   minsigma = x2[fitplace] 

   ;
   ; attempt to scale chi^2 to get errors correct
   ; 
  
   usethisdelta = 2.3

;   if (keyword_set()) then usethisdelta = 2.3 * minchi2 * deltachisq


;   if (keyword_set(npts)) then usethisdelta = sqrt(minchi2)

   if (keyword_set(npts)) then begin

;   minchi2point=min(chi2,pos)
;   deltasigma= 0.2*(x[pos]-x[pos-1])

;   if (x[pos] GT minsigma+deltasigma) then begin
;   chi2right=chi2[pos]
;   chi2left=chi2[pos-1]
;   endif 

;   if (x[pos] LT minsigma-deltasigma) then begin 
;   chi2right=chi2[pos+1]
;   chi2left=chi2[pos]
;   endif

;   if (x[pos] GE minsigma-deltasigma AND x[pos] LE minsigma+deltasigma) then begin 
;   chi2right=chi2[pos+1]
;   chi2left=chi2[pos-1]
;   endif

;   if (chi2right LT chi2left) then diffchi2= chi2left-minchi2 $
;   else diffchi2= chi2right-minchi2


   if (keyword_set(chi2scale)) then begin
   if (chi2scale * deltachisq GT 1) then $
   usethisdelta =  2.3 * chi2scale * deltachisq 
   endif else begin 
   if (minchi2 * deltachisq GT 1) then $
   usethisdelta = 2.3 * minchi2 * deltachisq 
   endelse

   endif	

;   if (keyword_set(npts)) then begin
;     scale =  minchi2/npts     		
;     if (scale GT 1.0) then usethisdelta = deltachisq * scale
;   print,'scale',scale,usethisdelta
;   endif

   range = where(chifit - minchi2 LT usethisdelta)

   if (range[0] EQ -1) then return

   uppersigma = x2[max(range)] - minsigma
   lowersigma = minsigma - x2[min(range)] 

   errsigma = 0.5*(uppersigma + lowersigma)


   if (keyword_set(doplot)) then begin
     if (doplot EQ 1) then begin
       wset,1
       plot, x, chi2-minchi2, ps=1, yr=[-10,100], $
          title='Chi2: diff (crosses), realspace (diamonds)'
     endif
     if (doplot EQ 2) then djs_oplot, x, chi2-minchi2, ps=4    
     djs_oplot, x2, chifit-minchi2
     djs_oplot, x2[range], (chifit-minchi2)[range],color='red'


   endif

;    if (keyword_set(npts)) then begin
;    minchi2= minchi2*npts
;    stop
;   endif	
   return
end    
   



