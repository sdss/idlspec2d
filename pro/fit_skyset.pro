; Inputs: xarc, yarc, arclambda, xsky, ysky, skylambda, skycoeff
; Outputs: goodlines, wset
; optional input keywords ymin, ymax, func

pro fit_skyset, xarc, yarc, arclambda, xsky, ysky, skylambda, skycoeff, $
 goodlines, wset, ymin=ymin, ymax=ymax, func=func

   if (keyword_set(func) EQ 0) then func = 'legendre'
   if (keyword_set(ymin) EQ 0) then ymin = wset.xmin
   if (keyword_set(ymax) EQ 0) then ymax = wset.xmax

   narctrace=(size(xarc))[1]
   nskytrace=(size(xsky))[1]
   if (nskytrace NE narctrace) then $
    message, 'Different number of traces found in arcs and skies'
   ntrace = narctrace
   narclines = (size(xarc))[2]
   nskylines = (size(xsky))[2]

   ncoeff = (size(wset.coeff))[1]
   npix=2048
 
    ymid = 0.5*(wset.xmax + wset.xmin)
    yrange = wset.xmax - wset.xmin
    xx = dindgen(2048)
    pixarray = 2.0d0*dindgen(npix)/(npix-1) - 1.0d0

    if (func EQ 'legendre')  then f_pixarray = flegendre(pixarray,ncoeff)
    if (func EQ 'chebyshev') then f_pixarray = fchebyshev(pixarray,ncoeff)

    if (func EQ 'legendre')  then function_name = 'flegendre'
    if (func EQ 'chebyshev') then function_name = 'fchebyshev'

    x = [[xarc], [xsky]]
    lambda = [arclambda, skylambda]

;wavnorm = 2.0d*xnew/(npix-1) - 1.0d
   wavnorm = 2.0d*x/(ymax-ymin) - 1.0d

   nline = n_elements(lambda)

   goodlines = lonarr(nline,ntrace) + 1

   for i=0, ntrace-1 do begin

     done = 0

     while (done EQ 0) do begin
       
       use = where(goodlines[*,i] NE 0)
       ncurarc = long(total(goodlines[0:narclines-1,i]))

       res = arcsky_fit((wavnorm[i,use])[*], lambda[use], ncurarc, $
               ncoeff, skycoeff, function_name=function_name, yfit=yfit)
       diff = yfit - lambda[use]


;
;   Take lines within 20 km/s - throw out 1 at a time
;
       bad = where(abs(diff) GT 3.0d-5, nbad)
       if (nbad EQ 0) then done = 1 $ 
       else begin
         maxdiff = max(abs(diff),badplace)
         goodlines[use[badplace],i] = 0
       endelse 
     endwhile

      wset.coeff[*,i] = res[0:ncoeff-1]
      wset.coeff[0:skycoeff-1,i] = $
      wset.coeff[0:skycoeff-1,i] + res[ncoeff:*]

      newfit = flegendre((wavnorm[i,use])[*], ncoeff) # wset.coeff[*,i]
       
;      for j=ncurarc,n_elements(use) - 1 do $
;         splog, 'Line ', j, yfit[j], lambda[use[j]], diff[j], newfit[j]

      print, format='($, ".",i4.4,a5)',i,string([8b,8b,8b,8b,8b])
   endfor   

   return
end
