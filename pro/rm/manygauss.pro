function manygauss, xval, pp, nline=nline

  if n_elements(nline) eq 0 then nline = n_elements(pp)/3L

   yval = 0.d
   for iline=0, nline-1 do $
    yval = yval + onegauss(xval, pp[iline*3:iline*3+2])

   return, yval
end

