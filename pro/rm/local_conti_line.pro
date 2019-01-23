; local power-law continuum + multi-gaussian lines
; xval in units of alog(wave)
; Used by rm_local_qsofit

function local_conti_line, xval, pp

   npar=n_elements(pp)
   nline=(npar - 2)/3L
   f_line = manygauss(xval, pp[0:3*nline - 1])
   f_pl = pp[3*nline]*(exp(xval)/3000.0)^pp[3*nline+1]  ; power-law continuum

   yval = f_pl + f_line

   return, yval
end

