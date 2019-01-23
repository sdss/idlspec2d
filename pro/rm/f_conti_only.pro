; return the combination of three continuum components:
; PL, Balmer continuum and polynominal.
; Called by rm_qsofit

function f_conti_only, xval, pp

   f_pl = pp[0]*(xval/3000.0)^pp[1]  ; power-law continuum
   f_conti_BC = balmer_conti(xval, pp[2:4])  ; Balmer continuum
   f_poly = f_poly_conti(xval, pp[5:*])
   yval = f_pl + f_conti_BC + f_poly

   return, yval
end

