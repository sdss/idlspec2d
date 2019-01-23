function f_poly_conti, xval, pp

  xval2 = xval - 3000.
  yval = 0*xval2
  for i=0L, n_elements(pp) - 1L do $
    yval = yval + pp[i]*xval2^(i+1)

  return, yval
end

