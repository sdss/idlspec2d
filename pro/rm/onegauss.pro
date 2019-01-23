function onegauss, xval, pp, dp

   term1 = exp( - (xval - pp[1])^2 / (2.D * pp[2]^2))
   yval = pp[0] * term1 / (sqrt(2.D*!pi) * pp[2])
   if (arg_present(dp)) then begin
      message, 'Derivatives not supported yet!!!???'
; These derivatives are for a different parameterization...
;      dp0 = term1
;      dp1 = yval * ((xval - pp[1]) / pp[2]^2)
;      dp2 = yval * ((xval - pp[1])^2 / pp[2]^3)
;      dp = [dp0,dp1,dp2]
   endif

   return, yval
end

