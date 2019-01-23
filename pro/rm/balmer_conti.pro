; Balmer continuum from the model of Dietrich+02

function balmer_conti, xval, pp
; xval = input wavelength, in units of A
; pp=[norm, Te, tau_BE] -- in units of [--, K, --]

lambda_BE = 3646.D  ; A

bbflux = planck(xval,pp[1])  ; in units of ergs/cm2/s/A

tau = pp[2]*(xval/lambda_BE)^3
result = pp[0]*bbflux*(1. - exp(-tau))

ind = where(xval gt lambda_BE)
if ind[0] ne -1 then result[ind]=0.D

return, result
end

