;+
; NAME: 
;   general_sn
; 
; PURPOSE:
;    Computes the S/N versus flux (10^(0.4*(22.5-mag))) 
;    It is used by fitsn_jb() 
;    The model assumes signal = flux and 
;    noise = flux + sky (with the sky value as constant)
; 
; CALLING_SEQUENCE:
;    general_sn( x, pars)
; 
; INPUTS:
;    x - flux in nanomaggies = 10^(0.4*(22.5-mag))
;    pars - two parameters of model
;
; OUTPUTS:
;    signal to noise
;
; REVISION HISTORY:
;    14-Jun-2015 Written by J. Bautista
;- 
;----------------------------------------------------------------
function general_sn, x, pars
	a = pars[0]
	b = pars[1]
	return, a*x/sqrt(x+b)
end
