; use the recipes in Shen++(2016, the velocity shift paper) to get zsys and stat+sys uncertainties


function get_zsys, obs_lam, obs_lam_err, line_use=line_use, logL=logL

   cs = 2.9979246d5
 
   line=['Hbeta_br', 'OIII5007', 'CaII3934', 'OII3728', 'NeV3426', 'MgII', 'MgII_br', 'CIII', 'HeII1640', 'CIV_br', 'SIIV_OIV']
   line=strupcase(line)
   wave=[4862.68, 5008.24, 3934.78, 3728.48, 3426.84, 2798.75, 2798.75, 1908.73, 1640.42, 1549.06, (1396.76 + 1402.06)*0.5]

   ; average offset from CaII, if applicable, after luminosity-trend correction
   avgoff = [-109., -48., 0., 8., -160., -57., -57., -229., 0., -27., -28.]
   ; average offste from CaII, before luminosity-trend correction
   avgoff2 = [-109., -48., 0., 8., -160., -57., -57., -229., -167., -365., -288.]

   ; average luminosity trend; only HeII1640, CIV_br and SIIV_OIV requires a luminosity correction
   ; and logL must be logL1700
   logL0 = [44., 44., 0., 44.5, 44.5, 44.5, 45., 45., 45., 45.]
   arr = [-88, -36, 0., 9., -121., -75., -75., -272., -231., -242., -123.]
   barr = [-4., -53., 0., -0.04, -2., 2., 2., -0.2, -282., -438., -345.]
   ; Eqn 1 in Shen++: v = a + b(logL - logL0)

   ; intrinsic velocity scatter wrt CaII in km/s (Table 3 of Shen++)
   sigv_int = [400., 56., 0., 46., 119., 205., 205., 233., 242., 415., 477.]


   ind = where(line eq strupcase(line_use) )
   if ind[0] ne -1 then begin

       if line[ind] eq 'HEII1640' or line[ind] eq 'CIV_BR' or line[ind] eq 'SIIV_OIV' then begin
          if logL gt 0 then v_add = arr[ind] + barr[ind]*(logL - logL0[ind]) $
            else v_add = avgoff2[ind] ; if luminosity is unspecified, use the average vel offset
       endif else v_add = 0.

       ; this is the total average voff from systemic (CaII)
       v_tot = avgoff[ind] + v_add
       rest_lam = wave[ind] * (1.D + v_tot/cs)

       ; now compuate zsys based on the specific line
       zsys = (obs_lam / rest_lam) - 1.D
        
       zsys_err = sqrt( (sigv_int[ind]/cs*(1. + zsys))^2 + ( (obs_lam_err > 0.)/rest_lam)^2  )

   endif else begin
      zsys = -1. & zsys_err = -1.

   endelse


   return, [zsys, zsys_err]

end
