;
; NAME:
;   template
;
; PURPOSE:
;   Given N template stars the code generates NT templates
;   by combining random sampling of the N template stars
;
; CALLING SEQUENCE:
;  template, starflux, startemp, ntemplate = ntemplate
;
; INPUTS:
;   starflux   - Template spectra [npix, N]
;
;
; OPTIONAL KEYWORDS:
;   ntemplate  - the number (NT) of generated templates
;
; OUTPUTS:
;   startemp   - The combined template spectra [npix, NT]
;
;   starfact   - the random numbers (1:10) which multiply the
;                template stars in order to generate
;                the combined templates.
;                It is an array [N,NT] 
;
; REVISION HISTORY:
;
;    2-Aug-2000  written by M. Bernardi
;------------------------------------------------------------------------------
Pro template, starflux, startemp, starfact, ntemplate = ntemplate

   npix=(size(starflux))[1]
   nstar= (size(starflux))[2]
 
   if (NOT keyword_set(ntemplate)) then ntemplate = 1

   startemp=fltarr(npix,ntemplate)
   starfact=intarr(nstar,ntemplate)

   for i=0,ntemplate-1 do begin 

      nrand=long(randomu(105*i,nstar,/uniform)*10)
      test=nrand[0]*starflux[*,0]

      for k=1,nstar-1 do begin
          test = test+nrand[k]*starflux[*,k]
      endfor

      startemp[*,i] = test
      starfact[*,i] = nrand 
   endfor 

return
end