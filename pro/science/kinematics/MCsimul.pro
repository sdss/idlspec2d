;+
; NAME:
;   MCsimul
;
; PURPOSE:
;   For each template star the code generates Monte Carlo simulations 
;   of galaxy-spectra as function of S/N (10 to 90), and
;   velocity dispersion (100 to 280 km/s). At each fixed S/N and 
;   velocity dispersion it generates 16 random realizations.
;
; CALLING SEQUENCE:
;   MCsimul, starflux, filename
;
; INPUTS:
;   starflux   - Template spectra [npix, nstar]
;
; OUTPUTS:
;   filename   - a .fit file containing the output 
;                array [npix,nveldisp,nSN,nrand] (npix, 10, 9, 16) 
;                with the simulated spectra
;
; REVISION HISTORY:
;
;   16-Mar-2000  Written by M. Bernardi & S. Burles, FNAL
;    2-Aug-2000  Modified by M. Bernardi
;------------------------------------------------------------------------------
Pro MCsimul, starflux, filename

   npix=(size(starflux))[1]
   nstar= (size(starflux))[2]
   x=findgen(npix)

   for i=0,nstar-1 do begin 

       res=80.0/70.0 
       y=starflux[*,i] 
       simgalflux=fltarr(npix,10,9,16)

       for n=0,9 do begin 
            res=res+20.0/70.0 
            ybroad=gauss_smooth(x,y,res,x) 
            sn=0.0
            print, i, n

            for nn=0,8 do begin 
               sn=sn+10.0 
               meanstar=mean(ybroad) 

                   for nnn=0,15 do begin 
                     ynoise=ybroad + randomu(100*nnn,npix,/normal)*meanstar/sn
                     simgalflux[*,n,nn,nnn]=ynoise
                   endfor 
            endfor 
        endfor 

     mwrfits, simgalflux, filename

   endfor

return
end
