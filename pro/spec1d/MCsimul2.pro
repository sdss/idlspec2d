
; code to generate Monte Carlo simulations of velocity dispersions
;  as a function of template, S/N, and sigma

Pro MCsimul2, filename

x=findgen(3918)

readidlout,starflux,sig=starsig,wave=starwave

for i=0,31 do begin 
  res=80.0/70.0 
  y=starflux[*,i] 
  simgalflux=fltarr(3918,10,9,16)

  for n=0,9 do begin 
    res=res+20.0/70.0 
    ybroad=gauss_smooth(x,y,res,x) 
    sn=0.0
    print, i, n

    for nn=0,8 do begin 
      sn=sn+10.0 
      meanstar=mean(ybroad) 

      for nnn=0,15 do begin 
        ynoise=ybroad + randomu(100*nnn,3918,/normal)*meanstar/sn
; veldisp, ynoise, starsig, starflux, starsig, results & $
; simsigma[,i,n,nn,nnn] ???
        simgalflux[*,n,nn,nnn]=ynoise
     endfor 
    endfor 
  endfor 

  mwrfits, simgalflux, filename
endfor

return
end
