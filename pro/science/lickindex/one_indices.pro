Pro one_indices, galfluxbest, galwavebest, galsigbest, bestredshift, pp, absline, errabsline, snoise, m = m , q = q, galsmooth = galsmooth

; pp=0  without continuum 
; pp=1  continuum obtained splining only the pseudo continuum regions
; pp=2  continuum obtained splining all the wavelength range
; pp=3  continuum obtained with a straight line

indexl=[4143.375,4143.375,4223.500,4282.625,4370.375,4453.375,4515.500,4635.250,4848.875,4979.000,5070.375,5155.375,5161.375,5247.375,5314.125,5390.250,5698.375,5778.375,5878.625,5938.375,6191.375]
indexr=[4178.375,4178.375,4236.000,4317.625,4421.625,4475.875,4560.500,4721.500,4877.625,5055.250,5135.375,5197.875,5193.875,5287.375,5354.125,5417.750,5722.125,5798.375,5911.125,5995.875,6273.875]
contl=[4081.375,4245.375,4085.125,4245.375,4212.250,4242.250,4267.625,4320.125, 4360.375,4444.125,4447.125,4478.375,4505.500,4561.750,4612.750,4744.000,4828.875,4877.625,4947.750,5055.250,4896.375,5302.375,4896.375,5302.375,5143.875,5192.625,5234.875,5287.375,5306.625,5355.375,5379.000,5417.750,5674.625,5724.625,5767.125,5799.625,5862.375,5923.875,5818.375,6040.375,6068.375,6374.375]
contr=[4118.875,4285.375,4097.625,4285.375,4221.000,4252.250,4283.875,4336.375, 4371.625,4456.625,4455.875,4493.375,4515.500,4580.500,4631.500,4757.750,4848.875,4892.625,4979.000,5066.500,4958.875,5367.375,4958.875,5367.375,5162.625,5207.625,5249.875,5319.875,5317.875,5365.375,5390.250,5427.750,5698.375,5738.375,5777.125,5813.375,5877.375,5949.875,5850.875,6105.375,6143.375,6416.875]

lindl= alog10(indexl)
lindr= alog10(indexr)
indl=alog10(contl)
indr=alog10(contr)

npix=(size(galfluxbest))[1]

x=galwavebest-alog10(1+bestredshift)


if (pp EQ 1 OR pp EQ 2) then begin

galivar = 1.0/galsigbest^2
galivar[where(finite(galivar) EQ 0)] = 0

mask = bytarr(npix)
 
  for k=0,41 do begin & $
    inside = where(x GE indl[k] AND x LE indr[k], ninside) & $
    if (ninside GT 0) then mask[inside] = 1 & $
  endfor 
  galerr = galivar * mask 

if (pp EQ 1) then galer=galerr
if (pp EQ 2) then galer=galivar

fullbkpt = slatec_splinefit(x, galfluxbest, coeff, nbkpt=5, invvar=galer,upper=1.5, lower=0.5)
continuum = slatec_bvalu(x, fullbkpt, coeff)
galsub=galfluxbest/continuum
galnoise = galsigbest/continuum
endif

xv=x
res=8.6/10^(xv)*0.43429 
if (pp EQ 0 OR pp EQ 3) then begin
y=galfluxbest
noise=1/galsigbest^2
endif
if (pp EQ 1 OR pp EQ 2) then begin
y=galsub
noise=1/galnoise^2
endif
inf = where(finite(noise) EQ 0)
if (inf[0] NE -1) then noise[inf]=0
;galsmooth= gauss_smooth_wave(xv,y,res,xv)
galsmooth=y

bcont=fltarr(21)
rcont=bcont
ind=bcont
dw=bcont
m=bcont
q=bcont
xmeanl=bcont
xmeanr=bcont
absline=bcont
errabsline=bcont
snoiseb=bcont
snoiser=bcont
snoise=bcont
acont=bcont
ew=bcont
magline=bcont
errew=bcont
errmagline=bcont
errin=bcont
errcont=bcont
errcontb=bcont
errcontr=bcont
tnb=bcont
tnr=bcont
tni=bcont

k=0
for nn=0,20 do begin 
  bw = where(x GE indl[k] AND x LE indr[k])
  rw = where(x GE indl[k+1] AND x LE indr[k+1])
  iw = where(x GE lindl[nn] AND x LE lindr[nn])

 if (rw[0] NE -1 AND bw[0] NE -1 AND iw[0] NE -1) then begin

    tnb[nn]=total(noise[bw])
    tnr[nn]=total(noise[rw])
    tni[nn]=total(noise[iw])

  if (tnb[nn] EQ 0.0 OR tnr[nn] EQ 0.0 OR tni[nn] EQ 0.0) then begin	
    absline[nn]=-1.0
 
  endif else begin

    bcont[nn]=total(galsmooth[bw]*noise[bw])/tnb[nn]  
    rcont[nn]=total(galsmooth[rw]*noise[rw])/tnr[nn] & $
    ind[nn]=total(galsmooth[iw]*noise[iw])/tni[nn] & $
    errcontb[nn]=sqrt(1/total(noise[bw]))
    errcontr[nn]=sqrt(1/total(noise[rw]))	
    errin[nn]=sqrt(1/total(noise[iw]))
    errcont[nn]=sqrt(errcontb[nn]^2+errcontr[nn]^2)

    dw[nn]=10^(lindr[nn])-10^(lindl[nn]) & $
    snoiseb[nn]= bcont[nn]/mean(1/sqrt(noise[bw]))
    snoiser[nn]= rcont[nn]/mean(1/sqrt(noise[rw]))
    snoise[nn]=	ind[nn]/mean(1/sqrt(noise[iw]))

    if (pp EQ 3) then begin & $
      xmeanl[nn]=10^((indr[k]+indl[k])/2) & $
      xmeanr[nn]=10^((indr[k+1]+indl[k+1])/2) & $
      m[nn]=(rcont[nn]-bcont[nn])/(xmeanr[nn]-xmeanl[nn]) & $
      q[nn]= bcont[nn]- m[nn]*xmeanl[nn] & $
      errm=errcont[nn]/(xmeanr[nn]-xmeanl[nn])
      errq=sqrt(errcontb[nn]^2+(errm*xmeanl[nn])^2)
      err=sqrt((10^(x[iw])*errm)^2+errq^2)
      contin=m[nn]*10^(x[iw])+q[nn]	 
      acont[nn]=total(contin/err^2)/total(1/err^2) 
      errcont[nn]=sqrt(1/total(1/err^2))
      ind[nn]=total(y[iw]/galsigbest[iw]^2)/total(1/galsigbest[iw]^2)
      errin[nn]=sqrt(1/total(1/galsigbest[iw]^2))

;      galnew=y[iw]/(m[nn]*10^(x[iw])+q[nn]) & $
;      galnewerr=galsigbest[iw]/(m[nn]*10^(x[iw])+q[nn])	
;      ind[nn]=total(galnew/galnewerr^2)/total(1/galnewerr^2) & $       
;      errin[nn]=sqrt(1/total(1/galnewerr^2))

      errcontb[nn]=errcontb[nn]/bcont[nn]
      errcontr[nn]=errcontr[nn]/rcont[nn]
      bcont[nn]=1 & $
      rcont[nn]=1 & $
   endif
 
 if (pp NE 3) then begin
acont[nn]=(bcont[nn]+rcont[nn])/2
errcont[nn]=sqrt(errcontb[nn]^2+errcontr[nn]^2)
 endif

; EDITED by hogg
;ew[nn]=dw[nn]*(acont[nn]-ind[nn])/acont[nn]
;errew[nn]= ew[nn]* $
;           sqrt((errcont[nn]^2+errin[nn]^2)/(acont[nn]-ind[nn])^2+ $
;                 errcont[nn]^2/acont[nn]^2)
ew[nn]=dw[nn]*(1.0-(ind[nn]/acont[nn]))
errew[nn]= (dw[nn]*ind[nn]/acont[nn])* $
           sqrt((errin[nn]/ind[nn])^2+(errcont[nn]/acont[nn])^2)

magline[nn]=-2.5*alog10(ind[nn]/acont[nn])
errmagline[nn]=2.5/2.303*sqrt(errin[nn]^2/ind[nn]^2+errcont[nn]^2/acont[nn]^2)


if (nn EQ 0 OR nn EQ 1 OR nn EQ 10 OR nn EQ 11 OR nn EQ 19 OR nn EQ 20) then begin
absline[nn]=magline[nn] 
errabsline[nn]=errmagline[nn]
endif else begin
absline[nn]=ew[nn] 
errabsline[nn]=errew[nn]
endelse

  endelse
 endif

k=k+2 

endfor
return

end