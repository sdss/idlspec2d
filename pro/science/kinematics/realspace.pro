function gconv, x, fwhm, edge_wrap=edge_wrap

  binfactor=1
  sigma=fwhm/(2.355*binfactor)
  ksize=round(4*sigma+1)*2
  xx = findgen(ksize)-ksize/2

  kernel=exp(-xx^2/(2*sigma^2))
  kernel=kernel/total(kernel)

  sm = convol(x, kernel, edge_wrap=edge_wrap)

  return, sm
end 



function matfit, eig, gal, sigma, inv=inv

  neig=(size(eig))[2]
  ep=eig/((fltarr(neig)+1)#sigma)

  galp=gal/sigma

  G = galp#ep

  E = transpose(ep)#ep
  inv = invert(E)
  C = inv##G
  print,C

return, c
end


function comp_fit, eig, lam, gal, lam1, smooth_sig=smooth_sig

; Generate Fourier components
  nlam=n_elements(lam)
  u = (findgen(nlam)+0.5)/nlam
  neig=(size(eig))[2]

  base = eig
  if keyword_set(smoothsig) then begin 
    for i=0,neig-1 do base[*,i]=gconv(eig[*,i], smooth_sig*2.355)
  endif
  e0=1
  for k=1,5 do base=[[base],[sin(u*!pi*k)*e0],[cos(u*!pi*k)*e0]]

  C = matfit(base, gal, fltarr(nlam)+1, inv=inv)

  fit = eig#C[0:neig-1]  
;  fit = base#transpose(C)

  return, fit
end


pro doug, width

lam0=4000
lam1=5500

readcol,'sp4',lam,s1,s2,s3,s4
readcol,'sp321',lam,s5,s6,s7
readcol,'spbrg',lamg,gal

nlam=n_elements(lamg)
sigma = fltarr(nlam)+1  ; uniform weighting
sigma=sigma+(sigma EQ 0)


eig=double([[s1],[s2],[s3],[s4],[s5],[s6],[s7]])
neig=(size(eig))[2]


for i=0,neig-1 do eig[*,i]=gconv(eig[*,i], width)



w=where((lam gt lam0) and (lam lt lam1),nlam)
eig=eig[w,*]
lam=lam[w]

w=where((lamg gt lam0) and (lamg lt lam1),nlam)
gal=double(gal[w])
sigma=sigma[w]

temp = comp_fit(eig, lam, gal, lam1)

u = (findgen(nlam)+0.5)/nlam
k=1
e0=temp
base=eig
for k=1,8 do base=[[base],[sin(u*!pi*k)*e0],[cos(u*!pi*k)*e0]]
neig=(size(base))[2]



ep=base/((fltarr(neig)+1)#sigma)

galp=gal/sigma

G = galp#ep

E = transpose(ep)#ep
C = invert(E)##G
print,C

fake = base#transpose(c)

f2=eig[*,0:6]#transpose(c[*,0:6])


chi2 = total(galp^2) - transpose(G)##C
print,'chi2',chi2

plot, lam, gal, /yno
oplot, lam, fake,color=155
oplot, lam, (fake-gal)+1, color=185

print,'stdev', stdev(fake-gal)

return
end
















