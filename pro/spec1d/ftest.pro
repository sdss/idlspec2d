PRO ftest_obs

    fftfreq = getidlfreq(npixbell)
    twopiei = 2.0 * !Pi * complex(0.0,1.0)
     phase = exp( - twopiei * fftfreq * fitcen)
  q = dcomplex(fluxfft)/starshift
  fq = fft(q, /inv)
  fq = shift(fq, 10)
  
  a = [2000, 10.0, 0.8]

  res=gaussfit(findgen(8192),double(fq),a,nterms=3)
  print, a

;  3746.25      9.82097     0.812598
;  3742.79      9.83615     0.814107
;
; *phase   3827.23      9.60355     0.754213

  return
END

PRO gauss_periodic_funct, x, a, g, pder

  print, 'gauss_periodic_fn: a, ', a
  amp = a[0]
  cen = a[1]
  sig = a[2]

  n = n_elements(x)

  pder = dblarr(n, 3)
  pd   = dblarr(n, 3)
  g = 0
  FOR i=-3, 3 DO BEGIN 
      z = (x-cen-i*n)/sig
      ez = exp(-z^2/2.d)
      g = g+amp*ez
      pd[*, 0] = ez
      pd[*, 1] = amp*ez*z/sig
      pd[*, 2] = pd[*, 1]*z
      pder = pder+pd
  ENDFOR 

  return
END



PRO foo, sig

  sig = double(sig)
  x = dindgen(8192)-4096
  g = exp(-x^2/(2*sig^2))
  f = shift(fft(g), 4096)
  plot, x, abs(f)

  fsig = 0.5d/!dpi*8192/sig
  print, 'fsig: ', fsig
  est = double([max(abs(f))-min(abs(F)), 100, fsig])

;  res=gaussfit(dindgen(8192),abs(dcomplex(f)),a,nterms=3, estimates=est)
  res = curvefit(x, abs(f), x-x+1.d, est, $
                 FUNCTION_name='gauss_periodic_funct', tol=1.d-19)

;  res = lmfit(x, abs(f), est, FUNCTION_name='gauss_periodiclmfit', /double, itmax=5)

  print, a[2]*sig/8192.*2.d*!dpi

  fg = gauss_periodic( x, [1., 0., fsig])
  oplot, x, fg/(sqrt(2.d*!dpi)*fsig), line=2, color=140

  oplot, x, res, /line

;  plot, x, fg/(sqrt(2.d*!dpi)*fsig)-abs(f), ps=3

  print,stdev(fg/(sqrt(2.d*!dpi)*fsig)-abs(f))/max(abs(f))
  return
END


; pass fwhm in pixels
FUNCTION gausskernal, fwhm

  n = (fix(fwhm*3)*2+1) > 5 ;always odd
  x = findgen(n)-n/2
  sig = fwhm/sqrt(alog(2)*8)
  g = exp(-(x/sig)^2/2.)
  g = g/total(g)

  return, g
END




PRO ftest, result2

  galfwhm = 3.0
  templatefwhm = 2.5

  npix = 3918
  iseed = !pi*2.
  raw = randomn(iseed, npix)
  ker = gausskernal(templatefwhm)
  smraw = convol(raw, ker, /edge_wrap)

  template = smraw;+randomn(iseed, npix)*0.1

  ker = gausskernal(galfwhm)
  gal = shift(convol(raw, ker, /edge_wrap), 3)

  ngal = 250
  gal=gal#(fltarr(ngal)+1)
  FOR i=0, ngal-1 DO gal[*, i]=gal[*, i]+randomn(iseed, npix)*0.157

  templatesig = template*0.+1
  galsig = gal*0.+1

  veldisp, gal, galsig, template, templatesig, result2, $
    sigmast=0.025, /nodiff
  truedisp = sqrt(abs(galfwhm^2-templatefwhm^2))
  truesig2 = (truedisp/sqrt(8.*alog(2)))^2
  print, 'True vel disp (sigma^2) = ', truesig2


  bsig2 = result2.sigma_quotient  ; b^2 (b=vel dispers)
  plothist, bsig2, bin=.010
  oplot, truesig2*[1, 1], [0, 1e6], line=2
  print, mean(bsig2)/truesig2

  return

END



PRO brg, result, z
zap = 1
  listpath = '~/idlspec2d/misc/'
  readcol, listpath+'brg.list', num, plt, fiber, mjd, zchic
;  w = where((plt GE  301) AND (plt LT 303))
  w = where((plt GT 302) AND plt LT  306, ct)
;  w = w[0:9] &  ct=n_elements(w)
;  w = [88] & ct=1
  plt = plt[w]
  mjd = mjd[w]
  fiber = fiber[w]
  zchic = zchic[w]

  readspec, 306, 250, flux=template,wave=wave,flerr=templatesig,plug=plug
  readspec, plt, fiber, flux=galflux,wave=galwave,flerr=galsig,plug=galplug
  bad = bytarr(ct)
  FOR i=0, ct-1 DO BEGIN 
      IF stdev(galflux[*, i]) EQ 0 THEN bad[i] = 1
  ENDFOR 
  w = where(bad EQ 0)
  galflux = galflux[*, w]
  galwave = galwave[*, w]
  galsig  = galsig[*, w]
  galplug = galplug[*, w]
  plt = plt[w]
  mjd = mjd[w]
  fiber = fiber[w]
  zchic = zchic[w]
      

  IF keyword_set(zap) THEN BEGIN 
      
      linemask = (galwave LT 5585) AND (galwave GE 5574)
      galflux = djs_maskinterp(galflux, linemask, iaxis=0)
      
  ENDIF 


  veldisp, galflux, galsig, galwave, template, templatesig, wave, result, $
    sigmast=0.05, maxsig=6, /nodiff

  z = 10.^(result.z/10000.)-1. 
  sig = sqrt(result.sigma_quotient)
  sig2err = (result.sigma_quotient_err)
  mag = galplug.mag[2]
  sigcc = result.sigma_cc

  goodz = where(abs(z-zchic) LT 0.002, ct)
  badz = where(abs(z-zchic) GT 0.002, ct)
  print, 'Bad Fiber numbers: ', fiber[badz]
stop
  return

END


PRO ntest, noisenum, noisedenom

  npix = 100000
  iseed = !pi
  x = double(randomn(iseed, npix))
  y = double(randomn(iseed, npix))
  x = (x-mean(x))/stdev(x)
  y = (y-mean(y))/stdev(y)
  r = sqrt(x^2+y^2)
  plothist, r, bin=.01
  print, mean(r)
  print, stdev(r)
  z = dcomplex(x, y)

  znum = z*noisenum+complex(2, 2)
  zden = z*noisedenom+complex(0, 1)
  q = znum/shift(zden, 1)

  print, mean(abs(q)), !pi/2
  plothist, float(q^0.5),xr=[-1, 3],bin=.01
  plothist, float(q^1.),xr=[-1, 3],bin=.01, /overplot, /line
  print, 'mean(q), median(q): ', mean(float(q)), median(float(q))

  print, 'mean(q^0.5)', mean(q^0.5), stdev(float(q^0.5))
  print, 'exp(mean(alog(q))) ', exp(mean(alog(q))) 

  res = fltarr(51)
  x0 = findgen(51)/10.
;  FOR i=0, 50 DO BEGIN 
;      r = sqrt((x-x0[i])^2+y^2)
;      print, mean(r)
;      res[i] = mean(r)
;  ENDFOR 

  plot, x0, res^2-x0^2

  model=exp(-(x0/1.25)^2/2.)*(!pi/2-1)+1
  plot, x0, res^2-x0^2-model
stop
  return
END

