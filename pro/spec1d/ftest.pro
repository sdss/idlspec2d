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




PRO ftest

  galfwhm = 3.5
  templatefwhm = 2.5

  npix = 3918
  iseed = !pi
  raw = randomn(iseed, npix)
  ker = gausskernal(templatefwhm)
  smraw = convol(raw, ker, /edge_wrap)

  template = smraw+randomn(iseed, npix)*0.1

  ker = gausskernal(galfwhm)
  gal = shift(convol(raw, ker, /edge_wrap), 3)
  gal = gal+randomn(iseed, npix)*0.2
 
  templatesig = template*0.+1
  galsig = gal*0.+1

  veldisp, gal, galsig, template, templatesig, result2, /doplot, sigmast=0.025
  truedisp = sqrt(abs(galfwhm^2-templatefwhm^2))
  print, 'True vel disp = ', truedisp/sqrt(8.*alog(2))


  return

END

