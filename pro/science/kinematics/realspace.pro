function gconv, x, sigma, edge_wrap=edge_wrap

; special case for no smoothing
  if sigma eq 0 then return, x

  binfactor=1
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
;  print,C

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


pro testsetup
; Read in stellar eigenspectra
  readcol,'sp4',lam,s1,s2,s3,s4
  readcol,'sp321',lam,s5,s6,s7
  readcol,'spbrg',lamg,gal

; Generate fake inverse variance
  nlam = n_elements(lamg)
  ivar = fltarr(nlam)+1        ; uniform weighting
  ivar = ivar+(ivar EQ 0)

; Build eigenspectrum array
  eig  = double([[s1],[s2],[s3],[s4],[s5],[s6],[s7]])

  save, lam, eig, lamg, gal, ivar, file='quicktest.sav'


  restore, 'sample_MAIN_sept00.dat'
  i = lindgen(100)
  galflux = galflux[*, i]
  galsig  = galsig[*, i]
  galwave = galwave[*, i]
  sample  = sample[i]


  save, galflux, galsig, galwave, sample, $
    file='small_sample.sav'




  return
end 

pro testrealspace, doplot=doplot
  

  fname = findfile('quicktest.sav', count=ct)
  if ct eq 0 then testsetup
  restore, 'quicktest.sav'  ;lam, eig, lamg, gal, ivar
  restore, 'small_sample.sav'

  N = 200
  b = fltarr(N)
  berr = fltarr(N)
  chi2 = fltarr(N)

  testsigma = (findgen(20)+0)/2.

  galivar = 1./((galsig^2) + (galsig eq 0))
  fac = (10.^(1e-4)-1)*299792   ; 69.0458

  i = 0
  
  gal = galflux[*, i]
  ivar = galivar[*, i]
  wavg = galwave[*, i]
  lamg = 10^wavg
  z = sample[i].z
  combine1fiber, wavg, gal, ivar, newloglam=wavg+alog10(1+z), $
    newflux=newflux, newivar=newivar
  iseed = !pi

  FOR i=0, N-1 DO BEGIN 
     
     nfac = i/50. ; noise factor
     noise = randomn(iseed, n_elements(newflux))*nfac
     nivar = 1./(1./ivar + nfac^2)
     
     ans = realspace(lam, eig, lamg, newflux+noise, nivar, $
                     testsigma=testsigma, $
                     lamrange=lamrange, doplot=doplot, broadlam=broadlam, $
                     broadarr=broadarr, dof=dof)

     minchi2  = ans[0]
     minsigma = ans[1]
     errsig   = ans[2]
;     dof = n_elements(gal)
     b[i] = minsigma*fac
     berr[i] = errsig*fac
     chi2[i] = minchi2/dof

     print, i, ' Chi2/dof', minchi2, '/', long(dof),  $
       '  sigma:', minsigma, ' pix', $
       minsigma*fac, ' +-', errsig*fac, ' km/s', $
       format='(I5,A,F10.2,A,I5,A,F6.3,A,F8.1,A,F6.1,A)'

  ENDFOR 

  save, b, berr, chi2, file='test.sav'

  plot, b, ps=7
  errplot, b-berr, b+berr


  return
end




pro callrealspace, doplot=doplot
  

  fname = findfile('quicktest.sav', count=ct)
  if ct eq 0 then testsetup
  restore, 'quicktest.sav'  ;lam, eig, lamg, gal, ivar
;  restore, 'small_sample.sav'
  restore, 'sample_MAIN_sept00.dat'

  N = (size(galflux))[2]
  b = fltarr(N)
  berr = fltarr(N)
  chi2 = fltarr(N)

  testsigma = (findgen(20)+0)/2.

  galivar = 1./((galsig^2) + (galsig eq 0))
  fac = (10.^(1e-4)-1)*299792   ; 69.0458

  FOR i=0, N-1 DO BEGIN 
     gal = galflux[*, i]
     ivar = galivar[*, i]
     wavg = galwave[*, i]
     lamg = 10^wavg
     z = sample[i].z
     combine1fiber, wavg, gal, ivar, newloglam=wavg+alog10(1+z), $
       newflux=newflux, newivar=newivar
     
     
     ans = realspace(lam, eig, lamg, newflux, newivar, $
                     testsigma=testsigma, $
                     lamrange=lamrange, doplot=doplot, broadlam=broadlam, $
                     broadarr=broadarr, dof=dof)

     minchi2  = ans[0]
     minsigma = ans[1]
     errsig   = ans[2]
;     dof = n_elements(gal)
     b[i] = minsigma*fac
     berr[i] = errsig*fac
     chi2[i] = minchi2/dof

     print, i, 'Chi2/dof', minchi2, '/', dof, ' sigma:', minsigma, ' pix ', $
       minsigma*fac, ' +-', errsig*fac, ' km/s'
  ENDFOR 

  save, b, berr, chi2, file='results.sav'

  return
end





;------------------------------------------------------------------------------
;+
; NAME:
;   realspace
;
; PURPOSE:
;   Perform a fit of broadened PCA templates to a galaxy spectrum
;   in order to measure velocity dispersion.  The PCA templates are
;   derived from stars which are assumed to have zero dispersion and
;   identical redshifts.  The fit is done in real (not Fourier) space
;   so that the uncertainties in the galaxy spectrum can be used. 
;
; CALLING SEQUENCE:
;   answers = realspace(lam, eig, lamg, gal, ivar, $
;               testsigma=, lamrange=, /doplot, broadarr= )
;
; INPUTS:
;   lam        - wavelengths for eigenmodes (eig) - log spacing
;   eig        - array [nwav, nspec] of eigenspectra (derived from
;                stars)
;   lamg       - wavelengths for galaxy 
;   gal        - galaxy spectrum
;   ivar       - galaxy inverse variance
;
; OPTIONAL KEYWORDS:
;   testsigma  - Array of sigma values to try
;   lamrange   - wavelength range (Angstroms) to use in fit
;   doplot     - Output diagnostic plots to Xwindow
;   broadarr   - array of pre-broadened templates (calculated if not passed)
;   broadlam   - lambda array (A) for broadarr
;
; OUTPUTS:
;   answers    - Four element array with:
;                [minchi2, minsigma, errsigma, bestalpha]
;                bestalpha is the normalization constant between
;                galaxy and star
;
; OPTIONAL OUTPUTS:
;   broadarr   - if undefined when passed, broadarr will be calculated
;                and returned for use in subsequent calls
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Sep-2000  Written by Doug Finkbeiner, UC Berkeley 
;                     with Daniel Eisenstein and David Schlegel. 
;-
;------------------------------------------------------------------------------
function realspace, lameig, eig, lamg, gal_in, ivar_in, $
                    testsigma=testsigma, broadlam=broadlam, $
                    lamrange=lamrange, doplot=doplot, broadarr=broadarr, $
                    dof=dof
  
; Set wavelength range to use for fit
  if n_elements(lamrange) eq 2 then begin 
     lam0 = lamrange[0] & lam1 = lamrange[1]
  endif else begin 
     lam0 = 4000 & lam1 = 5500
  endelse 

; Velocity dispersion is sigma, NOT FWHM.  Measured in pixels. 
  if not keyword_set(testsigma) then testsigma = (findgen(10)+5)/2.

  nsigma = n_elements(testsigma)
  if not keyword_set(broadarr) then begin 
; number of eigenmodes
     neig = (size(eig))[2]

; number of Fourier components
     nf = 16

; trim wavelength range of lameig
     w   = where((lameig gt lam0) and (lameig lt lam1),nlam)
     lam = lameig[w]            ; good for both eigenmodes and galaxy 
     broadlam = lam

; compute smoothed templates once and store in broadarr
     broadarr = dblarr(nlam, neig+nf, nsigma)

     u = (findgen(nlam)+0.5)/nlam
     e0 = 1  ; could be a spectrum

     for i=0, nsigma-1 do begin 
        base = eig
        for j=0, neig-1 do base[*,j] = gconv(eig[*,j], testsigma[i])
        base = base[w, *]
        ; add Fourier components 
        for k=1,nf/2 do base=[[base],[sin(u*!pi*k)*e0],[cos(u*!pi*k)*e0]]
        broadarr[*, *, i] = base
     endfor 

  endif 

  w    = where((lamg gt lam0) and (lamg lt lam1),nlam)
  gal  = double(gal_in[w])
  ivar = ivar_in[w]
  dof  = n_elements(gal)

;  temp = comp_fit(eig, lam, gal, lam1)

  chi2 = fltarr(nsigma)

  for i=0, nsigma-1 do begin 

     base = broadarr[*, *, i]
     chi2[i] = computechi2(gal, sqrt(ivar), base, acoeff=acoeff, dof=dof, yfit=yfit)
     fake = yfit
     C = acoeff
     
;     print,'i,sigma,chi2',i, testsigma[i], chi2[i]
     

     if keyword_set(doplot) then begin ; Make diagnostic plots
        plot,  broadlam, gal, /yno
        oplot, broadlam, fake,color=185, thick=3
        oplot, broadlam, (fake-gal)+1, color=215
     endif 
     
  endfor 

  findchi2min, testsigma, chi2, minchi2, minsigma, errsigma


; assign variables to return
  bestalpha = 0  ; deprecated

  return, [minchi2, minsigma, errsigma, bestalpha]
end








