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

  return
end 


pro callrealspace, width, doplot=doplot

  fname = findfile('quicktest.sav', count=ct)
  if ct eq 0 then testsetup
  restore, 'quicktest.sav'  ;lam, eig, lamg, gal, ivar

  ans = realspace(width, lam, eig, lamg, gal, ivar, testsigma=testsigma, $
                  lamrange=lamrange, doplot=doplot, broadarr=broadarr)

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
;   answers = realspace(width, lam, eig, lamg, gal, ivar, $
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
function realspace, width, lameig, eig, lamg, gal, ivar, testsigma=testsigma, $
                    lamrange=lamrange, doplot=doplot, broadarr=broadarr
  
; Set wavelength range to use for fit
  if n_elements(lamrange) eq 2 then begin 
     lam0 = lamrange[0] & lam1 = lamrange[1]
  endif else begin 
     lam0 = 4000 & lam1 = 5500
  endelse 

; Velocity dispersion is sigma, NOT FWHM.  Measured in pixels. 
  if not keyword_set(testsigma) then testsigma = (findgen(10)+5)/2.

; number of eigenmodes
  neig = (size(eig))[2]

  if not keyword_set(broadarr) then begin 
; number of Fourier components
     nf = 16

; trim wavelength range of lameig
     w   = where((lameig gt lam0) and (lameig lt lam1),nlam)
     lam = lameig[w]            ; good for both eigenmodes and galaxy 
        
; compute smoothed templates once and store in broadarr
     nsigma = n_elements(testsigma)
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
  gal  = double(gal[w])
  ivar = ivar[w]

;  temp = comp_fit(eig, lam, gal, lam1)


  for i=0, nsigma-1 do begin 

     base = broadarr[*, *, i]

     if 1 eq 1 then begin 
        chi2 = computechi2(gal, ivar, base, acoeff=acoeff, dof=dof, yfit=yfit)
        fake = yfit
        C = acoeff
     endif else begin 
        nn = (size(base))[2]
        ep=base/((fltarr(nn)+1)#ivar)
        
        galp=gal/ivar
        
; do fit
        
        G = galp#ep
        
        E = transpose(ep)#ep
        C = invert(E)##G
;  print,C
        
        fake = base#transpose(c)
        
        f2=eig[*,0:6]#transpose(c[*,0:6])
        
        
        chi2 = total(galp^2) - transpose(G)##C
     endelse 
     
     print,'i,sigma,chi2',i, testsigma[i], chi2
     
; Make diagnostic plots
     if keyword_set(doplot) then begin 
        
        plot, lam, gal, /yno
        oplot, lam, fake,color=185
        oplot, lam, (fake-gal)+1, color=215
        
     endif 
     
     print,'stdev', stdev(fake-gal)
     
  endfor 



; assign variables to return
  minchi2 = 0
  minsigma = 0
  errsigma = 0
  bestalpha = 0

  return, [minchi2, minsigma, errsigma, bestalpha]
end


; Legacy


; 12
;chi2      0.32753873
;stdev     0.015254566

; 8
;chi2      0.17706853
;stdev     0.011415969

; 4
;chi2      0.63045199
;stdev     0.021313765













