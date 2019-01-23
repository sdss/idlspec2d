;*****************************************************************************
; NAME:
; varlc
;
; AUTHOR:
; Jonathan Trump, UConn 2016
;
; PURPOSE:
; Calculate the intrinsic variance of a lightcurve, and the error in
; that intrinsic variance.  Reports results from three methods:
;   maximum likelihood (ML) - adapted from Keith Horne notes & fortran
;                             program
;   interquartile range (IQR) - following MacLeod et al. 2010
;   absolute average deviation (AAD) - e.g., vanden Berk et al. 2004
;
; INPUTS:
; Lightcurve files within directory "mergedlc_<filter>", where
; <filter> is hardcoded in the variable filter. Files must have first
; 3 columns: mjd/flux/fluxerr.
;
; OUTPUT:
; File named "varlc_<filter>.dat" with columns:
;   ID  iflux med(err)  rms_ml err   rms_iqr err   rms_aad err
;
;****************************************************************************


function getrms, fluxes, errors

;; Keith Horne's method to estimate intrinsic variability.
;;
;; V = sum( (xi-u)^2 gi^2 ) / sum(gi)
;;   V = intrinsic variance
;;   xi = flux measurement
;;   u = mean flux
;;   gi = 1 / (1+e^2/V)
;;   e = flux error
;;
;; Iterate to solve for V.  Note that e should be a combination of
;; the standard flux error plus the spectrophotometry error.

  if stddev(fluxes) eq 0 or n_elements(fluxes) le 1 then return,[0,-1]

  meanflux = median(fluxes)
;;   errors = sqrt(errors^2 + (0.04*meanflux)^2)  ;add 0.04% specphoterror

  varguess = (stddev(fluxes)^2 - median(errors)^2)>0.05

  repeat begin
     varold = varguess
     gi = 1 / (1+errors^2/varold)
     varnew = total( (fluxes-meanflux)^2 * gi^2 ) / total(gi)
     varguess = varnew
  endrep until abs(varold-varnew)/varnew le 1e-5 or varnew lt 1e-12

  varerror = 2d / ( 2*total( (fluxes-meanflux)^2 / (varnew+errors^2)^3 ) $
                            - total(1d/(varnew+errors^2)^2) )

  return, [sqrt(varnew>0), sqrt(varerror>0)/sqrt(varnew>1e-4)]

;  varerror2 = 2d * ( 2d*total((fluxes-meanflux)^2 * gi^3) $
;                    / total((fluxes-meanflux)^2 * gi^2) * total(gi) - total(gi^2) )
;
;  return, [sqrt(varnew>0), sqrt(varnew>0)/sqrt(varerror2>1e-4)]

end


function iqrerror, flux,error

  niter = 1000
  nv = n_elements(flux)+1
  iqr = fltarr(niter)

  for ee=0, niter-1 do begin
     ind = long(randomu(seed,nv+1)*(nv+1))<nv
     flux1 = flux[ind] + randomn(seed,nv+1)*error[ind]
     fluxsort = sort(flux1)

     iqr[ee] = flux1[fluxsort[3*nv/4]] - flux1[fluxsort[nv/4]]
     endfor

  return,stddev(0.74*iqr)
end


function getaad, fluxes, errors

  nn = n_elements(fluxes)

  aad = sqrt( !pi/2 * (avg(abs(fluxes-avg(fluxes)))^2 - avg(errors^2))>0 )

  aaderr = (aad eq 0) ? 0 : 0.5/aad * sqrt( !pi^2*avg(abs(fluxes-avg(fluxes)))^2 $
                                       * stddev(fluxes)^2/nn + stddev(errors^2)^2/nn)

  return,[aad,aaderr]

end


pro varlc

filter='mg2'
dir = 'mergedlc_'+filter+'/'
nqso = 850

files = findfile(dir+'rm*'+filter+'_cream.dat')
sdssdata = mrdfits('../target_fibermap.fits',1)
iflux = 10^(-0.4*(sdssdata[0:nqso-1].psfmag[3]+48.6)+29) ;uJy units

openw,1, 'varlc_'+filter+'.dat'    ;'_linefit.dat'
printf,1,'# ID  iflux med(err)  rms_ml err   rms_iqr err   rms_aad err'

for ii=0, n_elements(files)-1 do begin
   qq = long(strmid(files[ii],strpos(files[ii],'/')+3,3))
   readcol,files[ii],mjd,flux,err, /silent

   iqr=0.0
   iter=0
   maxiter=10
   repeat begin

      toterr = sqrt(err^2+iqr^2)
      linefit = linfit((mjd-6600),flux, measure_errors=toterr)
      flux0 = flux - (linefit[0] + linefit[1]*(mjd-6600))
   
 ;ML variance estimator
      var = getrms(flux0,err)

 ;IQR: range of 25-75% of cumulative distribution,
 ;     with correction (x0.74) to equal sigma for a Gaussian
      vsort = sort(flux0)
      nv = n_elements(flux0)-1
      iqr_obs = 0.74 * (flux0[vsort[3*nv/4]] - flux0[vsort[nv/4]])
      iqrerr_obs = iqrerror(flux0,err)
      iqr = sqrt(iqr_obs^2 - median(err)^2) > 0
      iqrerr = (iqr ne 0) ? 1/iqr * sqrt( (iqr_obs*iqrerr_obs)^2 $
                                          + medabsdev(err^2)^2/4/nv) : 0

 ;AAD: absolute average deviation (e.g. vanden Berk et al 2004)
      aad = getaad(flux0,err)

      iter++
   endrep until median(abs(toterr^2 - (err^2+iqr^2))/flux) le 1e-4 or iter gt maxiter

   printf,1,qq,iflux[qq], median(err), var, iqr,iqrerr, aad, $
          format='(I3,"  ",F6.2,7("  ",F6.2))'

   if ii mod 10 eq 0 then print,strtrim(ii,2),' done!'
endfor

close,1

end
