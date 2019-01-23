; calculate the RMS spectrum using the ML estimator (Keith)
; Note the current ML implementation by Jon Trump uses the median as the 
; average, instead of using the optimal mean as the average of the LC
; Input
;   RMID
;   epoch: [1,2,3, ...] -- array of epochs

pro rm_calc_rms_spec, rmid, epoch=epoch, result=result,_extra=extra, calibdir=calibdir

if n_elements(result) ne 0 then tmp=temporary(result)
if ~keyword_set(calibdir) then calibdir='wh_skysub/'

if n_elements(rmid) gt 1 then begin
  for jj=0, n_elements(rmid) - 1 do begin
    rm_calc_rms_spec, rmid[jj], epoch=epoch, result=result1,_extra=extra, calibdir=calibdir
    if n_elements(result) eq 0 then result=result1 else result=[result,result1]
  endfor
  return
endif

target_file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
tt = mrdfits(target_file,1,/silent)


plate=tt[rmid].plate & fiber=tt[rmid].fiberid & mjd=tt[rmid].mjd
if keyword_set(epoch) then begin
   plate=plate[epoch-1]
   fiber=fiber[epoch-1]
   mjd=mjd[epoch-1]
endif

rm_readspec, plate, fiber, mjd=mjd, wave=wave, flux=flux, flerr=err, calibdir=calibdir, _extra=extra

npix = (size(flux))[1]
nep = (size(flux))[2]

wave = wave[*,0]
result={wave:wave, medflux:dblarr(npix),rmsflux:dblarr(npix), rmserr:dblarr(npix)-1.D, nspec:lonarr(npix)}
medflux=dblarr(npix) & rmsflux=dblarr(npix) & rmserr=dblarr(npix)-1.D & nspec=lonarr(npix)

for i=0, npix - 1 do begin

   flux1 = flux[i, *] & err1 = err[i,*] 
  
   ind=where(err1 gt 0, ngood)
   nspec[i] = ngood

   if ngood gt 2 then begin
     flux1 = flux1[ind] & err1 = err1[ind]
     medflux[i] = median(flux1)
     rms = getrms(flux1, err1)
     rmsflux[i] = rms[0]
     rmserr[i] = rms[1]
   endif
endfor

result.medflux = medflux
result.rmsflux = rmsflux
result.rmserr = rmserr
result.nspec = nspec

end
