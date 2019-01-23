; 
pro get_all_rms_spec, rmid_arr, zsys_arr, result_all = result_all, nomask=nomask

file = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
target = mrdfits(file,1)

nobj = n_elements(rmid_arr)
plate = target[rmid_arr[0]].plate & fiber = target[rmid_arr[0]].fiberid 
mjd = target[rmid_arr[0]].mjd
ep_str = {plate:plate[0:31], fiberid:fiber[0:31], mjd:mjd[0:31]}
mk_rms_spec, rmid_arr[0], zsys_arr[0], ep_str = ep_str, result=result0, nomask=nomask

result_all = replicate(result0, nobj)
for i=1, nobj - 1 do begin

   plate = target[rmid_arr[i]].plate & fiber = target[rmid_arr[i]].fiberid 
   mjd = target[rmid_arr[i]].mjd
   ep_str = {plate:plate[0:31], fiberid:fiber[0:31], mjd:mjd[0:31]}
   mk_rms_spec, rmid_arr[i], zsys_arr[i], ep_str = ep_str, result=result0, nomask=nomask

   result_all[i] = result0

   splog, 'finished: ', i+1, '/', nobj
endfor

end


; make an rms spectrum for a given RM quasar

pro mk_rms_spec, rmid, zsys, diag=diag, ep_str=ep_str, result = result, nomask = nomask


; first readin the coadded spectrum
rm_readspec,0,rmid+1, mjd=56837, wave=lam0, flux=flux0, invvar=ivar0,/silent, calibdir='wh_skysub/'
lam0 = lam0 / (1. + zsys)
line_window = [ [1160, 1340], $
                [1360, 1446], $
                [1494, 1680], $
                [1830, 1976] ]
nline = (size(line_window))[2]
mask = ivar0
flux = flux0
if ~keyword_Set(nomask) then begin ; remove emission lines
  for i=0, nline - 1 do begin
    ind = where( lam0 ge line_window[0,i] and lam0 lt line_window[1,i] )
    if ind[0] ne -1 then mask[ind] = 0
  endfor
endif

; interpolate the ivar=0 pixels
ind = where(mask eq 0, complement = indd)
if ind[0] ne -1 then flux[ind] = interpol(flux[indd],lam0[indd],lam0[ind])
 ;spline(lam0[indd], flux[indd], lam0[ind])

if keyword_set(diag) then begin
  plot, lam0, flux0
  oplot, lam0, flux, color=cgcolor('red')
endif

; set up the common restframe wavelength range
wave = 1300. + findgen(1101)*1. ; this requires 1.85<z<3.55 for full coverage
flux_mean = interpol(flux, lam0, wave)
npix = n_elements(wave)
result = {rmid:rmid,zsys:zsys,wave:wave, flux_mean:flux_mean}

; now go through all the epochs for this RMID
if keyword_set(ep_str) then begin  ; ep_str is {plate,fiber,mjd} struct
   nep = n_elements(ep_str.plate)
   result = struct_addtags(result, ep_str)
   tmp = {flux_ep:dblarr(npix,nep), flux_ep_norm:dblarr(npix,nep), flux_ep_norm_med:dblarr(npix)}
   result = struct_addtags(result, tmp)

   arr = dblarr(npix,nep) & arr_norm = dblarr(npix,nep)
   for j=0, nep - 1 do begin
       plate = (ep_str.plate)[j] & fiber = (ep_str.fiberid)[j] & mjd = (ep_str.mjd)[j]
       rm_readspec,plate,fiber, mjd=mjd, wave=lam1, flux=flux1, invvar=ivar1,/silent, $
         calibdir='wh_skysub/'
       lam1 = lam1 / (1. + zsys)
       mask1 = ivar1
       flux = flux1
       if ~keyword_Set(nomask) then begin
         for i=0, nline - 1 do begin
           ind = where( lam1 ge line_window[0,i] and lam1 lt line_window[1,i] )
           if ind[0] ne -1 then mask1[ind] = 0
         endfor
       endif
       ; interpolate the ivar=0 pixels
       ind = where(mask1 eq 0, complement = indd)
       if ind[0] ne -1 then flux[ind] = interpol(flux[indd],lam1[indd],lam1[ind])
       arr[*, j] = interpol(flux, lam1, wave)
       arr_norm[*, j] = alog10(arr[*, j]/flux_mean)  ; divide by the mean flux
       ind_norm = where(wave gt 2180 and wave lt 2220) ; normalzie at 1740A
       arr_norm[*, j] = arr_norm[*, j] / median(arr_norm[ind_norm, j])       
   endfor
   result.flux_ep = arr & result.flux_ep_norm = arr_norm
   result.flux_ep_norm_med = median(arr_norm, dim=2)

endif

end
