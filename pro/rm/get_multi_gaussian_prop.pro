; get the peak, mode of a mixture of gauasians
; note that the xaxis is in units of ln(lambda)


function get_multi_gaussian_prop, para, diet=diet,plot_check = plot_check

n_gauss = n_elements(para)/3L

cen = dblarr(n_gauss) & sig = dblarr(n_gauss) 
area = dblarr(n_gauss)

for i=0L, n_gauss - 1L do begin

  ; flip the sign so the flux is always positive
  para[i*3L] = abs(para[i*3L]) 

  cen[i] = para[i*3L+1L] & sig[i] = para[i*3L + 2L]
  area[i] = para[i*3L]*exp(para[i*3L+1L])
endfor
area = total(area, /double)

; define an array
left = min(cen - 3.*sig) & right = max(cen + 3.*sig)
disp = 1d-5  ; 3 km/s
npix = round( (right - left)/disp )

xarr = left + dindgen(npix+1)*disp
yarr = manygauss(xarr, para, nline = n_gauss)
; add one more step to remove essentially zero fluxes
ind_non_zero = where(yarr gt 1d-6)
if ind_non_zero[0] ne -1 then begin
  xarr = xarr[ind_non_zero] & yarr = yarr[ind_non_zero]
endif

; now get the peak location
y_peak = max(yarr, ipeak)
x_peak = xarr[ipeak]
if area lt 1d-6 then begin
  x_peak = mean(cen)  ;
  x_cen_50per = x_peak
endif

if ipeak gt 0 and ipeak lt npix then begin
  fwhm_left = spline(yarr[0:ipeak-1L],xarr[0:ipeak - 1L], 0.5*y_peak)
  fwhm_right = spline(reverse(yarr[ipeak:*]), reverse(xarr[ipeak:*]), 0.5*y_peak)
  fwhm = fwhm_right - fwhm_left
  ; determine the centroid of pixels above 0.5*f_peak
  ind_cen=where(xarr ge fwhm_left and xarr le fwhm_right)
  ; remember that the array is uniform in alog(lambda) not uniform in lambda
  x_cen_50per=alog( total( (exp(xarr[ind_cen]))^2*yarr[ind_cen],/double) / total(exp(xarr[ind_cen])*yarr[ind_cen],/double) )

  if not keyword_set(diet) then begin
     ; full-width-at-third-maximum
     fwtm_left = spline(yarr[0:ipeak-1L],xarr[0:ipeak - 1L], y_peak/3.D)
     fwtm_right = spline(reverse(yarr[ipeak:*]), reverse(xarr[ipeak:*]), y_peak/3.D)
     fwtm = fwtm_right - fwtm_left

     ; full-width-at-quarter-maximum
     fwqm_left = spline(yarr[0:ipeak-1L],xarr[0:ipeak - 1L], 0.25*y_peak)
     fwqm_right = spline(reverse(yarr[ipeak:*]), reverse(xarr[ipeak:*]), 0.25*y_peak)
     fwqm = fwqm_right - fwqm_left
  endif
endif else begin
  fwhm_left = 0.d & fwhm_right = 0.d & fwhm = 0.d
  fwtm_left = 0.d & fwtm_right = 0.d & fwtm = 0.d
  fwqm_left = 0.d & fwqm_right = 0.d & fwqm = 0.d
endelse

if not keyword_set(diet) then begin
   ; now get the log_wave that separate the broad line into equal halves
   product = exp(xarr)*yarr*disp
   area_left = total(product, /cumulative, /double)
   ; print, area, min(area_left), max(area_left)
   area_diff = abs(area_left - 0.5*area)
   temp = min(area_diff, ind_min)
   x_half_flux = xarr[ind_min]
endif

if not keyword_set(diet) then begin
   result = [x_peak, fwhm, area, fwhm_left, fwhm_right, $
         y_peak, fwtm, fwtm_left, fwtm_right, fwqm, fwqm_left, fwqm_right, x_half_flux, x_cen_50per]
endif else result = [x_peak, fwhm, area, fwhm_left, fwhm_right, x_cen_50per]

if keyword_set(plot_check) then begin
  plot, xarr, yarr
  for i=0L, n_gauss - 1L do oplot, xarr, onegauss(xarr, para[3L*i:3L*i+2]), color=fsc_color('green')
  oplot, [x_peak, x_peak], [0, y_peak*1.05], linestyle = 2
  oplot, [fwhm_left, fwhm_right], 0.5*[y_peak, y_peak]
  oplot, [x_half_flux, x_half_flux], [0, y_peak*1.05], linestyle = 1
endif

return, result
end
