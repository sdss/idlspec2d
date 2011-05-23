;+
;
; NAME: trace_sparse_crude
;
; PURPOSE:
;  Semi-crude, semi-robust tracing of sparse flats.
;
; USAGE:
;   trace_sparse_crude, image, xpos, ypos [, csize=csize, $
;    maxpeaks=maxpeaks, pthresh=pthresh, fthresh=fthresh, $
;    dxfit=dxfit, sigma_guess=sigma_guess, maxdev=maxdev, $
;    ncoeff=ncoeff]
;
; ARGUMENTS:
;   image: flatfield image array with traces running
;     approx. vertically.  Best to * by (invvar gt 0.)
;     to mask known artifacts/bad columns from SDSSPROC.
;
; OPTIONAL ARGUMENTS:
;   csize: chunk size for dividing up image (def = 50)
;   maxpeaks: max # of traces to find (def = 40)
;   pthresh, fthresh: peaks must be fthresh times the
;    pthresh quantile in the row x-sec'n to be detected.
;    (def: pthresh = 0.975, fthresh = 1.0)
;   dxfit: fit peak pos'n with +/- dxfit pixels data (def = 7)
;   sigma_guess: Gaussian sigma guess (def = 1.)
;   maxdev: maximum pixel deviation for tracefit (def = 0.1)
;   ncoeff: poly order for tracefit (def = 5)
;
; WRITTEN:
;  A. Bolton @ Utah 2010may
;
;-

; Supporting pixel-integrated Gaussian function:
function gausspix_ampli, x, par
  return, par[0] * (gaussint((x - par[1] + 0.5d0) / par[2]) - gaussint((x - par[1] - 0.5d0) / par[2]))
end

pro trace_sparse_crude, image, xpos, ypos, csize=csize, maxpeaks=maxpeaks, pthresh=pthresh, $
 fthresh=fthresh, dxfit=dxfit, sigma_guess=sigma_guess, maxdev=maxdev, ncoeff=ncoeff

xpos = 0
ypos = 0

; Defaults:
if (not keyword_set(csize)) then csize = 50L
if (not keyword_set(maxpeaks)) then maxpeaks = 40L
;if (not keyword_set(cthresh)) then cthresh = 0.975
if (not keyword_set(pthresh)) then pthresh = 0.975
if (not keyword_set(fthresh)) then fthresh = 1.
if (not keyword_set(dxfit)) then dxfit = 7L
if (not keyword_set(sigma_guess)) then sigma_guess = 1.0
if (not keyword_set(maxdev)) then maxdev = 0.1
if (not keyword_set(ncoeff)) then ncoeff = 5L

nx = (size(image))[1]
ny = (size(image))[2]

; Divide into chunks:
nchunk = ny / csize
ylo = csize * lindgen(nchunk)
yhi = ylo + csize - 1
ychunk = 0.5 * (float(ylo) + float(yhi))
xpeak = fltarr(maxpeaks, nchunk)
npeakvec = replicate(-1L, nchunk)

; Loop over chunks:
splog, 'Peakfinding in chunk-wise medians.'
med_vectors = fltarr(nx, nchunk)
for j = 0L, nchunk-1 do begin
  med_csec = median(image[*,ylo[j]:yhi[j]], dimen=2)
  med_vectors[*,j] = med_csec
  med_srt = abs(med_csec[sort(abs(med_csec))])
  cval = fthresh * med_srt[long(pthresh * nx)]
  ipeak = where((med_csec gt shift(med_csec, -1)) and (med_csec ge shift(med_csec, 1)) $
                and (med_csec gt cval), npeak)
  npeakvec[j] = npeak
  npeak = npeak < maxpeaks
  if (npeak gt 0) then xpeak[0:npeak-1,j] = ipeak[0:npeak-1]
endfor

; Find the mode of the number of peaks, and call that the "real" number:
histo = histogram(npeakvec, locations=histbase)
maxhist = max(histo, wmax)
nreal = histbase[wmax]

; Reduce the arrays to those chunks containing all real peaks:
goodchunk = where(npeakvec eq nreal,  ngood)

; Only look at the "good" chunks:
ychunk = ychunk[goodchunk]
xpeak = xpeak[0:nreal-1,goodchunk]
ypeak = replicate(1., nreal) # ychunk
ylo = ylo[goodchunk]
yhi = yhi[goodchunk]
med_vectors = med_vectors[*,goodchunk]

; Sub-pixel fitting for peak locations:
xpeak_new = 0. * xpeak
splog, 'Fitting for sub-pixel peak locations.'
for j = 0L, ngood - 1 do begin
  for i = 0L, nreal - 1 do begin
    xlo = long(round(xpeak[i,j] - dxfit)) > 0
    xhi = long(round(xpeak[i,j] + dxfit)) < (nx - 1)
    data = med_vectors[xlo:xhi,j]
    xbase = findgen(n_elements(data))
    pguess = [total(data), mean(xbase), sigma_guess]
    err = 1. + 0. * xbase
    pfit = mpfitfun('gausspix_ampli', xbase, data, err, pguess, /quiet)
    xpeak_new[i,j] = xlo + pfit[1]
   endfor
endfor

; Fit traceset:
splog, 'Fitting smooth trace position functions.'
xy2traceset, transpose(ypeak), transpose(xpeak_new), tset, ncoeff=ncoeff, maxdev=maxdev

; Make new yset:
ypos = (findgen(max(yhi) - min(ylo) + 1) + min(ylo)) # replicate(1., nreal)

; Compute xset from it:
traceset2xy, tset, ypos, xpos

;xpos = xpos_ret

;stop

return
end

