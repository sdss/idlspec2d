;+
;
; NAME:
;  bolton_mlpca
;
; PURPOSE:
;  Bolton's maximum-likelihood PCA implementation for astronomy
;
; USAGE:
;  tempset = bolton_mlpca(objflux=objflux, objivar=objivar [, $
;             nbasis=nbasis , convtest=convtest, maxiter=maxiter, $
;             coeffs=coeffs])
;
; ARGUMENTS:
;  objflux: NPIX x NOBJ array of REST-FRAME ALIGNED object spectra
;  objivar: NPIX x NOBJ inverse variance array for objflux
;
; OPTIONAL INPUTS:
;  nbasis: number of basis spectra to compute (default is 5)
;  convtest: convergence param. to test for completion (default 1.d-10)
;  maxiter: maximum number of iterations (default is 10000)
;
; RETURNS:
;  tempset: an NPIX x NBASIS array of object template basis spectra
;  coeffs: the best-fit coeffs of the input pectrum sample,
;          such that tempset # coeffs gives the minimum chi2 model
;          describing objflux
;
; WRITTEN:
;  A. Bolton, U. of Utah, 2010 May
;
; COMMENTS:
;  Implements the iterative Maximum-Likelihood PCA algorithm
;  (no intercepts, uncorrelated errors) described by
;  Wentzell et al. 1997, Journal of Chemometrics, 11, 339
;
;-

function bolton_mlpca, objflux=objflux, objivar=objivar, $
 nbasis=nbasis, convtest=convtest, maxiter=maxiter, coeffs=coeffs

; Defaults:
if (not keyword_set(nbasis)) then nbasis = 5L
if (not keyword_set(convtest)) then convtest = 1.d-10
if (not keyword_set(maxiter)) then maxiter = 10000

; Transpose and double-precisionize the data.
; Transposition is done to make these look right as matrices
; to the singular-value decompositon routines.
data = double(transpose(objflux))
ivar = double(transpose(objivar))
; Zero out zero-ivar pixels:
data = data * (ivar gt 0.d0)

; For convenience:
nobj = (size(data))[1]
npix = (size(data))[2]
b_ones = replicate(1.d0, nbasis)
coeffs = replicate(0.d0, nbasis, nobj)
; The following assignment seeds the iterative process:
model = data

niter = 0L
convergence_param = 1.d12
while ((convergence_param gt convtest) and (niter lt maxiter)) do begin
   niter++

; Truncated SVD in the object space:
   la_svd, model, w, u, v, /double
   vhat = v[0:nbasis-1,*]

; Loop over pixels to find ML model:
   for i = 0L, npix-1 do begin
      ivhat = vhat * (b_ones # ivar[*,i])
      ialpha = ivhat # transpose(vhat)
      alpha = invert(ialpha)
      beta = ivhat # data[*,i]
      model[*,i] = transpose(vhat) # (alpha # beta)
   endfor

   chisq1 = total((data - model)^2 * ivar)

; Truncated SVD in the pixel space:
   la_svd, transpose(model), w, u, v, /double
   vhat = v[0:nbasis-1,*]

; Loop over objects to find ML model:
   for i = 0, nobj-1 do begin
      ivhat = vhat * (b_ones # ivar[i,*])
      ialpha = ivhat # transpose(vhat)
      alpha = invert(ialpha)
      beta = ivhat # transpose(data[i,*])
      coeffs[*,i] = alpha # beta
;     model[i,*] = transpose(vhat) # (alpha # beta)
   endfor
   model = transpose(coeffs) # vhat

   chisq2 = total((data - model)^2 * ivar)
   convergence_param = (chisq1 - chisq2) / chisq2
   print, niter, chisq2, convergence_param
endwhile

; Tidy up by making all the basis functions non-negative:
vtot = total(vhat, 2)
signflip = double(1L - 2L * long(vtot lt 0.d0))
vhat = vhat * (signflip # replicate(1.d0, npix))
coeffs = coeffs * (signflip # replicate(1.d0, nobj))

return, transpose(vhat)
end
