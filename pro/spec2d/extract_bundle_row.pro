function lowertriangularinvert, L

n = n_elements(L[*,0])


X = dblarr(n, n) ;X is the matrix inverse of L

for i = 0, n - 1 do begin

    X[i,i] = 1d / L[i,i]

    if i lt n - 1 then begin

        for j = i + 1, n - 1 do begin

            sum = 0d
            for k = i, j - 1 do sum = sum - L[k,j] * X[i,k]
            X[i,j] = sum / L[j,j]

        endfor

    endif

endfor

return, X
end

;+
; NAME:
;   extract_bundle_row
;
; PURPOSE:
;   Fit the fiber profiles and background in a single row with least
;   squares, chunking it up by bundle, using Gaussians plus polynomial
;   background.
;
;   Based upon extract_row.pro, which was a call to extract_row.c
;
; NOTES ON CONVERSION TO EXTRACT_BUNDLE_ROW (A. Bolton, Utah, 2011Aug):
;   The extraction guts of this program are entirely new,
;   written in pure IDL rather than in IDL-wrapped C.
;   However, I have attempted (1) to preserve most of the same
;   outward functionality in terms of inputs and outputs, at least
;   those that are relevant to the BOSS case, and (2) to retain
;   the same rejection-iteration behavior, so that the same sorts
;   of cases trigger rejection as before.  I did not parse all
;   the details of the rejection logic.  One thing I know is that
;   the rejection in current form would be more efficient if it were
;   done on each bundle separately, rather than on the entire row.
;   But that would require a fine-grained refactoring that would risk
;   breaking the rejection logic, so I didn't attempt it.
;
;   Some deprecated (and even retained) keyword input/output variables
;   may have unexpected behavior is used in a context other than the
;   current BOSS-implementation calls in EXTRACT_IMAGE.
;
;   I have retained a fair bit of the original code in commented-out
;   form, for forensic purposes.  Anything removed during the
;   conversion to bundle form is commented out with a ";#".
;
; CALLING SEQUENCE:
;; First four positional arguments inherited from extract_row, plus rdnoise:
;   ans = extract_row( fimage, invvar, rdnoise, xcen, sigma,
;; Next keyword arguments new to extract_bundle_row:
;; (skew and kurt not currently enabled)
;              [nperbun=, buffsize=, skew=, kurt=,
;; Next keyword retained but reinterpreted:
;              npoly=, 
;; Next keywords retained in more or less the same role (I think):
;              maxiter=, lowrej=, highrej=, niter=, reducedChi=, xvar=,
;              mask=, relative=, diagonal=, pixelmask=, reject=,, ymodel=
;; Next keyword arguments deprecated:
;              fscat=, proftype=, wfixed=, inputans=, iback=, bfixarr=,
;              fullcovar=, wfixarr=, squashprofile=, whopping=,
;              wsigma=, nBand=  ])
;
; INPUTS:
;   fimage     - Vector [nCol]
;   invvar     - Inverse variance [nCol]
;   rdnoise    - Readout noise [nCol]
;   xcen       - Initial guesses for X centers [nFiber]
;   sigma      - Sigma of gaussian profile; (scalar or [nFiber])
;
; OPTIONAL KEYWORDS (new):
;   nperbun    - Number of fibers per bundle (default is 20).
;   buffsize   - Pixel buffer size on the very outside of the image
;                (default is 8)
;   npoly      - order of polynomial background in bundle, default=2 (linear).
;
; OPTIONAL KEYWORDS (retained):
;   maxiter    - Maximum number of profile fitting iterations; default to 10
;   lowrej     - Negative sigma deviation to be rejected; default to 5
;   highrej    - Positive sigma deviation to be rejected; default to 5
;   reject     - Three-element array setting partial and full rejection
;                thresholds for profiles; default [0.2, 0.6, 0.6].
;                When there is less than REJECT[2] of the area is left,
;                  then drop fitting of all higher-order terms.
;                When there is less than REJECT[1] of the area is left,
;                  then the pixel is rejected (inverse variance is set to 0).
;                When there is less than REJECT[0] of the area is left,
;                  then assume that there's no fiber there, and don't fit
;                  for that fiber at all.
;   relative   - Set to use reduced chi-square to scale rejection threshold
;
; DEPRECATED KEYWORDS:
;   inputans   - Input fit, excluding background and whopping terms
;                [ncoeff*nFiber]
;                The array is sorted as follows:
;                  [ncoeff] values for fiber #0
;                   ...
;                  [ncoeff] values for fiber #(nFiber-1)
;   squashprofile - ???
;   nband      - Band-width of covariance matrix of fiber profiles: default 1
;   whopping   - X locations to center additional "whopping" terms to describe
;                the exponentail tails of flux near bright fibers; default
;                to -1, which means not to use any such terms.
;   wsigma     - Sigma width for exponential whopping profiles; default to 25
;
; MODIFIED INPUTS (OPTIONAL, retained):
;   xvar       - X values of fimage and invvar; default is findgen(NX).
;   mask       - Image mask: 1=good, 0=bad [NX]
;   pixelmask  - Bits set for each fiber due to extraction rejection
;                [nFiber]
;
; MODIFIED INPUTS (OPTIONAL, deprecated):
;   wfixed     - Array to describe which parameters to fix in the profile;
;                0=fixed, 1=float; default to [1].
;                The number of parameters to fit per fiber is determined
;                this way; e.g. nCoeff = n_elements(wfixed), so the default
;                is to fit only 1 parameter per fiber.  For the (default)
;                Gaussian profile, this is the height of the Gaussian.
;                Note that WFIXED is used to build the array WFIXARR.
;   iback      - 1D array of input background coeff 
;                (needed if fixed parameters are non-zero)
;   bfixarr    - 1D integer array to specify which terms of the background
;                coefficients to fix; 0=fixed, 1=float.
;   wfixarr    - 1D integer array to specify which parameters in the full fit
;                to fix; 0=fixed, 1=float.
;                The array is sorted as follows:
;                  [ncoeff] values for fiber #0
;                   ...
;                  [ncoeff] values for fiber #(nFiber-1)
;                  [npoly] values for the background polynomial terms
;                  [whoppingct] values for the whopping terms
;
; OUTPUTS (reinterpreted):
;   ans        - Output fit [nFiber]
;                (extracted Gaussian profile amplitudes)
;
; OPTIONAL OUTPUTS (retained):
;   ymodel     - Evaluation of best fit [nCol]
;   ybkg       - Background par of best fit [nCol]
;   diagonal   - 1D diagonal of covariance matrix
;   niter      - Number of rejection iterations performed
;   reducedChi - Reduced chi ???
;
; OPTIONAL OUTPUTS (deprecated):
;   fscat      - Scattered light contribution in each fiber [nFiber]
;   fullcovar  - 2D covariance matrix.  This is a symmetric matrix, and we
;                only fill the lower triangle.  Computing this increases CPU
;                time by a factor of 2 or 3.
;
; REVISION HISTORY:
;    8-Aug-1999  extract_row Written by Scott Burles, Chicago 
;    Sep 2010    converted to extract_bundle_row by Adam S. Bolton, Utah
;    Aug 2011    attempted documentation cleanup, Adam S. Bolton, Utah
;    Dec 2014    Added iterative inverse variance model update, S Bailey LBL
;-
;------------------------------------------------------------------------------
function extract_bundle_row, fimage, invvar, rdnoise, xcen, sigma, ymodel=ymodel, $
 ybkg=ybkg,fscat=fscat, proftype=proftype, wfixed=wfixed, inputans=inputans, $
 iback=iback, bfixarr=bfixarr, xvar=xvar, mask=mask, relative=relative, $
 squashprofile=squashprofile, diagonal=p, fullcovar=fullcovar, $
 wfixarr=wfixarr, npoly=npoly, maxiter=maxiter, $
 lowrej=lowrej, highrej=highrej, niter=niter, reducedChi=reducedChi, $
 whopping=whopping, wsigma=wsigma, pixelmask=pixelmask, reject=reject, $
 oldreject=oldreject, nband = nband, contribution=contribution, $
 buffsize=buffsize, skew=skew, kurt=kurt, chi2pdf=chi2pdf, $
 use_image_ivar=use_image_ivar, nbundles=nbundles, bundlefibers=bundlefibers

   ; Need 4 parameters
   if (N_params() LT 4) then $
    message, 'Wrong number of parameters'

   ntrace = n_elements(xcen)
   nx = n_elements(fimage)

   if (n_elements(sigma) NE ntrace) then begin
      sigma1 = sigma[0]
      sigma = xcen*0.0 + sigma1
   endif 

   ; extract_row keywords:
   if (n_elements(npoly) EQ 0) then npoly = 2L ; order of background is now per bundle!
   if (NOT keyword_set(nband)) then nband = 1L
   if (NOT keyword_set(maxiter)) then maxiter = 50
   if (NOT keyword_set(highrej)) then highrej = 15.0
   if (NOT keyword_set(lowrej)) then lowrej = 20.0
   if (NOT keyword_set(wfixed)) then wfixed = [1]
   if (NOT keyword_set(proftype)) then proftype = 1 ; Gaussian
   relative = keyword_set(relative) 
   squashprofile = keyword_set(squashprofile) 
   if (NOT keyword_set(wsigma)) then wsigma = 25.0
   ; extract_bundle_row keywords:
   if (NOT keyword_set(buffsize)) then buffsize = 8L

   ; Here we want a three element array where both are between 0 and 1, and the
   ; first is larger than the second.  
   ; The first threshold sets the minimum area required to perform 
   ;    single parameter profile fitting
   ; The second threshold is the minimum area required not to reject 
   ;    the pixel in the final extracted spectrum.
   ; The third parameter is the area required
   ;    in the profile fit containing good pixels to do a full fit.

   if (n_elements(reject) NE 3) then reject = [0.2, 0.6, 0.6]

   if (n_elements(pixelmask) NE ntrace $
    OR size(pixelmask,/tname) NE 'LONG') then $
    pixelmask = lonarr(ntrace)

   if (NOT keyword_set(whopping)) then whopping = -1
   if (whopping[0] EQ -1) then whoppingct = 0L $
    else whoppingct = n_elements(whopping)

   if (NOT keyword_set(xvar)) then xvar = findgen(nx) $
    else if (nx NE n_elements(xvar)) then $
     message, 'Number of elements in FIMAGE and XVAR must be equal'

   if (NOT keyword_set(mask)) then mask = bytarr(nx) + 1b $
    else if (nx NE n_elements(mask)) then $
     message, 'Number of elements in FIMAGE and MASK must be equal'

   ncoeff = n_elements(wfixed)

   if (nx NE n_elements(invvar)) then $
    message, 'Number of elements in FIMAGE and INVVAR must be equal'
    
   ; if rdnoise is just a scalar, convert to vector
   if n_elements(rdnoise) EQ 1 then $
    rdnoise = fltarr(nx) + rdnoise[0]

   ;----------
   ; Allocate memory for the C subroutine.

   ymodel = fltarr(nx)
   ybkg   = fltarr(nx)

   ;----------
   ; Test which points are good

   qgood = invvar GT 0.0 AND mask NE 0
   igood = where(qgood, ngood)

   ;----------
   ; Set the following variables before any possible RETURN statement.

   reducedChi = 0.0
   niter = 0
   ans = fltarr(ntrace)
   p = fltarr(ntrace)

   if (ngood EQ 0) then return, ans

   ;----------
   ; Check that XCEN is sorted in increasing order
   ; with separations of at least 3 pixels.

   junk = where(xcen[0:ntrace-2] GE xcen[1:ntrace-1] - 3, ct)
   if (ct GT 0) then $
     splog, 'XCEN is not sorted or not separated by greater than 3 pixels.'

   ;----------
   ; Build the fixed parameter array if it was not passed.

   wfixarr = lonarr(ntrace) + 1

   ; PUTTING IN LOOP OVER BUNDLES:
   ; Find number of bundles:
   nbun = nbundles
   ; Find breaks between the bundles
   ; and set limits of pixels to be associated with each bundle:
   t_lo = intarr(nbundles)
   for ibun=1, nbundles-1 do t_lo[ibun]=total(bundlefibers[0:ibun-1])
   t_hi = t_lo + bundlefibers -1
   xc_lo = xcen[t_lo]
   xc_hi = xcen[t_hi]
   if (nbun gt 1) then begin
      midpoints = 0.5 * (xc_lo[1:nbun-1] + xc_hi[0:nbun-2])
      jmax = floor(midpoints)
      jmin = jmax + 1
      jmin = [(floor(xc_lo[0]) - buffsize) > 0, jmin]
      jmax = [jmax, (ceil(xc_hi[nbun-1]) + buffsize) < (nx - 1)]
   endif else begin
      jmin = (floor(xc_lo) - buffsize) > 0
      jmax = (ceil(xc_hi) + buffsize) < (nx - 1)
   endelse

   ; JG
   chi2pdf=fltarr(ntrace)

   ; Mask outside the bundle-fitting zone:
   mask = mask * (xvar ge min(jmin)) * (xvar le max(jmax))

   ; The loop over bundles:
   totalreject = 0
   partial = lonarr(ntrace)
   fullreject = lonarr(ntrace)
   finished = 0

   ; JG invert ordering of loops : first bundles then iterative clipping
   ; per bundle

   ; JG : use only readout noise and apply mask by default (other option
   ; for arc fit because more robust)
   if(keyword_set(use_image_ivar)) then $ ; JG   
      workinvvar = (FLOAT(invvar * mask)) $ ; JG   
   else workinvvar = (FLOAT(mask*(rdnoise gt 0.)*(invvar gt 0.)/(rdnoise^2+(rdnoise eq 0.)))) ; JG   
   
   ; loop on bundle

   
   start = 0
   for ibun = 0L, nbun-1 do begin
       fiberbase = findgen(bundlefibers[ibun])
       nperbun = bundlefibers[ibun]
       endf = start + nperbun

       ; Make the extracting basis:
       workxcen = xcen[t_lo[ibun]:t_hi[ibun]]
       worksigma = sigma[t_lo[ibun]:t_hi[ibun]]
       workpar = (transpose([[workxcen], [worksigma]]))[*]
       workx = xvar[jmin[ibun]:jmax[ibun]]
       workbasis = gausspix(workx, workpar)
       workimage = fimage[jmin[ibun]:jmax[ibun]]
       
       workmodel = 0.*workimage

       ; Associate pixels with fibers:
       pixelfiber = round(interpol(fiberbase, workxcen, workx))
       if bundlefibers[ibun] eq 1 then pixelfiber= round(replicate(fiberbase, n_elements(workx)))
       
       bworkinvvar = workinvvar[jmin[ibun]:jmax[ibun]]
       bworkrdnoise = rdnoise[jmin[ibun]:jmax[ibun]]      

       workbkg   = 0.*workimage
       if (npoly gt 0) then begin
          ; Fit the background only on pixels on the edges of the bundle to
          ; avoid biasing the background if the PSF shape is inaccurate

          if (npoly gt 2) then begin
             splog, "ERROR in extract_bundle_row, npoly must be <=2 here"
             splog, "because we are fitting the polynomial only on the edges of a bundle"
             STOP
          endif
          
          bkg_invvar = (bworkinvvar gt 0 )*( (workx lt (workxcen[0]-3*worksigma[0]) ) or (workx gt (workxcen[nperbun-1]+3*worksigma[nperbun-1]) ) )
          workxpoly   = 2. * (workx - min(workx)) / (max(workx) - min(workx)) - 1.
          bkg_basis   = flegendre(workxpoly, npoly)
          bkg_itbasis = transpose(bkg_basis * (bkg_invvar # replicate(1., npoly)))
          icovar      = bkg_itbasis # bkg_basis
          beta        = bkg_itbasis # workimage
          L           = icovar
          la_choldc, L, status=cholstat
          if (cholstat eq 0) then begin                                
             bkg_coeffs = la_cholsol(L, beta)
             workbkg    = bkg_basis # bkg_coeffs ; this is the array of background values
             ybkg[jmin[ibun]:jmax[ibun]] = workbkg ; save this
          endif
          
       endif



       number_of_pixel_rejected_in_bundle_row = 0

       finished = 0
       niter = 0

       saved_workbasis=workbasis

       while(finished NE 1) do begin 
           
           ; need this because altered below
           workbasis=saved_workbasis

           ; Determine the good area fraction for each fiber,
           ; limiting to pixels "associated with" that fiber:
           goodarea = 0. * fltarr(nperbun)
           for ijk = 0L, nperbun-1 do begin
               totalarea = total((pixelfiber eq ijk) * (workbasis[*,ijk]))
               totalgood = total((pixelfiber eq ijk) * (bworkinvvar gt 0.) * (workbasis[*,ijk]))
               if (totalarea gt 0.) then goodarea[ijk] = totalgood / totalarea
           endfor
           
        
           ; Populate the rejection masks:
           partial[t_lo[ibun]:t_hi[ibun]] = goodarea lt reject[2]
           fullreject[t_lo[ibun]:t_hi[ibun]] = goodarea lt reject[1]
           
           ; Remove fitting for any fiber with area less than reject[0]:
           use_component = where(goodarea ge reject[0], n_use)
           
           ; Don't bother fitting at all if there are no unmasked fibers:
           if (n_use eq 0) then begin 
              break
           endif 
           
           nfit_this = n_elements(use_component)
           workbasis = workbasis[*,use_component]
           
           ; The extraction steps:
           itworkbasis = transpose(workbasis * (bworkinvvar # replicate(1., nfit_this)))
           icovar = itworkbasis # workbasis
           beta = itworkbasis # (workimage - workbkg)
           workidx = lindgen(nfit_this)
           


         
           ; JG : Evaluate the variance of the
           ; JG : fluxes where we used only readoutnoise.
           ; JG : this estimator is unbiased, it is optimal at low S/N and suboptimal
           ; JG : at high S/N. however for a number of electrons per fiber per row <
           ; JG : 1000, the errors are only 10% larger than the optimal case
           ;
           ; JG : here, icovar is not anymore the inverse of the covariance matrix
           ; JG : because the estimator is not optimal. we need to reevaluate the errors
           ;
           ; JG : '#' is the matrix multiplication symbol
           ; JG : the fluxes are in the workcoeffs along with the 'background' polynomial coefficients.
           ;
           ; JG : workcoeffs = la_cholsol(icovar, beta)
           ; JG : workcoeffs = covar # beta
           ; JG : workcoeffs = covar # itworkbasis # workimage
           ;
           ; JG : var_workcoeffs = (covar # itworkbasis) # workimage # workimage_t # transpose(covar # itworkbasis)
           ; JG : var_workcoeffs = (covar # itworkbasis) # var_workimage # transpose(covar # itworkbasis)
           ; JG : var_workcoeffs = (covar # itworkbasis) # (rdnoise^2+workmodel) # transpose(covar # itworkbasis)
           ;
           ; JG : first need to invert icovar after la_cholsol ... and as in python
           ; JG : there is no idl routine for this, we write it above (lowertriangularinvert)
           if ( n_use gt 1 ) then begin
           ; JG : cholesky decomposition of icovar
               L = icovar
               la_choldc, L, status=cholstat
               ; JG : now icovar = transpose(L)#L
           
               if (cholstat eq 0) then begin
                   ; JG : solve linear system
                   workcoeffs = la_cholsol(L, beta)
                   ; JG : set upper terms of lower triangular matrix L to zero
                   n = n_elements(L[*,0])
                   for i = 0,n - 2 Do L[i+1:*,i] = 0
                   ; JG : compute inverse of L = Li
                   Li = lowertriangularinvert(L)
                   ; JG : because icovar=transpose(L)#L , covar=Li # transpose(Li)
                   covar = Li # transpose(Li)
                   ; JG : the pixel model
                   workmodel = workbasis # workcoeffs + workbkg
                   ; JG : compute diagonal matrix of pixel values including this time readout
                   ; JG : noise and poisson noise
                   pixel_covvar = DIAG_MATRIX(bworkrdnoise^2+workmodel*(workmodel gt 0))
                   ; JG : here is the covariance of the parameters
                   tmp_matrix = covar # itworkbasis
                   true_parameter_covar = tmp_matrix # pixel_covvar # transpose(tmp_matrix)
                   ; JG : keep only the diagonal
                   true_parameter_var = DIAG_MATRIX(true_parameter_covar)
                   ; JG : keep only the part that correspond to fluxes (discard background terms)
                   extractvar = true_parameter_var[workidx,workidx]
                   extractivar = (extractvar gt 0.)/(extractvar+(extractvar eq 0.))
               
               endif else begin
                   splog , "cholesky failed"
                   splog , "JG", size(use_component)
                   splog , "JG", nfit_this
                   splog , "JG", size(workbasis)
                   splog , "JG", size(itworkbasis)
                   splog , "JG", size(bworkinvvar)
                   splog , "JG", size(icovar)

                   ; JG : here the extraction has failed
                   workcoeffs = replicate(0., nfit_this)
                   workmodel = workbasis # workcoeffs + workbkg
                   extractivar[*] = 0.
               endelse
           endif else begin ; n_use=1, la_cholsol failed with that
               workcoeffs = beta/icovar
               covar = 1/icovar
               workmodel = workbasis # workcoeffs + workbkg
               pixel_covvar = DIAG_MATRIX(bworkrdnoise^2+workmodel*(workmodel gt 0))
               tmp_matrix = covar # itworkbasis
               true_parameter_covar = tmp_matrix # pixel_covvar # transpose(tmp_matrix)
               true_parameter_var = DIAG_MATRIX(true_parameter_covar)
               extractvar = true_parameter_var[workidx,workidx]
               extractivar = (extractvar gt 0.)/(extractvar+(extractvar eq 0.))
           endelse


           ; JG : now start outlier rejection
           ; JG : for chi2 cut replace workinvar by expected variance in pixel
           ; JG : here we assume the optimal variance
           ; JG : add a generous uncertainty on psf shape that allows to lower
           ; JG : the chi2 rejection threshold (tested with nsig=4)
           ; JG : we are only working with one bundle here
           if(keyword_set(use_image_ivar)) then $
               diffinvvar=FLOAT(bworkinvvar gt 0.)/(bworkrdnoise^2+workmodel*(workmodel gt 0.)+(0.05*workmodel)^2*(workmodel gt 0.)) $
           else diffinvvar=FLOAT(bworkinvvar gt 0.)*invvar[jmin[ibun]:jmax[ibun]]

           diffs = (workimage - workmodel) * sqrt(diffinvvar)
           chi2_array = diffs^2
           chisq = total(diffs^2)
           ndata =  total(diffs^2 gt 0)
           errscale = 1.0
           if ((relative) and (ndata gt 0)) then errscale = sqrt((chisq/ndata) > 1.0)
           
           finished = 1

           goodi = where(diffinvvar GT 0, ngoodi)
           if (ngoodi gt 0) then begin 
               indchi = diffs[goodi] / (highrej*errscale)
               neg = where(indchi LT 0, nneg)
               if (nneg gt 0) then $
                 indchi[neg] = -indchi[neg] * highrej / lowrej

               ; JG : indchi is positive, if it's greater than one we discard the pixel
               ; JG : actually we only mask the worst outlier pixel and then refit
               ; JG : because this outlier could have altered significantly the model
               worstdiff = max(indchi,worst)
               if ( worstdiff GT 1.0 ) then begin
                   pixel_index=goodi[worst]
                   ; print ," JG: DEBUG outlier pixel for bundle",ibun,"index=",pixel_index,"CCD col=",jmin[ibun]+pixel_index,"val=",workimage[pixel_index],"model=",workmodel[pixel_index],"iter=",niter
                   finished = 0
                   bworkinvvar[pixel_index] = 0
                   mask[jmin[ibun]+pixel_index] = 0
                   totalreject = totalreject + 1
                   number_of_pixel_rejected_in_bundle_row = number_of_pixel_rejected_in_bundle_row + 1
               endif
           endif

           niter = niter + 1
           if (niter EQ maxiter) then finished = 1

       ; JG : end of iteration loop on one bundle
       endwhile


       ; Unpack results, accounting for masked fibers:
       fullcoeffs = replicate(0., nperbun)
       fullexivar = replicate(0., nperbun)
       if (n_use gt 0) then begin 
          fullcoeffs[use_component] = workcoeffs
          fullexivar[use_component] = extractivar
          ; Populate the model:
          ymodel[jmin[ibun]:jmax[ibun]] = workmodel
          
          ; JG : 'ans' is the extracted flux
          ans[t_lo[ibun]:t_hi[ibun]] = fullcoeffs[0:nperbun-1]
          ; JG : 'p' is sqrt(inverse_variance), see extract_bundle_image.pro
          p[t_lo[ibun]:t_hi[ibun]] = sqrt(fullexivar[0:nperbun-1])
          
          ; JG : reduced chi2
          for ijk = 0L, nperbun-1 do begin
             sum_chi2   = total( chi2_array * (pixelfiber eq ijk) )
             sum_ndata  = total( (chi2_array gt 0) * (pixelfiber eq ijk) )
             sum_nparam = 1
             ndf=((sum_ndata-sum_nparam)>1)
             chi2pdf[start+ijk]=sum_chi2/ndf
          endfor
       endif
   start = endf

   ; JG : end of loop on bundles
   endfor


; Apply fullreject mask to output flux and inverse-error:
p = p * (fullreject eq 0)
ans = ans * (fullreject eq 0)

; Add bits to PIXELMASK
pixelmask = pixelmask OR (pixelmask_bits('PARTIALREJECT') * fix(partial))
pixelmask = pixelmask OR (pixelmask_bits('FULLREJECT') * fix(fullreject))


return, ans

end

