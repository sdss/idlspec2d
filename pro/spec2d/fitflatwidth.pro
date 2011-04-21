; NAME:
;   fitflatwidth
;
; PURPOSE:
;   Fit a traceset to the first-order corrected width of the flat field
;
; CALLING SEQUENCE:
;   widthset = fitflatwidth(flux, fluxivar, ansimage, [ fibermask, $
;    ncoeff=, sigma=, medwidth=, mask=, inmask=, double= ])
;
; INPUTS:
;   flux       - flat-field extracted flux
;   fluxivar   - corresponding inverse variance
;   ansimage   - output from extract image which contains parameter values
;
; OPTIONAL INPUTS:
;   fibermask  - nTrace bit mask, which marks bad fibers
;   ncoeff     - Order of legendre polynomial to apply to width vs. row;
;                default to 5.
;   sigma      - The SIGMA input to EXTRACT_IMAGE when determining ANSIMAGE;
;                default to 1.0 pix.  This can be a scalar, an [NFIBER] vector,
;                or an [NROW,NFIBER] array.
;
; OUTPUTS:
;   widthset   - Traceset structure containing fitted coefficients
;
; OPTION INPUTS:
;   mask      - the mask
;   inmask    - inmask for the call to xy2traceset
;   double    - input for the call to xy2traceset
;
; OPTIONAL OUTPUTS:
;   medwidth  - Median dispersion widths in each of the 4 quadrants
;               of the CCD, ordered LL,LR,UL,UR.
;
; COMMENTS:
;   The widths are forced to be the same as a function of row number
;   for all 20 fibers in each fiber bundle.
;
;   Used to fill flatstruct.widthset, which can then be applied
;   to object extraction (known profile widths).
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   xy2traceset
;
; REVISION HISTORY:
;   01-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function fitflatwidth, flux, fluxivar, ansimage, fibermask, $
 ncoeff=ncoeff, sigma=sigma, medwidth=medwidth, mask=mask, $
 inmask=inmask, double=double

   if (NOT keyword_set(ncoeff)) then ncoeff = 5
   if (NOT keyword_set(sigma)) then sigma = 1.0

   ntrace = (size(flux,/dimen))[1]
   nrow = (size(flux,/dimen))[0]
   numbundles = ntrace/20

   ;new changes so using the middle half of image
   flux2=flux[nrow/4.:3*nrow/4.-1,*]
   fluxivar2=fluxivar[nrow/4.:3*nrow/4.-1,*]
   ansimage2=ansimage[*,nrow/4.:3*nrow/4.-1]
   nrow2 = (size(flux2,/dimen))[0]

   ;----------
   ; Generate a mask of good measurements based only upon the fibermask.

   if (NOT keyword_set(mask)) then mask = (flux2 GT 0) * (fluxivar2 GT 0) 

   if (keyword_set(fibermask)) then begin
      badflats = where(fibermask NE 0)
      if (badflats[0] NE -1) then mask[*,badflats] = 0
   endif

   ;----------
   ; Determine the widths from the output array from EXTRACT_IMAGE.

   igood = where(mask)
   widthterm = transpose(ansimage2[lindgen(ntrace)*2+1,*])
   width = make_array(size=size(flux2), /float)
   if (igood[0] NE -1) then $
    width[igood] = (1 + widthterm[igood] / flux2[igood])

   ndim = size(sigma, /n_dimen)
   if (n_elements(sigma) EQ 1) then begin
      width = width * sigma[0]
   endif else if (ndim EQ 1) then begin
      for itrace=0, ntrace-1 do $
       width[*,itrace] = width[*,itrace] * sigma[itrace]
   endif else if (ndim EQ 2) then begin
      if (n_elements(sigma) NE n_elements(width)) then $
       message, 'Dimensions of SIGMA and WIDTH do not agree'
      width = width * sigma
   endif else begin
      message, 'Unsupported number of elements for SIGMA'
   endelse

   ;----------
   ; Compute the widths in each of 4 quandrants on the CCD

   ;medwidth = [ median(width[0:nrow2/2-1,0:ntrace/2-1]), $
   ;             median(width[0:nrow2/2-1,ntrace/2:ntrace-1]), $
   ;             median(width[nrow2/2:nrow-1,0:ntrace/2-1]), $
   ;             median(width[nrow2/2:nrow-1,ntrace/2:ntrace-1]) ];

   ;splog, 'Median spatial widths = ' $
   ; + string(medwidth,format='(4f5.2)') + ' pix (LL LR UL UR)'

   ;----------
   ; Perform median across bundles on good arclines only
   ; somewhat tedious, but it works

   width = reform(width,nrow2,20,numbundles)
   mask = reform(mask,nrow2,20,numbundles)
   width_bundle = fltarr(nrow2,numbundles)

   for irow=0, nrow2-1 do begin
      for j=0, numbundles-1 do begin
         ss = where(mask[irow,*,j])
         if (ss[0] NE -1) then $
          width_bundle[irow,j] = djs_median(width[irow,ss,j])
      endfor
   endfor

   width_final = rebin(width_bundle, nrow2, ntrace, /sample)

   ;----------
   ; Turn the widths back into a traceset.

   ; Generate the corresponding mask that is the same within each
   ; bundle, and marked as good if at least 25% of the points are unmasked
   if (n_elements(inmask) NE 0) then begin
      mask_bundle = rebin(float(inmask), nrow2, numbundles) GE 0.25
      mask_final = rebin(mask_bundle, nrow2, ntrace, /sample)
   endif else begin
      mask_final = width_final GT 0
   endelse

   xy2traceset, findgen(nrow2) # replicate(1,ntrace), $
    width_final, widthset, ncoeff=ncoeff, xmin=xmin, xmax=xmax, $
    inmask=mask_final, double=double

   ;----------
   ; Compute the widths in each of 4 quandrants on the CCD
   ; as the median of the unmasked pixels


   ;commented out 4-19 ;;;;;;;;;;;;;   
   ;old method using quadrant
   ;traceset2xy, widthset, xx, width_fit
   ;x1 = [0,0,nrow/2,nrow/2]
   ;x2 = [nrow/2-1,nrow/2-1,nrow-1,nrow-1]
   ;y1 = [0,ntrace/2,0,ntrace/2]
   ;y2 = [ntrace/2-1,ntrace-1,ntrace/2-1,ntrace-1]
   ;medwidth = fltarr(4)
   ;for i=0,3 do begin
   ;   indx = where(mask_final[x1[i]:x2[i],y1[i]:y2[i]],ct)
   ;   if (ct GT 0) then $
   ;    medwidth[i] = $
   ;     median([ (width_fit[x1[i]:x2[i],y1[i]:y2[i]])[indx] ])
   ;endfor
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;matt modify for 4x4 grid, use quadrupole terms
   traceset2xy, widthset, xx, width_fit
   x1 = [0,nrow2/4,nrow2/4,3*nrow2/4]
   x2 = [nrow2/4-1,3*nrow2/4-1,3*nrow2/4-1,nrow2-1]

   y1 = [ntrace/4,0,3*ntrace/4,ntrace/4]
   y2 = [3*ntrace/4-1,ntrace/4-1,ntrace-1,3*ntrace/4-1]
   medwidth = fltarr(4)
   for i=0,3 do begin
      indx = where(mask_final[x1[i]:x2[i],y1[i]:y2[i]],ct)
      if (ct GT 0) then $
       medwidth[i] = $
        median([ (width_fit[x1[i]:x2[i],y1[i]:y2[i]])[indx] ])
   endfor


   splog, 'Median spatial widths = ' $
    + string(medwidth,format='(4f5.2)') + ' pix (L B T R)';left bottom top right
   return, widthset
end
;------------------------------------------------------------------------------
