;+
; NAME:
;   combine1fiber
;
; PURPOSE:
;   Combine several spectra of the same object, or resample a spectra.
;
; CALLING SEQUENCE:
;   combine1fiber, inloglam, objflux, [ objivar, finalmask=, indisp=, $
;    newloglam=, newflux=, newivar=, andmask=, ormask=, newdisp=, $
;    nord=, binsz=, bkptbin=, maxsep= ]
;
; INPUTS:
;   inloglam       - Wavelengths in log10-Angstroms [NPIX,NSPEC]
;   objflux        - Flux [NPIX,NSPEC]
;
; REQUIRED KEYWORDS:
;   newloglam      - Wavelengths for output evaluation, also in log10-Angstroms
;                    [NPIX,NSPEC]
;
; OPTIONAL INPUTS:
;   objivar        - Inverse variance [NPIX,NSPEC]
;   finalmask      - Pixel mask [NPIX,NSPEC]
;   indisp         - Dispersion values [NPIX,NSPEC]
;   binsz          - Bin separation for INLOGLAM; if not set, then default
;                    to INLOGLAM[1]-INLOGLAM[0]. [NPIX,NSPEC]
;   nord           - Order of spline fit; default to 3.
;   bkptbin        - Break point binning; default to 1.2 * BINSZ.
;   maxsep         - Maximum separation between input wavelengths.  The spline
;                    fit is split into pieces, with the breaks wherever this
;                    spacing is exceeded.  Default to 2.0 * BINSZ.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   newflux        - Resampled flux.
;   newivar        - Resampled inverse variance.
;   andmask        - Resampled mask. For each mask bit, set that bit only if
;                    every input spectrum at this wavelength has that bit set
;                    (e.g., this is a logical AND).
;   ormask         - Resampled mask. For each mask bit, set that bit if any
;                    of the input spectra at this wavelength has that bit set
;                    (e.g., this is a logical OR).
;   newdisp        - Resampled dispersion values.
;
; COMMENTS:
;   One can pass this routine a single spectrum to be fit by a spline and
;   re-sampled, in which case all the inputs (such as FLUX) are 1-dimensional
;   arrays.  Or, one can pass it several spectra, in which case these inputs
;   are 2-dimensional arrays.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_laxisgen()
;   djs_maskinterp()
;   djs_median()
;   pixelmask_bits()
;
; REVISION HISTORY:
;   02-Jan-2000  Written by D. Schlegel; modified from COMBINE2DOUT
;-
;------------------------------------------------------------------------------
pro combine1fiber, inloglam, objflux, objivar, $
 finalmask=finalmask, indisp=indisp, $
 newloglam=newloglam, newflux=newflux, newivar=newivar, $
 andmask=andmask, ormask=ormask, newdisp=newdisp, $
 nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep

   ;----------
   ; Check that dimensions of inputs are valid

   npix = n_elements(inloglam)
   nfinalpix = n_elements(newloglam)

   if (npix EQ 0 OR nfinalpix EQ 0 OR n_params() LT 2) then $
    message, 'INLOGLAM, OBJFLUX and NEWLOGLAM are all required'
   if (n_elements(objflux) NE npix) then $
    message, 'Dimensions of INLOGLAM and OBJFLUX do not agree'
   if (keyword_set(objivar)) then $
    if (n_elements(objivar) NE npix) then $
     message, 'Dimensions of INLOGLAM and OBJIVAR do not agree'
   if (keyword_set(finalmask)) then $
    if (n_elements(finalmask) NE npix) then $
     message, 'Dimensions of INLOGLAM and FINALMASK do not agree'
   if (keyword_set(indisp)) then $
    if (n_elements(indisp) NE npix) then $
     message, 'Dimensions of INLOGLAM and INDISP do not agree'

   ;----------
   ; Set defaults

   if (NOT keyword_set(binsz)) then begin
      if (npix EQ 1) then binsz = 1 $
       else binsz = inloglam[1] - inloglam[0]
   endif
   if (NOT keyword_set(nord)) then nord = 3
   if (NOT keyword_set(bkptbin)) then bkptbin = 1.2 * binsz
   if (NOT keyword_set(maxsep)) then maxsep = 2.0 * binsz

   ndim = size(inloglam, /n_dimen)
   if (ndim EQ 1) then $
    specnum = fltarr(n_elements(inloglam)) $ ; Set specnum=0 for all elements
   else $
    specnum = djs_laxisgen( size(inloglam,/dimens), iaxis=1)

   ; Use fullcombmask for modifying the pixel masks in the original input files.
   fullcombmask = bytarr(npix)

   newflux = fltarr(nfinalpix)
   if (arg_present(newivar)) then newivar = fltarr(nfinalpix)
   if (arg_present(newdisp)) then newdisp = fltarr(nfinalpix)

   if (keyword_set(objivar)) then begin
      nonzero = where(objivar GT 0.0, ngood)
   endif else begin
      nonzero = lindgen(npix)
      ngood = npix
   endelse

   if (ngood EQ 0) then begin

      splog, 'No good points'
      andmask = pixelmask_bits('NODATA')
      ormask = pixelmask_bits('NODATA')
      return

   endif else begin

      ; Now let's break sorted wavelengths into groups where
      ; pixel separations are larger than maxsep

      isort = nonzero[sort(inloglam[nonzero])]
      wavesort = inloglam[isort]

      padwave = [min(wavesort) - 2.0*maxsep, wavesort, $
       max(wavesort) + 2.0*maxsep]

      ig1 = where(padwave[1:ngood] - padwave[0:ngood-1] GT maxsep, nstart)
      ig2 = where(padwave[2:ngood+1] - padwave[1:ngood] GT maxsep, nend)
      if (nstart NE nend) then $
       message, 'ABORT: Grouping tricks did not work!'

      for igrp=0, nstart-1 do begin

         ss = isort[ig1[igrp] : ig2[igrp]]
         bkpt = 0

         if (keyword_set(objivar)) then begin
            fullbkpt = slatec_splinefit(inloglam[ss], objflux[ss], coeff, $
             nord=nord, eachgroup=1, maxiter=20, upper=3.0, lower=3.0, $
             bkspace=bkptbin, bkpt=bkpt, invvar=objivar[ss], mask=bmask, $
             /silent)
         endif else begin
            fullbkpt = slatec_splinefit(inloglam[ss], objflux[ss], coeff, $
             nord=nord, eachgroup=1, maxiter=20, upper=3.0, lower=3.0, $
             bkspace=bkptbin, bkpt=bkpt, mask=bmask, $
             /silent)
         endelse

         inside = where(newloglam GE min(bkpt) $
          AND newloglam LE max(bkpt), numinside)

         if (numinside EQ 0) then begin
            splog,'WARNING: No wavelengths inside breakpoints'
         endif else if (total(abs(coeff)) EQ 0.0) then begin
            splog,'WARNING: All B-spline coefficients have been set to zero!'
         endif else begin         

            newflux[inside] = slatec_bvalu(newloglam[inside], fullbkpt, coeff)

            splog, 'Masked ', fix(total(1-bmask)), ' of', $
             n_elements(bmask), ' pixels'

            ;----------
            ; Determine which pixels should be masked based upon the spline fit.
            ; Set the COMBINEREJ bit.

            ireplace = where(bmask EQ 0)

            if (ireplace[0] NE -1) then begin 
               ; The following would replace the original flux values
               ; of masked pixels with b-spline evaluations.
;               objflux[ss[ireplace]] = $
;                slatec_bvalu(inloglam[ss[ireplace]],fullbkpt,coeff)

               if (keyword_set(finalmask)) then $
                finalmask[ss[ireplace]] = finalmask[ss[ireplace]] OR $
                 pixelmask_bits('COMBINEREJ')
            endif

         endelse
         fullcombmask[ss] = bmask

      endfor

      ;---------------------------------------------------------------------
      ; Combine inverse variance and pixel masks.

      if (arg_present(andmask)) then $
       andmask = lonarr(nfinalpix) - 1 ; Start with all bits set in AND-mask.
      if (arg_present(ormask)) then $
       ormask = lonarr(nfinalpix)

      for j=0, max(specnum) do begin
         these = where(specnum EQ j)

         if (these[0] NE -1) then begin

            inbetween = where(newloglam GE min(inloglam[these]) AND $
                              newloglam LE max(inloglam[these]))
            if (inbetween[0] NE -1) then begin

               if (arg_present(newivar)) then begin
                  ; Conserve inverse variance by doing a linear interpolation
                  ; on that quantity.

                  result = interpol(objivar[these] * fullcombmask[these], $
                   inloglam[these], newloglam[inbetween])

                  ; Grow the fullcombmask below to reject any new sampling
                  ; containing even a partial masked pixel.

                  smask = interpol(float(fullcombmask[these]), $
                   inloglam[these], newloglam[inbetween])
                  ibad = where(smask LT 1.0)
                  if (ibad[0] NE -1) then result[ibad] = 0

                  newivar[inbetween] = newivar[inbetween] + result
               endif

               lowside = fix((inloglam[these]-newloglam[0])/binsz)
               highside = lowside + 1

               if (arg_present(andmask) AND keyword_set(finalmask)) then begin
                  andmask[lowside] = andmask[lowside] AND finalmask[these]
                  andmask[highside] = andmask[highside] AND finalmask[these]
               endif

               if (arg_present(ormask) AND keyword_set(finalmask)) then begin
                  ormask[lowside] = ormask[lowside] OR finalmask[these]
                  ormask[highside] = ormask[highside] OR finalmask[these]
               endif

               if (arg_present(newdisp)) then begin
                  ; Combine the dispersions in the dumbest way possible
                  newdisp[inbetween] = interpol(indisp[these], $
                   inloglam[these], newloglam[inbetween])
               endif
            endif

         endif
      endfor
      splog, 'Medians:', string(format='(f7.2)', $
       djs_median(objflux[these]))

   endelse

   ;----------
   ; Replace NaN's in combined spectra; this should really never happen

   if (keyword_set(newivar)) then $
    inff = where(finite(newflux) EQ 0 OR finite(newivar) EQ 0) $
   else $
    inff = where(finite(newflux) EQ 0)

   if (inff[0] NE -1) then begin
      splog, 'WARNING: NaNs in combined spectra ', N_elements(inff)
      newflux[inff] = 0.0
      if (keyword_set(newivar)) then newivar[inff] = 0.0
   endif

   if (keyword_set(newivar)) then begin
      ;----------
      ; Interpolate over masked pixels, just for aesthetic purposes

       newflux = djs_maskinterp(newflux, newivar EQ 0, /const)

      ;----------
      ; Set the NODATA mask bit wherever there is no good data

      ibad = where(newivar EQ 0)
      if (ibad[0] NE -1) then begin
         if (keyword_set(andmask)) then $
          andmask[ibad] = andmask[ibad] AND pixelmask_bits('NODATA')
         if (keyword_set(ormask)) then $
          ormask[ibad] = ormask[ibad] AND pixelmask_bits('NODATA')
      endif
   endif

   ;----------
   ; Replace values of -1 in the AND mask with 0's

   if (keyword_set(andmask)) then $
    andmask = andmask * (andmask NE -1)

   return
end
;------------------------------------------------------------------------------
