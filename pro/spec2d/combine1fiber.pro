;+
; NAME:
;   combine1fiber
;
; PURPOSE:
;   Combine several spectra of the same object, or resample a spectra.
;
; CALLING SEQUENCE:
;   combine1fiber, fullwave, fullspec, fullivar, fullpixelmask, fulldisp, $
;    finalwave, bestflux, bestivar, andmask, ormask, bestdisp, $
;    nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep
;
; INPUTS:
;   fullwave       - Wavelengths in log10-Angstroms
;   fullspec       - Flux
;   fullivar       - Inverse variance
;   fullpixelmask  - Pixel mask
;   fulldisp       - Dispersion values
;   finalwave      - Wavelengths for output evaluation, also in log10-Angstroms
;
; REQUIRED KEYWORDS:
;   binsz          - Bin separation for FULLWAVE
;
; OPTIONAL KEYWORDS:
;   nord           -
;   bkptbin        -
;   maxsep         -
;
; OUTPUTS:
;   bestflux       -
;   bestivar       -
;   andmask        -
;   ormask         -
;   bestdisp       -
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
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
pro combine1fiber, fullwave, fullspec, fullivar, fullpixelmask, fulldisp, $
 finalwave, bestflux, bestivar, andmask, ormask, bestdisp, $
 nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep

; ??? Return fullcombmask for modifying the masks in the original input files.

   if (NOT keyword_set(nord)) then nord = 3
   if (NOT keyword_set(maxsep)) then maxsep = 2.0 * binsz
   if (NOT keyword_set(bkptbin)) then bkptbin = 1.2 * binsz

   specnum = djs_laxisgen( size(fullwave,/dimens), iaxis=1)

   fullcombmask = bytarr(n_elements(fullspec))

   nfinalpix = N_elements(finalwave)
   bestflux = fltarr(nfinalpix)
   bestivar = bestflux*0.0
   bestdisp = bestflux*0.0

   nonzero = where(fullivar GT 0.0, ngood)

   if (ngood EQ 0) then begin

      splog, 'No good points'
      andmask = pixelmask_bits('NODATA')
      ormask = pixelmask_bits('NODATA')
      return

   endif else begin

      ; Now let's break sorted wavelengths into groups where
      ; pixel separations are larger than maxsep

      isort = nonzero[sort(fullwave[nonzero])]
      wavesort = fullwave[isort]

      padwave = [min(wavesort) - 2.0*maxsep, wavesort, $
       max(wavesort) + 2.0*maxsep]

      ig1 = where(padwave[1:ngood] - padwave[0:ngood-1] GT maxsep, nstart)
      ig2 = where(padwave[2:ngood+1] - padwave[1:ngood] GT maxsep, nend)
      if (nstart NE nend) then $
       message, 'ABORT: Grouping tricks did not work!'

      for igrp=0, nstart-1 do begin

         ss = isort[ig1[igrp] : ig2[igrp]]
         bkpt = 0

         fullbkpt = slatec_splinefit(fullwave[ss], fullspec[ss], coeff, $
          nord=nord, eachgroup=1, maxiter=20, upper=3.0, lower=3.0, $
          bkspace=bkptbin, bkpt=bkpt, invvar=fullivar[ss], mask=bmask, /silent)

         inside = where(finalwave GE min(bkpt) $
          AND finalwave LE max(bkpt), numinside)

         if (numinside EQ 0) then begin
            splog,'WARNING: No wavelengths inside breakpoints'
         endif else if (total(abs(coeff)) EQ 0.0) then begin
            splog,'WARNING: All B-spline coefficients have been set to zero!'
         endif else begin         

            bestflux[inside] = slatec_bvalu(finalwave[inside],fullbkpt,coeff)

            splog, 'Masked ', fix(total(1-bmask)), ' of', $
             n_elements(bmask), ' pixels'

            ;-----------------------------------------------------------------
            ;  Here replace original flux values of masked pixels with b-spline
            ;  evaluations.

            ireplace = where(bmask EQ 0)

            if (ireplace[0] NE -1) then begin 
               fullspec[ss[ireplace]] = $
                slatec_bvalu(fullwave[ss[ireplace]],fullbkpt,coeff)

               fullpixelmask[ss[ireplace]] = fullpixelmask[ss[ireplace]] OR $
                pixelmask_bits('COMBINEREJ')
            endif

         endelse
         fullcombmask[ss] = bmask

      endfor

      ;---------------------------------------------------------------------
      ; Combine inverse variance and pixel masks.

      andmask = lonarr(nfinalpix) - 1 ; Start with all bits set in AND-mask.
      ormask = lonarr(nfinalpix)

      for j=0, max(specnum) do begin
         these = where(specnum EQ j)
         if (these[0] NE -1) then begin

            inbetween = where(finalwave GE min(fullwave[these]) AND $
                              finalwave LE max(fullwave[these]))
            if (inbetween[0] NE -1) then begin

               ; Conserve inverse variance by doing a linear interpolation
               ; on that quantity.

               result = interpol(fullivar[these] * fullcombmask[these], $
                fullwave[these], finalwave[inbetween])

               ; Grow the fullcombmask below to reject any new sampling
               ; containing even a partial masked pixel.

               smask = interpol(float(fullcombmask[these]), $
                fullwave[these], finalwave[inbetween])
               ibad = where(smask LT 1.0)
               if (ibad[0] NE -1) then result[ibad] = 0

               bestivar[inbetween] = bestivar[inbetween] + result

               lowside = fix((fullwave[these]-finalwave[0])/binsz)
               highside = lowside + 1
               andmask[lowside]  = andmask[lowside] AND fullpixelmask[these]
               andmask[highside] = andmask[highside] AND fullpixelmask[these]
               ormask[lowside]   = ormask[lowside] OR fullpixelmask[these]
               ormask[highside]  = ormask[highside] OR fullpixelmask[these]

               ; Combine the dispersions in the dumbest way possible

               bestdisp[inbetween] = interpol(fulldisp[these], $
                fullwave[these], finalwave[inbetween])

            endif

         endif
      endfor
      splog, 'Medians:', string(format='(f7.2)', $
       djs_median(fullspec[these]))

   endelse

   ;----------
   ; Replace NaN's in combined spectra; this should really never happen

   inff = where(finite(bestflux) EQ 0 OR finite(bestivar) EQ 0)
   if (inff[0] NE -1) then begin
      splog, 'WARNING: NaNs in combined spectra ', N_elements(inff)
      bestflux[inff] = 0.0
      bestivar[inff] = 0.0
   endif

   ;----------
   ; Interpolate over masked pixels, just for aesthetic purposes

   bestflux = djs_maskinterp(bestflux, bestivar EQ 0, /const)

   ;----------
   ; Set the NODATA mask bit wherever there is no good data

   ibad = where(bestivar EQ 0)
   if (ibad[0] NE -1) then begin
      andmask[ibad] = andmask[ibad] AND pixelmask_bits('NODATA')
      ormask[ibad] = ormask[ibad] AND pixelmask_bits('NODATA')
   endif

   ;----------
   ; Replace values of -1 in the AND mask with 0's

   andmask = andmask * (andmask NE -1)

   return
end
;------------------------------------------------------------------------------
