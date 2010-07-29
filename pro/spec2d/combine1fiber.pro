;+
; NAME:
;   combine1fiber
;
; PURPOSE:
;   Combine several spectra of the same object, or resample a single spectrum.
;
; CALLING SEQUENCE:
;   combine1fiber, inloglam, objflux, [ objivar, finalmask=, indisp=, skyflux=,$
;    newloglam=, newflux=, newivar=, andmask=, ormask=, newdisp=, newsky=, $
;    nord=, binsz=, bkptbin=, maxsep=, _EXTRA=KeywordsForReject, /verbose ]
;
; INPUTS:
;   inloglam       - Wavelengths in log10-Angstroms [NPIX,NSPEC]
;   objflux        - Flux [NPIX,NSPEC]
;
; REQUIRED KEYWORDS:
;   newloglam      - Wavelengths for output evaluation, also in log10-Angstroms
;                    [NNEWPIX]
;
; OPTIONAL INPUTS:
;   objivar        - Inverse variance [NPIX,NSPEC]
;   finalmask      - Pixel mask [NPIX,NSPEC]
;   indisp         - Dispersion values [NPIX,NSPEC]
;   skyflux        - Sky flux vectors [NPIX,NSPEC]
;   binsz          - Bin separation for INLOGLAM; if not set, then default
;                    to INLOGLAM[1]-INLOGLAM[0].
;   nord           - Order of spline fit; default to 3.
;   bkptbin        - Break point binning; default to 1.2 * BINSZ.
;   maxsep         - Maximum separation between input wavelengths.  The spline
;                    fit is split into pieces, with the breaks wherever this
;                    spacing is exceeded.  Default to 2.0 * BINSZ.
;   _EXTRA         - Keywords for DJS_REJECT().
;   verbose        - If set, then output messages about bad break points and
;                    masked data points.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   finalmask      - Modified from its input by setting the COMBINEREJ bit
;                    for deviant pixels in individual spectra that have
;                    been rejected.
;   objivar        - Modified from itts input by setting to zero wherever
;                    the COMBINEREJ bit has been set in FINALMASK.
;   newflux        - Resampled flux [NNEWPIX].
;   newivar        - Resampled inverse variance [NNEWPIX].
;   andmask        - Resampled mask. For each mask bit, set that bit only if
;                    every input spectrum at this wavelength has that bit set
;                    (e.g., this is a logical AND) [NNEWPIX].
;   ormask         - Resampled mask. For each mask bit, set that bit if any
;                    of the input spectra at this wavelength has that bit set
;                    (e.g., this is a logical OR) [NNEWPIX].
;   newdisp        - Resampled dispersion values [NNEWPIX].
;   newsky         - Resampled sky flux [NNEWPIX].
;
; COMMENTS:
;   One can pass this routine a single spectrum to be fit by a spline and
;   re-sampled, in which case all the inputs (such as FLUX) are 1-dimensional
;   arrays.  Or, one can pass it several spectra, in which case these inputs
;   are 2-dimensional arrays.
;
;   There's also some code in here to grow masked regions by another pixel
;   on either end if the region is more than 3 pixels wide.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   aesthetics()
;   bspline_iterfit()
;   bspline_valu()
;   djs_laxisgen()
;   djs_median()
;   pixelmask_bits()
;   splog
;
; REVISION HISTORY:
;   02-Jan-2000  Written by D. Schlegel; modified from COMBINE2DOUT
;-
;------------------------------------------------------------------------------
pro combine1fiber, inloglam, objflux, objivar, $
    finalmask=finalmask, indisp=indisp, skyflux=skyflux, $
    newloglam=newloglam, newflux=newflux, newivar=newivar, $
    andmask=andmask, ormask=ormask, newdisp=newdisp, newsky=newsky, $
    nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
    _EXTRA=KeywordsForReject, verbose=verbose

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

  EPS = (machar()).eps ; Use this to avoid some round-off errors

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
  newmask = lonarr(nfinalpix)
  if (arg_present(newivar)) then newivar = fltarr(nfinalpix)
  if (arg_present(newdisp)) then newdisp = fltarr(nfinalpix)
  if (arg_present(newsky)) then newsky = fltarr(nfinalpix)
  if (arg_present(newdisp) OR arg_present(newsky)) then $
    newdispweight = fltarr(nfinalpix)

  if (keyword_set(objivar)) then begin
    nonzero = where(objivar GT 0.0, ngood)
  endif else begin
    nonzero = lindgen(npix)
    ngood = npix
  endelse

  ; Create the ORMASK if we need ANDMASK, since we set the NODATA
  ; bad pixels in the ORMASK first.
  if (arg_present(andmask)) then $
    andmask = lonarr(nfinalpix)
  if (arg_present(ormask) OR arg_present(andmask)) $
    then ormask = lonarr(nfinalpix)

  if (ngood EQ 0) then begin

    ; In this case of no good points, set the NODATA bit everywhere.
    ; Also, if NOPLUG is set in the first input bit-mask, assume it
    ; should be set everywhere in the output bit masks.
    ; No other bits are set.

    if (keyword_set(verbose)) then splog, 'No good points'
    bitval = pixelmask_bits('NODATA')
    if (keyword_set(finalmask)) then bitval = bitval OR $
      (pixelmask_bits('NOPLUG') * (finalmask[0] AND pixelmask_bits('NOPLUG')))
    if (keyword_set(andmask)) then andmask = andmask OR bitval
    if (keyword_set(ormask)) then ormask = ormask OR bitval

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

    for igrp=0L, nstart-1 do begin

      ss = isort[ig1[igrp] : ig2[igrp]]
      bkpt = 0
      bmask = 0

      if (n_elements(ss) GT 2) then begin
        if (keyword_set(objivar)) then $
          sset = bspline_iterfit(inloglam[ss], objflux[ss], $
          nord=nord, /groupbadpix, requiren=1, $
          bkspace=bkptbin, bkpt=bkpt, invvar=objivar[ss], outmask=bmask, $
          _EXTRA=KeywordsForReject, /silent) $
        else $
          sset = bspline_iterfit(inloglam[ss], objflux[ss], $
          nord=nord, /groupbadpix, requiren=1, $
          bkspace=bkptbin, bkpt=bkpt, outmask=bmask, $
          _EXTRA=KeywordsForReject, /silent)
        if (total(abs(sset.coeff)) EQ 0.0) then begin
          sset = 0
          if (keyword_set(verbose)) then $
            splog, 'WARNING: All B-spline coefficients have been set to zero!'
        endif
      endif else begin
        bmask = bytarr(n_elements(ss)) ; All set to zero (rejected)
        sset = 0
        if (keyword_set(verbose)) then $
          splog,'WARNING: All B-spline coefficients have been set to zero!'
      endelse

      inside = where(newloglam GE min(inloglam[ss])-EPS $
        AND newloglam LE max(inloglam[ss])+EPS, numinside)

      ; It is possible for NUMINSIDE to be zero, if the input data points
      ; span an extremely small wavelength range, within which there are
      ; no output wavelengths.
      if (keyword_set(sset) AND numinside GT 0) then begin
        newflux[inside] = bspline_valu(newloglam[inside], sset, $
          mask=bvalumask)
        goodvalu = where(bvalumask)
        if goodvalu[0] NE -1 then newmask[inside[goodvalu]] = 1

        if (keyword_set(verbose)) then $
          splog, 'Masked ', long(total(1-bmask)), ' of', $
          n_elements(bmask), ' pixels'

        ;----------
        ; Determine which pixels should be masked based upon the spline fit.
        ; Set the COMBINEREJ bit.

        ireplace = where(bmask EQ 0)

        if (ireplace[0] NE -1) then begin
          ; The following would replace the original flux values
          ; of masked pixels with b-spline evaluations.
          ;               objflux[ss[ireplace]] = $
          ;                bspline_valu(inloglam[ss[ireplace]], sset)

          ; Set the inverse variance of these pixels to zero.
          if (keyword_set(objivar)) then $
            objivar[ss[ireplace]] = 0.0

          if (keyword_set(finalmask)) then $
            finalmask[ss[ireplace]] = finalmask[ss[ireplace]] OR $
            pixelmask_bits('COMBINEREJ')
        endif

      endif
      fullcombmask[ss] = bmask

    endfor

    ;---------------------------------------------------------------------
    ; Combine inverse variance and pixel masks.

    ; Start with all bits set in AND-mask.
    if (arg_present(andmask)) then andmask[*] = -1L

    for j=0L, max(specnum) do begin
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
            ibad = where(smask LT 1.0 - EPS)
            if (ibad[0] NE -1) then result[ibad] = 0

            newivar[inbetween] = newivar[inbetween] + $
              result * newmask[inbetween]

          endif

          lowside = floor((inloglam[these]-newloglam[0])/binsz)
          highside = lowside + 1L

          if (arg_present(andmask) AND keyword_set(finalmask)) then begin
            andmask[lowside] = andmask[lowside] AND finalmask[these]
            andmask[highside] = andmask[highside] AND finalmask[these]
          endif

          if (arg_present(ormask) AND keyword_set(finalmask)) then begin
            ormask[lowside] = ormask[lowside] OR finalmask[these]
            ormask[highside] = ormask[highside] OR finalmask[these]
          endif

          ; Combine the dispersions + skies in the dumbest way possible
          if (arg_present(newdisp) OR arg_present(newsky)) then begin
            newdispweight[inbetween] = newdispweight[inbetween] + result
            if (arg_present(newdisp)) then $
              newdisp[inbetween] = newdisp[inbetween] + $
              interpol(indisp[these], inloglam[these], $
              newloglam[inbetween]) * result
            if (arg_present(newsky)) then $
              newsky[inbetween] = newsky[inbetween] + $
              interpol(skyflux[these], inloglam[these], $
              newloglam[inbetween]) * result
          endif
        endif

      endif
    endfor

    if (arg_present(newdisp)) then $
      newdisp  = newdisp / (newdispweight + (newdispweight EQ 0))
    if (arg_present(newsky)) then $
      newsky  = newsky / (newdispweight + (newdispweight EQ 0))

  endelse

  ;----------
  ; Grow regions where 3 or more pixels are rejected together ???

  if (keyword_set(newivar)) then begin
    badregion = where(smooth(newivar,3) EQ 0,nbad)
    if (keyword_set(verbose)) then $
      splog, 'WARNING: Growing bad pixel region, ',nbad,' pixels found.'
    if badregion[0] NE -1 then begin
      newivar[(badregion - 2) > 0] = 0.0
      newivar[(badregion + 2) < (nfinalpix - 1)] = 0.0
    endif
  endif

  ;----------
  ; Replace NaN's in combined spectra; this should really never happen

  if (keyword_set(newivar)) then $
    inff = where(finite(newflux) EQ 0 OR finite(newivar) EQ 0) $
  else $
    inff = where(finite(newflux) EQ 0)

  if (inff[0] NE -1) then begin
    splog, 'WARNING: ',N_elements(inff),' NaNs in combined spectra.'
    newflux[inff] = 0.0
    if (keyword_set(newivar)) then newivar[inff] = 0.0
  endif

  if (keyword_set(newivar)) then begin
    ;----------
    ; Fix up the masked pixels.
    ;
    newflux = aesthetics(newflux,newivar) ;,method='mean')
    goodpts = where(newivar GT 0)
    if goodpts[0] NE -1 then begin
      minglam = min(newloglam[goodpts])
      maxglam = max(newloglam[goodpts])

      ibad = where(newloglam LT minglam OR newloglam GT maxglam)
      if (ibad[0] NE -1) then begin
        if (keyword_set(ormask)) then $
          ormask[ibad] = ormask[ibad] OR pixelmask_bits('NODATA')
      endif
    endif

  ;----------
  ; Set the NODATA mask bit wherever there is no good data
  ;
  ;    ibad = where(newivar EQ 0)
  ;    if (ibad[0] NE -1) then begin
  ;      if (keyword_set(andmask)) then $
  ;        andmask[ibad] = andmask[ibad] OR pixelmask_bits('NODATA')
  ;      if (keyword_set(ormask)) then $
  ;        ormask[ibad] = ormask[ibad] OR pixelmask_bits('NODATA')
  ;    endif

  endif

  ;----------
  ; Replace values of -1 in the AND mask with 0's

  if (keyword_set(andmask)) then $
    andmask = andmask * (andmask NE -1L)

  ;----------
  ; Copy the NODATA bad pixels from the ORMASK to the ANDMASK

  if (keyword_set(andmask)) then $
    andmask = andmask OR (ormask AND pixelmask_bits('NODATA'))

  return
end
;------------------------------------------------------------------------------
