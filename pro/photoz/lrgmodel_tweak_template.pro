;+
; NAME:
;   lrgmodel_tweak_template
;
; PURPOSE:
;   Compute the corrections to the input LRG template used for LRG_PHOTOZ().
;
; CALLING SEQUENCE:
;   lrgmodel_tweak_template, pflux, pflux_ivar, zz, [ weights=, $
;    maxiter=, ageguess=, metalguess=, agerange=, metalrange=, coeff=, $
;     _EXTRA= ]
;
; INPUTS:
;   pflux          - Object fluxes in the 5 SDSS filters [5,NOBJ]
;   pflux_ivar     - Inverse variances for FLUX [5,NOBJ]
;   zz             - Spectroscopic redshifts [NOBJ]
;
; OPTIONAL INPUTS:
;   weights        - Weights for each object in the minimization
;                    (this is multiplied by the chi^2 of each point)
;   maxiter        - Maximum number of iterations; default to 40
;   ageguess       - Starting guess for AGEBURST; default to 2.5 Gyr
;   metalguess     - Starting guess for ZMETAL; default to Z=0.2
;   agerange       - Range of possible values for AGEBURST; default
;                    to [0.0, 6.0] Gyr
;   metalrange     - Range of possible values for ZMETAL; default
;                    to [0.008, 0.05]
;   _EXTRA         - Keywords for LRGMODEL_PHOTOZ(), such as ABCORRECT,
;                    EXTINCTION, ABFUDGE, ADDERR, FILTERLIST
;
; OUTPUTS:
;   coeff          - Best-fit coefficients for tweaking input LRG template
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine minimizes the parameters used for the Bruzual-Charlot
;   models in computing photometric redshifts.  There are two parameters
;   on which we minimize:
;      AGEBURST - age of the Universe at the burst [Gyr]
;      ZMETAL   - metallicity of the burst
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   lrgmodel_photoz()
;   mpfit()
;
; INTERNAL SUPPORT ROUTINES:
;  lrgmodel_tweak_fn()
;
; REVISION HISTORY:
;   10-Dec-2003  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
forward_function mpfit, lrgmodel_tweak_fn

;------------------------------------------------------------------------------
; Return a vector of all the chi values
function lrgmodel_tweak_fn, coeff

   common com_lrgmodel_tweak_fluxes, pflux, pflux_ivar, zz, $
    sqweights, KeywordsForPhotoz

   ; Compute the best-fit redshift for each object
   zfit = lrgmodel_photoz(pflux, pflux_ivar, $
    ageburst=coeff[0], zmetal=coeff[1], chi2=chi2, $
    _EXTRA=KeywordsForPhotoz)

   ; Return the redshift errors as the "chi" values
   chivec = zfit - zz
   if (keyword_set(sqweights)) then chivec = chivec * sqweights

   return, chivec
end
;------------------------------------------------------------------------------
pro lrgmodel_tweak_template, pflux1, pflux_ivar1, zz1, weights=weights, $
 maxiter=maxiter, ageguess=ageguess1, metalguess=metalguess1, $
 agerange=agerange1, metalrange=metalrange1, _EXTRA=EXTRA

   common com_lrgmodel_tweak_fluxes, pflux, pflux_ivar, zz, $
    sqweights, KeywordsForPhotoz

   ; Set defaults
   if (NOT keyword_set(maxiter)) then maxiter = 40
   if (keyword_set(ageguess1)) then ageguess = ageguess1 $
    else ageguess = 2.5d0
   if (keyword_set(metalguess1)) then metalguess = metalguess1 $
    else metalguess = 0.02
   if (keyword_set(agerange1)) then agerange = agerange1 $
    else agerange = [0.0, 6.0]
   if (keyword_set(metalrange1)) then metalrange = metalrange1 $
    else metalrange = [0.008, 0.05]
   if (keyword_set(EXTRA)) then KeywordsForPhotoz = EXTRA $
    else KeywordsForPhotoz = 0

   ; Set variables in common blocks
   zz = zz1
   if (keyword_set(weights)) then sqweights = sqrt(weights)

   ; Discard any objects where the baseline photo-z is discrepent by
   ; more than 0.10, and discard any low-redshift objects with z > 0.10.
   pflux = pflux1
   pflux_ivar = pflux_ivar1
   zfit = lrgmodel_photoz(pflux, pflux_ivar, _EXTRA=KeywordsForPhotoz)
;   ibad = where(abs(zfit - zz) GT 0.10 OR zz LT 0.10, nbad)
; Only discard the low-redshift objects...
   ibad = where(zz LT 0.10, nbad)
   if (nbad GT 0) then begin
      print, 'Discard ', nbad, ' objects with low and/or discrepent photo-z'
      pflux[*,ibad] = 0
      pflux_Ivar[*,ibad] = 0
   endif

   ; Call MPFIT to iterate on the solution for the template
   parinfo = {value: 0.D, fixed: 0, limited: [0b,0b], limits: [0.d0,0.d0]}
   parinfo = replicate(parinfo, 2)
   parinfo.value = [ageguess, metalguess]
   parinfo[0].limited = [1b, 1b]
   parinfo[0].limits = agerange
   parinfo[1].limits = metalrange

   ftol = 1d-20
   gtol = 1d-20
   xtol = 1d-20
   coeff = mpfit('lrgmodel_tweak_fn', parinfo=parinfo, perror=perror, $
    maxiter=maxiter, ftol=ftol, gtol=gtol, xtol=xtol, $
    niter=niter, status=status)

   print, 'STATUS = ', status
   print, 'Best-fit AGEBURST = ', coeff[0]
   print, 'Best-fit ZMETAL = ', coeff[1]
   print, 'Errors = = ', perror

   return
end
;------------------------------------------------------------------------------
