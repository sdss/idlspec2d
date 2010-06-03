;+
; NAME:
;   lrgmodel_train
;
; PURPOSE:
;   Calling script for LRGMODEL_TWEAK_TEMPLATE for training the Bruzual-
;   Charlot models for use with photo-z's.
;
; CALLING SEQUENCE:
;   lrgmodel_train, [ public=, /regenerate, twodf=, zrange=, dzweight=, $
;    /absdev, _EXTRA=]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   public         - If set, then use the requested version of the spAll
;                    file rather than all objects; for example, 'EDR', 'DR1',
;                    'DR2', or simply /PUBLIC to select all public data.
;   regenerate     - If set, then regenerate the trimmed spAll file.
;   twodf          - Suffix for which 2dF files to use; default to the two
;                    set of files specified by ['2003A', '2003B']
;   zrange         - Trim to objects in this redshift range; default [0.1,1.0].
;   dzweight       - Re-weight the galaxies such that all objects within each
;                    redshift interval DZWEIGHT share unit weight; default
;                    to 0.02; set to 0 to disable this re-weighting.
;   absdev         - If set, then minimize the absolute value of the redshift
;                    deviations, rather than the more conventional square of
;                    those deviations.  This should be more robust.
;   _EXTRA         - Keywords for LRGMODEL_PHOTOZ(), such as ABCORRECT,
;                    EXTINCTION, ABFUDGE, ADDERR, FILTERLIST.
;
; OUTPUTS:
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
;   This routine first tries a coarse grid of AGEBURST,ZMETAL.
;   The best-fit from that grid is used for the initial guess
;   that's passed to MPFIT() for a full minimization.
;
; BUGS:
;   Note that the reported chi^2 values are really just sums of the
;   squares of redshift errors (multiplied by the weights).
;
; DATA FILES:
;   Both SDSS redshifts and 2dF redshifts are used, reading the files:
;      $BOSS_SPECTRO_REDUX/spAll.fits
;      $BOSS_SPECTRO_REDUX/spAll-$PUBLIC.fits (if PUBLIC is specified)
;      $SDSS_2DF/data/calibObj-2003A.fits
;      $SDSS_2DF/data/calibObj-2003B.fits
;      $SDSS_2DF/data/catalogue2003A.fits
;      $SDSS_2DF/data/catalogue2003B.fits
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   hogg_mrdfits()
;   lrgmodel_photoz()
;   lrgmodel_tweak_template()
;   mrdfits()
;   splog
;
; REVISION HISTORY:
;   18-Dec-2003  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro lrgmodel_train, public=public, regenerate=regenerate, $
 twodf=twodf, zrange=zrange, dzweight=dzweight, absdev=absdev, $
 _EXTRA=KeywordsForPhotoz

   if (n_elements(twodf) EQ 0) then twodf = ['2003A', '2003B']
   if (n_elements(zrange) EQ 0) then zrange = [0.10, 0.60]
   if (n_elements(dzweight) EQ 0) then dzweight = 0.02

   splog, file='lrgmodel_train.log'

   ;----------
   ; Read the data files

   ; Read the data from spAll
   spall = lrgmodel_read_spall(public=public, regenerate=regenerate)

   ; Read the data from SDSS-2dF
   if (keyword_set(twodf)) then begin
      for i=0, n_elements(twodf)-1 do $
       spall = lrgmodel_append_twodf(spall, twodf[i])
   endif

   ; Trim to only objects in the redshift range of interest
   if (keyword_set(zrange)) then begin
      itrim = where(spall.z GE zrange[0] AND spall.z LE zrange[1], ntrim)
      splog, 'Trimming from ', n_elements(spall), ' to ', ntrim, ' objects'
      spall = spall[itrim]
   endif

   ;----------
   ; Select weights for each galaxy that uniformly weights in
   ; redshift bins

   if (keyword_set(dzweight)) then begin
      nobj = n_elements(spall)
      weights = fltarr(nobj)
      for z1=0.0, max(spall.z), dzweight do begin
         ii = where(spall.z GE z1 AND  spall.z LT z1+dzweight, ct)
         if (ct GT 0) then weights[ii] = 1./ct
      endfor
      ; Set the mean weight to unity
      weights = weights / mean(weights)
   endif

   ;----------
   ; Start with a coarse grid in AGEBURST,ZMETAL
   ; in order to get a good starting point

   ageburst = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
   zmetal = [0.008, 0.014, 0.02, 0.025, 0.03, 0.04]
;ageburst = 4.0 ; Test ???
;zmetal = 0.03 ; Test ???

   ageguess = 0
   metalguess = 0
   chi2best = 0
   for iage=0, n_elements(ageburst)-1 do begin
      for imetal=0, n_elements(zmetal)-1 do begin
         zfit = lrgmodel_photoz(spall.modelflux, spall.modelflux_ivar, $
          ageburst=ageburst[iage], zmetal=zmetal[imetal], $
          _EXTRA=KeywordsForPhotoz)
         chivec = (zfit - spall.z) / 0.01 ; This is an arbitrary normalization
         if (keyword_set(weights)) then chivec = chivec * sqrt(weights)
         if (keyword_set(absdev)) then thischi2 = total(abs(chivec)) $
          else thischi2 = total(chivec^2)
         splog, 'For AGE=', ageburst[iage], ' ZMETAL=', zmetal[imetal], $
          ' chi2=', thischi2
         if (chi2best EQ 0 OR thischi2 LT chi2best) then begin
            ageguess = ageburst[iage]
            metalguess = zmetal[imetal]
            chi2best = thischi2
         endif
      endfor
   endfor

   splog, 'Initial guess AGE=', ageguess, ' ZMETAL=', metalguess

   ;----------
   ; Now do the full minimization problem

   lrgmodel_tweak_template, spall.modelflux, spall.modelflux_ivar, $
    spall.z, weights=weights, $
    ageguess=ageguess, metalguess=metalguess, absdev=absdev, $
    coeff=coeff, _EXTRA=KeywordsForPhotoz

   splog, 'Best guess = ', coeff
   splog, /close

   return
end
;------------------------------------------------------------------------------
