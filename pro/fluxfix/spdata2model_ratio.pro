;+
; NAME:
;   spdata2model_ratio
;
; PURPOSE:
;   Construct flux correction vectors by ratioing the observed spectra of
;   standard stars to Kurucz models of the appropriate spectral type and
;   smoothing the result.
;
; CALLING SEQUENCE:
;    corvector = spdata2model_ratio(loglam, stdflux, stdivar, stdmask, 
;      stdinfo, corvivar=, cormed=, /norm)
;
; INPUTS:
;   loglam  -- wavelength array in log10(Angstroms) [npix]
;   stdflux -- array of standard star fluxes [npix, nstd]
;   stdivar -- inverse variance of standard star fluxes [npix, nstd]
;   stdmask -- ormask of standard star spectra [npix, nstd] 
;              (used to mask sky residuals)
;   stdinfo -- structure containing information about which Kurucz model
;              is to be used with each standard star [nstd]. (This is the
;              output of "stype_standard".)
; 
; OPTIONAL INPUT:
;   norm   --  normalize the vectors before returning them (between 
;              5700 - 6300 A -- avoiding the last 200 pixels before the
;              dichroic)
;
; OUTPUT:
;   Vectors which represent the smoothed ratio of (data/model) for each
;   standard star.  [npix, nstd]
;
; OPTIONAL OUTPUT:
;   corvivar - inverse variance corresponding to each flux correction vector.  
;              [npix, nstd]
;   cormed   - median of each flux correction vector between 5700 and 6300 A
;              (but avoiding the 200 pixels nearest the dichroic) [nstd]
;
; COMMENTS:
;   For each standard star the best fit Kurucz model (as determined by 
;   "stype_standard") is restored, redshifted, and linearly interpolated to 
;   match the wavelength grid of the data.  Each model is then reddened
;   using the SDF reddening at the RA/DEC of the standard star and the 
;   extinction curve used by SFD (O'donnell).  
;
; BUGS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   divideflat
;   djs_maskinterp
;   djs_median()
;   ext_odonnell
;   filter_thru()
;   kurucz_restore
;   linterp
;   skymask() 
; 
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   12-Aug-2003  Create by C. Tremonti, Steward Observatory
;-
;------------------------------------------------------------------------------

function spdata2model_ratio, loglam, stdflux, stdivar, stdmask, stdinfo, $
         corvivar = corvivar, cormed = cormed, norm = norm

   ;--------------
   ; Read in Kurucz model files

   kurucz_restore, kwave, kflux, kindx = kindx

   ;-------------------
   ; Mask out bad pixels and regions dominated by sky-sub residuals

   ;stdivar = skymask(stdivar, 0, stdmask, ngrow=3)  ; Does this help???
   stdflux = djs_maskinterp(stdflux, stdivar EQ 0, iaxis=0, /const)

   ;-----------------
   ; Compute the flux correction vector from the ratio of model/data

   cspeed = 2.99792458e5
   npix = n_elements(stdflux[*,0])
   nstd = n_elements(stdflux[0,*])
   corvector = fltarr(npix, nstd)
   corvivar = fltarr(npix, nstd)
   cormed = fltarr(nstd)
   wave = 10.0^loglam 

   for istd=0, nstd-1 do begin
     model_index = (where(kindx.model eq stdinfo[istd].model))[0]
     kwave_full = kwave*(1 + stdinfo[istd].v_off/cspeed)
     kflux_full = kflux[*,model_index]
     linterp, kwave_full, kflux_full, wave, kfluxi

     ;------------
     ; Get extinction from SFD maps
   
     A_v = 3.1 * stdinfo[istd].e_bv_sfd
     a_odonnell = ext_odonnell(wave, 3.1)
     red_kflux = kfluxi * exp(-1 * a_odonnell * A_v / 1.086) 

     ;-------------
     ; Get zeropoint from phot fiber mag 

     ; choose guiding center as either r band (2) or g (1)
     scalefactor = 10.0^(-0.4 * (stdinfo[istd].mag[2] - $
                                 stdinfo[istd].red_model_mag[2]))
     red_kflux = red_kflux * scalefactor / 1e-17

     ;-----------
     ; Divide star flux in counts by model flux in erg/s/cm^2/A

     fluxvect = stdflux[*,istd]
     fluxvivar = stdivar[*,istd]
     divideflat, fluxvect, invvar=fluxvivar, red_kflux, minval = 0.1
   
     ;-----------
     ; Smooth the flux vector to reduce noise (do we need to smooth invar?)

     fluxvect = djs_median(fluxvect, width = 75, boundary = 'reflect')
     corvector[*,istd] = smooth(fluxvect, 25, /NAN)
     corvivar[*, istd] = fluxvivar

     ;-----------
     ; Normalize in the dichroic region but avoiding the exact edges

     norm_indx = where(wave gt 5700 and wave lt 6300 and $
                       wave lt max(wave) - 200 and wave gt min(wave) + 200)

     cormed[istd] = djs_median(corvector[norm_indx,istd])
    
     if keyword_set(norm) then begin
       corvector[*,istd] = corvector[*,istd] / cormed[istd]
       corvivar[*, istd] = corvivar[*,istd] * cormed[istd]^2
     endif
   endfor

   return, corvector
end
;------------------------------------------------------------------------------
