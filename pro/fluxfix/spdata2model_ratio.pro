;+
;-
;------------------------------------------------------------------------------
function spdata2model_ratio, loglam, stdflux, stdivar, stdmask, stdinfo, $
         corvivar = corvivar

   ;--------------
   ; Read in Kurucz model files

   kurucz_restore, kwave, kflux, kindx = kindx

   ;-------------------
   ; Mask out bad pixels and regions dominated by sky-sub residuals

   stdivar = skymask(stdivar, stdmask, ngrow=5)
   stdflux = djs_maskinterp(stdflux, stdivar EQ 0, iaxis=0, /const)

   ;-----------------
   ; Compute the flux correction vector from the ratio of model/data

   cspeed = 2.99792458e5
   npix = n_elements(stdflux[*,0])
   nstd = n_elements(stdflux[0,*])
   corvector = fltarr(npix, nstd)
   corvivar = fltarr(npix, nstd)
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
     scalefactor = 10.0^(-0.4 * (stdinfo[istd].mag[1] - $
                                 stdinfo[istd].red_model_mag[1]))
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

   endfor

   return, corvector
end
;------------------------------------------------------------------------------
