;+
; NAME:
;   stype_standard
;
; PURPOSE:
;   Spectral type the standard stars and save the result as a FITS file.
;   The input spectra are normalized and then compared to a grid of 
;   Kurucz models produced by the SPECTRUM code (R.0. Gray & C. J. Corbally,
;   1994, AJ, 107, 742) and convolved to SDSS resolution. Since nearly all 
;   of the interesting spectral features are in the blue, data from the 
;   blue camera only can be used.
;
; CALLING SEQUENCE:
;   stype_standard, loglam, flux, invvar, plugmap, outname 
;
; INPUTS:
;   loglam     - Wavelength array of input spectra in log10(Angstroms) [npix]
;   flux       - Array of standard star spectra [npix, nstar]
;   invvar     - Inverse variance of standard star spectra [npix,nstar]
;   plugmap    - Plugmap corresponding to the input standard stars [nstar]
;   outname    - Name of output file. It must be of the from
;                ?????-pppp-mmmm-s.fits where, pppp is the plateid, 
;                mmmmm is the mjd and s is the spectrograph ID.
;                e.g. spStd-0519-52283-1.fits 
;
; OUTPUT:  A FITS binary table containing the information on the best fit 
;          models.  Diagnostic plots are also produced.  Each spectophoto 
;          standard is shown normalized with the best fit spectrum plotted 
;          over top in red.  
;
; COMMENTS:  A file containing the Kurucz model data is required.  It is 
;            called "kurucz_stds_interp.fit" and it should reside in 
;            IDLSPEC2D_DIR/etc.) The 0th HDU contains the flux in
;            ergs/s/cm^2/A -- the absolute value of the flux arbitrary.
;            The spectra have been convolved to SDSS resolution (approximately)
;            and rebinned to dloglam = 1e-4.  The 1st HDU contains the 
;            normalized flux.  The 2nd HDU contains information about each 
;            model such as effective temperature, surface gravity, and 
;            metallicity.  The wavelength information is in the header.
;
; BUGS:
;            The SFD dust maps are needed to estimate the foreground extinction 
;            to add to the models.  To acomodate this a new environment
;            variable "DUST_DIR" is required.  In the future the dust 
;            maps will be part of IDL_SPEC2D.
;
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   correct_dlam
;   divideflat
;   djs_filepath()
;   djs_maskinterp()
;   djs_median()
;   djs_oplot
;   djs_plot
;   fibermask_bits()
;   fileandpath()
;   mrdfits()
;   mwrfits
;   sxpar()
;   traceset2xy
;
;
; INTERNAL SUPPORT ROUTINES
;   qgoodfiber()
;   kurucz_match()
;
; REVISION HISTORY:
;   28-Sep-2002  Written by C. Tremonti
;-
;------------------------------------------------------------------------------
function qgoodfiber, fibermask
   qgood = ((fibermask AND fibermask_bits('NOPLUG')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADTRACE')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADFLAT')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADARC')) EQ 0) $
       AND ((fibermask AND fibermask_bits('MANYBADCOLUMNS')) EQ 0) $
       AND ((fibermask AND fibermask_bits('NEARWHOPPER')) EQ 0) $
       AND ((fibermask AND fibermask_bits('MANYREJECTED')) EQ 0)
   return, qgood
end

;------------------------------------------------------------------------------
; Find the Kurucz model which best matches a star

pro kurucz_match, wave, nflux, nivar, nkflux, kindx, fiber, $
    plottitle = plottitle

   
   ;---------------------
   ; Compute the chisq  

   nummod = n_elements(nkflux[0,*]) 
   chi2 = fltarr(nummod)
   nivar[where(wave gt 5570 and wave lt 5590)] = 0 ; Avoid 5577 sky line

   for imod = 0, nummod - 1 do begin
  ; Use only blue side
     goodpix = where(wave gt 3850 and wave lt 5800 and nivar ne 0) 
 
     chi2[imod] = total((nflux[goodpix] -  nkflux[goodpix,imod])^2 $
                         * nivar[goodpix], /NAN) 
   endfor 
  
;stop 
   ;-----------------
   ; Store model info in structure that was passed in 

   minchi = min(chi2, chiindx)
   fiber.model = kindx[chiindx].model
   fiber.chi2 = minchi
   fiber.feh = kindx[chiindx].feh
   fiber.teff = kindx[chiindx].teff
   fiber.g = kindx[chiindx].g
   fiber.model_mag = kindx[chiindx].mag
 
   ;--------------
   ; Plot normalized star spectra and best fit model

   if keyword_set(plottitle) then begin
     origpmulti = !P.MULTI
     !P.MULTI = [0, 1, 2]
     gpix = where(wave gt 3810 and wave lt 9190)
     djs_plot, wave[gpix], nflux[gpix], psym=10, xr=[3800,5000], $
               yr=[0,1.5], /xs, title = plottitle, $
               xtitle = 'Wavelength (A)', ytitle = 'Normalized Flux'
     djs_oplot, wave[gpix], nkflux[gpix,chiindx], color='red', psym=10

     xyouts, 3850, 1.35, 'Fiber: ' + string(fiber.fiberid, format='(I3)')
     xyouts, 4050, 1.35, 'S/N = ' + string(fiber.sn, format='(F4.1)')
     xyouts, 4300, 1.35, 'Model: ' + fiber.model
     xyouts, 4750, 1.35, 'Chi2 = ' + string(fiber.chi2, format='(F6.2)')

     djs_plot, wave[gpix], nflux[gpix], psym=10, xr=[5000,6800], yr=[0,1.5], $
               /xs, xtitle = 'Wavelength (A)', ytitle = 'Normalized Flux'
     djs_oplot, wave[gpix], nkflux[gpix,chiindx], color='red', psym=10

     !P.MULTI = origpmulti
   endif
end

;------------------------------------------------------------------------------
pro stype_standard, loglam, flux, invvar, plugmap, outname, $
    nonsdss = nonsdss, smoothpix = smoothpix 

   dloglam = 1.0d-4 
   cspeed = 2.99792458e5
   wave = 10.0^loglam 

   ;---------
   ; Extract Plate and MJD and spectrograph ID from output name
   words = strsplit(outname, '-', /extract)
   plate = words[1]
   mjd = words[2]
   side = words[3]

   nsphoto = n_elements(plugmap)

   ;------------- 
   ; Normalize the spectrum

   nflux = rectify(flux, invvar, nivar = ninvvar, /mask, wave = wave)
   nflux = djs_maskinterp(nflux, ninvvar EQ 0, iaxis=0, /const)

   ;-------------
   ; Set up structure to hold fstar data

   fiber = {fiber1, plateid:0, mjd: 0L, fiberid: 0, side: 0, $
            ra: 0.0, dec: 0.0, mag: [0.0, 0.0, 0.0, 0.0, 0.0], $
            xfocal: 0.0, yfocal: 0.0, sn: 0.0, v_off: 0.0, $
            model: ' ', chi2: 0.0, teff: 0.0, g: 0.0, feh: 0.0, $
            e_bv_sfd: 0.0, model_mag: [0.0, 0.0, 0.0, 0.0, 0.0], $
            red_model_mag: [0.0, 0.0, 0.0, 0.0, 0.0]}

   stdstar = make_array(value = fiber, dim = nsphoto)
   
   if not keyword_set(nonsdss) then begin 
     stdstar[*].plateid = plate
     stdstar[*].mjd = mjd
     stdstar[*].side = side
     stdstar.fiberid = plugmap.fiberid
     stdstar.ra = plugmap.ra
     stdstar.dec = plugmap.dec
     stdstar.mag = plugmap.mag
     stdstar.xfocal = plugmap.xfocal
     stdstar.yfocal = plugmap.yfocal
     stdstar.sn = djs_median(nflux * sqrt(ninvvar), 1)

     ;------------------
     ; Get extinction from SFD maps

     dustdir = getenv('DUST_DIR') 
     if keyword_set(dustdir) then begin
       dustdir = dustdir + '/maps/'
       glactc, stdstar.ra, stdstar.dec, 2000, l, b, 1, /degre
       stdstar.e_bv_sfd = dust_getval(l, b, ipath = dustdir, /interp)
     endif else begin
       splog, 'Environment variable DUST_DIR not set -- no path to SFD maps'
       splog, 'WARNING: No dust correction will be applied to standard stars'
     endelse
   endif
 
   ;--------------
   ; Read in Kurucz model files

   kurucz_restore, kwave, kflux, nkflux = nkflux, kindx = kindx, $
                   smoothpix = smoothpix
   nk = n_elements(kindx)

   ; Compute the offset of the models from the data in pixels 
   model_offset = (loglam[0] - alog10(kwave[0])) / dloglam

   ; Get rid of bad pixels so they doesn't screw up the velocity determination 
   bad = where(finite(nflux) ne 1 or finite(ninvvar) ne 1)
   if bad[0] ne -1 then nflux[bad] = 1 
   if bad[0] ne -1 then ninvvar[bad] = 0 

   ;------------------
   ; Use fiducial model with T_eff = 6000, Log(gravity) = 4, and 
   ; log(Fe/H) = -1 to determine the velocity offset of the stars in pixels

   fidmodel = where(kindx.teff eq 6000 and kindx.g eq 4 and $
                    kindx.feh eq -1.0)
   nkfluxi = nkflux[*, fidmodel] 
  
   bad = where(kwave lt min(wave) + 200 or kwave gt max(wave) - 200) 
   starmask = intarr(n_elements(kwave)) + 1
   starmask[bad] = 0  
  
   zans = zcompute(nflux, ninvvar, nkfluxi, starmask, $
                   poffset = model_offset, pmin = -20.0, pmax = 20.0)   

  ; convert from pix to km/s
  stdstar.v_off = (10.^(dloglam * zans[*].z) -1 )*cspeed

  ;----------------
  ; Find the best matched Kurucz model for each star & store result in
  ; structure  

  for iobj=0, nsphoto-1 do begin
    fiber = stdstar[iobj] 
   
    rkwave = kwave * (1 + stdstar[iobj].v_off/cspeed) 
    rkflux = fltarr(n_elements(wave), nk)
    for kobj = 0, nk - 1 do $
       rkflux[*,kobj] = interpol(nkflux[*,kobj], rkwave, wave)

    kurucz_match, wave, nflux[*,iobj], ninvvar[*,iobj], $
                  rkflux, kindx, fiber, plottitle = outname
     
    stdstar[iobj] = fiber
  endfor 

  ;--------------
  ; Save the output as a FITS binary table 
  mwrfits, stdstar, outname, /create

  return
end
;------------------------------------------------------------------------------
