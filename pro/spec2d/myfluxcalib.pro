
;------------------------------------------------------------------------------
function fluxfit, loglam, objflux, objivar, color=color, mags=mags

   wave = 10^loglam
   logmin = min(loglam)
   logmax = max(loglam)
   fitivar = objivar

   ;----------
   ; Read the spectrum of an F8 star

   f8file = filepath('f8vspline.dat', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, f8file, f8wave, f8flux

   ;----------
   ; Divide the input spectrum by that of the F8 star

   intrinspl = spl_init(alog10(f8wave), f8flux)
   f8spline = spl_interp(alog10(f8wave), f8flux, intrinspl, loglam)
   fitflux = objflux / f8spline

   ;----------
   ; Mask out around absorption lines

   absfile = filepath('f8v.abs', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, absfile, absmin, absmax

   for i=0, n_elements(absmin)-1 do begin
      ipix = where(wave GT absmin AND wave LT absmax)
      if (ipix[0] NE -1) then fitivar[ipix] = 0
   endfor

   ;----------
   ; Select break points for spline

   if (color EQ 'b') then bkptfile = filepath('blue.bkpts', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   if (color EQ 'r') then bkptfile = filepath('red.bkpts', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   readcol, bkptfile, allbkpts
   allbkpts = alog10(allbkpts) ; Convert to log-10 Angstroms

   ibk = where(allbkpts GE logmin AND allbkpts LE logmax, nbk)
   if (nbk LT 4) then $
    message, 'Error selecting break points'
   allbkpts = [logmin, allbkpts[ibk[1:nbk-2]], logmax]

   ;----------
   ; Do the spline fit

   sset = bspline_iterfit(loglam, fitflux, $
    invvar=fitivar, nord=3, bkpt=allbkpts, upper=3, lower=3, $
    maxrej=0.05*n_elements(objflux))

   ;----------
   ; Scale the flux
;       at v=0, flux at 5556A is 956 photons/cm^2/A/s
;               = 3.52e-9 ergs/cm^2/s/A
;       scale to r', with 10^(-r'/2.5)
;       and return in units to 1e-17 ergs/cm/s/A
;       so factor in exponent is 10^((21.37 - r')/2.5)
;
;       AB magnitude, the scaling this assumes
;       the AB magnitude of BD+17 is 9.4 at 5560
;
;       we then use f_nu = 10^-0.4*(AB+48.6)
;       and then f_lambda = f_nu * c / (lambda)^2
; 
;       c is 3.0e18 Ang/s
;       lambda is in Ang

   scalefac = 10.^((21.37 - mags[2]) / 2.5)
   sset.coeff = sset.coeff / scalefac

; ???
;junk = bspline_valu(loglam, sset)
;splot, 10^loglam, fitflux,xr=[5500,6500]  
;soplot,10^loglam,junk*scalefac,color='red'
;stop
   return, sset
end

;------------------------------------------------------------------------------
function qgoodfiber, fibermask
   qgood = ((fibermask AND fibermask_bits('NOPLUG')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADTRACE')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADFLAT')) EQ 0) $
       AND ((fibermask AND fibermask_bits('BADARC')) EQ 0) $
       AND ((fibermask AND fibermask_bits('MANYBADCOLUMNS')) EQ 0) $
       AND ((fibermask AND fibermask_bits('MANYREJECTED')) EQ 0)
   return, qgood
end

;------------------------------------------------------------------------------
pro myfluxcalib, bsciname, rsciname, bcalibset, rcalibset

   dloglam = 1.0d-4 ; ???

   ;----------
   ; Read in the flux, errors, mask, wavelengths and plug-map

;   bflux = mrdfits(bsciname,0)
bflux = mrdfits(bsciname,0) * mrdfits(bsciname,6) ; ???
;   bivar = mrdfits(bsciname,1)
bivar = mrdfits(bsciname,1) * mrdfits(bsciname,6) ; ???
   bmask = mrdfits(bsciname,2)
   bwset = mrdfits(bsciname,3)
   bplug = mrdfits(bsciname,5)
   traceset2xy, bwset, junk, bloglam

   rflux = mrdfits(rsciname,0)
rflux = mrdfits(rsciname,0) * mrdfits(rsciname,6) ; ???
;   rivar = mrdfits(rsciname,1)
rivar = mrdfits(rsciname,1) * mrdfits(rsciname,6) ; ???
   rmask = mrdfits(rsciname,2)
   rwset = mrdfits(rsciname,3)
   rplug = mrdfits(rsciname,5)
   traceset2xy, rwset, junk, rloglam

   ;----------
   ; Make a map of the size of each pixel in delta-(log10-Angstroms),
   ; and re-normalize the flux to ADU/(dloglam)

   npix = (size(bflux,/dimens))[0]
   dlogimg = [ bloglam[1,*] - bloglam[0,*], $
    0.5 * (bloglam[2:npix-1,*] - bloglam[0:npix-3,*]), $
    bloglam[npix-1,*] - bloglam[npix-2,*] ]
   dlogimg = abs(dlogimg)
   divideflat, bflux, bivar, (dlogimg/dloglam), minval=0

   npix = (size(bflux,/dimens))[0]
   dlogimg = [ rloglam[1,*] - rloglam[0,*], $
    0.5 * (rloglam[2:npix-1,*] - rloglam[0:npix-3,*]), $
    rloglam[npix-1,*] - rloglam[npix-2,*] ]
   dlogimg = abs(dlogimg)
   divideflat, rflux, rivar, (dlogimg/dloglam), minval=0

   ;----------
   ; Decide if a fiber is good

   bgoodfiber = qgoodfiber(bmask[0,*])
   rgoodfiber = qgoodfiber(rmask[0,*])

   ;----------
   ; Assign scores to each object based upon its color relative
   ; to the color BD+17 4708.

   bd17mag = [10.56, 9.64, 9.35, 9.25, 9.23]
   bd17color = bd17mag[0:3] - bd17mag[1:4]

   colordiff = bplug.mag[0:3] - bplug.mag[1:4] 
   for i=0, 3 do $
    colordiff[i,*] = colordiff[i,*] - bd17color[i]
   colordiff = sqrt(total(colordiff^2,1))

   colorscore = bplug.mag[2] + 40.0 * colordiff ; Lower score is better

   ;----------
   ; Select a single flux-calibration star

   indx = where(strtrim(bplug.objtype) EQ 'SPECTROPHOTO_STD' $
    AND bgoodfiber AND rgoodfiber)
   if (indx[0] EQ -1) then begin
      indx = where(strstrim(bplug.objtype) EQ 'REDDEN_STD' $
       AND bgoodfiber AND rgoodfiber)
      if (indx[0] EQ -1) then $
       message, 'No SPECTROPHOTO or REDDEN stars for flux calibration' $
      else $
       splog, 'WARNING: Must use REDDEN std instead of SPECTROPHOTO for fluxing'
   endif

   bestscore = min(colorscore[indx], ibest)
   ibest = indx[ibest]

   splog, 'Spectrophoto star is fiberid = ', bplug[ibest].fiberid
   splog, 'Spectrophoto mag = ', bplug[ibest].mag
   splog, 'Spectrophoto score = ', bestscore

   ; Reverse sort the blue spectrum so that it's in ascending order ???
   bcalibset = fluxfit(reverse(bloglam[*,ibest]), reverse(bflux[*,ibest]), $
    reverse(bivar[*,ibest]), color='b', mags=bplug[ibest].mag)
   rcalibset = fluxfit(rloglam[*,ibest], rflux[*,ibest], $
    rivar[*,ibest], color='r', mags=rplug[ibest].mag)

   return
end
;------------------------------------------------------------------------------
