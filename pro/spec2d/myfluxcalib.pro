;------------------------------------------------------------------------------
; Re-sort the 2nd plug-map to match that of the first file, so
; that PLUGMAP2[INDX] = PLUGMAP1
function resortplugmap, plugmap1, plugmap2

   nfiber = n_elements(plugmap1)
   indx = lonarr(nfiber)

   for ifiber=0, nfiber-1 do begin
      adist = djs_diff_angle(plugmap1.ra, plugmap1.dec, $
       plugmap2[ifiber].ra, plugmap2[ifiber].dec, units='degrees')
      indx[ifiber] = $
       (where(adist LT 2./3600. AND strtrim(plugmap1.objtype,2) NE 'NA'))[0]
   endfor

   return, indx
end

;------------------------------------------------------------------------------
function fluxfit, loglam, objflux, objivar, color=color, mags=mags

   ;----------
   ; Change inputs such that wavelengths are in ascending order.
   ; Need to do this to avoid bug currently in the B-spline code ???

   if (loglam[1] LT loglam[0]) then begin
      loglam = reverse(loglam)
      objflux = reverse(objflux)
      objivar = reverse(objivar)
   endif

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
;stop ; ???
;splot,10^loglam,fitflux
;soplot,10^loglam,bspline_valu(loglam,sset),color='red'

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
;   filename   - Name(s) of input file names, typically one blue and
;                one red file.  The pluggings need not be the same,
;                but there must be at least one good SPECTROPHOTO or
;                REDDEN standard in common.
;   calibfile  - Name(s) of output calibration files, one per FILENAME

pro myfluxcalib, filename, calibfile, colors=colors, adderr=adderr, $
    rset=rset, bset=bset

   dloglam = 1.0d-4 ; ???

   nfile = n_elements(filename)
   if (n_elements(calibfile) NE nfile) then $
    message, 'Dimensions of FILENAME and CALIBFILE do not agree'

   ;----------
   ; Assign scores to each object based upon its color relative
   ; to the color BD+17 4708.

   bd17mag = [10.56, 9.64, 9.35, 9.25, 9.23]
   bd17color = bd17mag[0:3] - bd17mag[1:4]

   ;----------
   ; Read the plug maps and masks and decide which star to use

   for ifile=0, nfile-1 do begin
      thismask = mrdfits(filename[ifile],2)
      thisplug = mrdfits(filename[ifile],5)
      goodfiber = qgoodfiber(thismask[0,*])

      if (ifile EQ 0) then begin
         ; Save the first plug-map for comparison to the others
         plugmap = thisplug
         plugindx = lindgen(n_elements(thisplug))

         ; Save the first good-fiber mask
         goodmask = goodfiber

         ; Assign scores based upon colors and magnitudes
         colordiff = thisplug.mag[0:3] - thisplug.mag[1:4] 
         for i=0, 3 do $
          colordiff[i,*] = colordiff[i,*] - bd17color[i]
         colordiff = sqrt(total(colordiff^2,1))
         colorscore = thisplug.mag[2] + 40.0 * colordiff ; Lower score is better
      endif else begin
         ; Re-sort the plug-map to match that of the first file
         indx = resortplugmap(plugmap, thisplug)
         plugindx = [ [plugindx], [indx] ]

         ; Multiply this good-fiber mask with the others,
         ; setting to zero if the object does not exist in this plugging.
         goodmask = goodmask * goodfiber[indx>0] * (indx NE -1)
      endelse
   endfor

   ;----------
   ; Select a single flux-calibration star
   ; It's index is IBEST in the first file, or
   ;   where(plugindx[*,ifile] EQ ibest) in the other files.

   indx = where(strtrim(plugmap.objtype) EQ 'SPECTROPHOTO_STD' AND goodmask)
   if (indx[0] EQ -1) then begin
      indx = where(strstrim(plugmap.objtype) EQ 'REDDEN_STD' AND goodmask)
      if (indx[0] EQ -1) then $
       message, 'No SPECTROPHOTO or REDDEN stars for flux calibration' $
      else $
       splog, 'WARNING: Must use REDDEN std instead of SPECTROPHOTO for fluxing'
   endif

   bestscore = min(colorscore[indx], ibest)
   ibest = indx[ibest]

   splog, 'Spectrophoto star is fiberid = ', plugmap[ibest].fiberid, $
    ' in file ', filename[0]
   splog, 'Spectrophoto mag = ', plugmap[ibest].mag
   splog, 'Spectrophoto score = ', bestscore

   ;----------

   for ifile=0, nfile-1 do begin
      ;----------
      ; Read in the flux, errors, and wavelengths

      objflux = mrdfits(filename[ifile],0)
      objivar = mrdfits(filename[ifile],1)
      wset = mrdfits(filename[ifile],3)
      traceset2xy, wset, 0, loglam

      ;----------
      ; Add an additional error term equal to ADDERR of the flux.

      if (keyword_set(adderr)) then begin
         gmask = objivar NE 0 ; =1 for good points
         objivar = 1.0 / ( 1.0/(objivar + (1-gmask)) $
          + (adderr * (objflux>0))^2 ) * gmask
      endif

      ;----------
      ; Make a map of the size of each pixel in delta-(log10-Angstroms),
      ; and re-normalize the flux to ADU/(dloglam)

      correct_dlam, objflux, objivar, wset, dlam=dloglam

      ;----------
      ; Do the actual fit

      jbest = (where(plugindx[*,ifile] EQ ibest))[0]
      calibset = fluxfit(loglam[*,jbest], objflux[*,jbest], $
       objivar[*,jbest], color=colors[ifile], mags=plugmap[ibest].mag)

      if colors[ifile] EQ 'r' then rset = calibset
      if colors[ifile] EQ 'b' then bset = calibset

;stop ; ???
;set_plot,'x'
;junk=bspline_valu(loglam[*,jbest], calibset)
;splot,loglam[*,jbest],objflux[*,jbest]
;splot, loglam[*,jbest], objflux[*,jbest]/junk

      ;----------
      ; Create header cards describing the fit range

      hdr = ['']
      indx = where(objivar[*,jbest] NE 0)
      wavemin = 10.^min(loglam[indx,jbest])
      wavemax = 10.^max(loglam[indx,jbest])
      sxaddpar, hdr, 'WAVEMIN', wavemin
      sxaddpar, hdr, 'WAVEMAX', wavemax

      mwrfits, 0, calibfile[ifile], hdr, /create
      mwrfits, calibset, calibfile[ifile]
   endfor

   return
end
;------------------------------------------------------------------------------
