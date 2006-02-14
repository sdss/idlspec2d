;+
; NAME:
;   nebularsky
;
; PURPOSE:
;   Fit nebular emission lines in SDSS sky + galaxy spectra
;
; CALLING SEQUENCE:
;   nebularsky, [ plate, mjd=, lambda=, fitrange=, zlimits=, siglimits=, $
;    npoly=, fitflux=, lwidth=, outfile=, /debug ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   plate      - List of plates; default to using the PLATELIST procedure
;                to select all reduced plates without quality set to 'bad'
;   mjd        - MJD corresponding to each PLATE
;   lambda     - Wavelengths for sky emission lines [vacuum Ang]
;   fitrange   - Fitting region; default to [6650,6800] vacuum Ang
;   zlimits    - Redshift limits for all emission lines; default to
;                [-300,300]/3e5
;   siglimits  - Velocity dispersion limits for all emission lines; default to
;                [50,200] km/s
;   npoly      - Number of polynomial terms for sky spectrum; default to 3
;                (quadratic)
;   fitflux    - Subtract out either the 'synflux' or 'lineflux' spectrum;
;                default to 'lineflux', which is the velocity-dispersion
;                broadened Elodie templates plus emission lines.
;   lwidth     - Full width for masking around possible galaxy emission lines;
;                default to 0.002 in log-wavelenghth (about 1382 km/s)
;   outfile    - Output file; default to 'nebular.fits'
;   debug      - If set, then make debugging plots, and wait for keystroke
;                after each plot
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine creates the output file, and then appends to it one
;   plate at a time.
;
;   All wavelengths are in vacuum, and velocities are barycentric.
;
; EXAMPLES:
;   nebularsky, 231, mjd=51456
;
; BUGS:
;   Discard galaxies that are too bright ???
;   Make use of the instrumental response at each line center ???
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/emlines.par
;
; PROCEDURES CALLED:
;   copy_struct_inx
;   linebackfit()
;   mwrfits_chunks
;   pixelmask_bits()
;   readonespec
;   readspec
;   soplot
;   splog
;   splot
;   struct_selecttags()
;
; REVISION HISTORY:
;   12-Jan-2006  Written by A. West & D. Schlegel, Berkeley
;-
;------------------------------------------------------------------------------
pro nebularsky, plate, mjd=mjd1, lambda=lambda1, fitrange=fitrange1, $
 zlimits=zlimits1, siglimits=siglimits1, fitflux=fitflux1, lwidth=lwidth1, $
 npoly=npoly1, outfile=outfile1, debug=debug

   if (keyword_set(lambda1)) then lambda = lambda1 $
    else lambda = [6718.294, 6732.678]
   if (keyword_set(fitrange1)) then fitrange = fitrange1 $
    else fitrange = [6650,6800]
   if (keyword_set(zlimits1)) then zlimits = zlimits1 $
    else zlimits = [-300.,300.]/3e5
   if (keyword_set(siglimits1)) then siglimits = siglimits1 $
    else siglimits = [50.,200.]
   if (n_elements(npoly1) GT 0) then npoly = npoly1 $
    else npoly = 3
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile = 'nebular.fits'
   if (keyword_set(fitflux1)) then fitflux = strlowcase(fitflux1) $
    else fitflux = 'lineflux'
   if (keyword_set(lwidth1)) then lwidth = lwidth1 $
    else lwidth = 0.002
   if (fitflux NE 'synflux' AND fitflux NE 'lineflux') then $
    message, 'Invalid string for FITFLUX'

   res_all = 0
   nline = n_elements(lambda)

   zindex = fltarr(n_elements(lambda))
   windex = fltarr(n_elements(lambda))
   zguess = fltarr(n_elements(lambda))

   csize = 1.6
   select_tags = ['PLATE','MJD','FIBERID','PLUG_RA','PLUG_DEC']

   ;----------
   ; Read line lists and convert to vacuum

   linefile = filepath('emlines.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   yanny_read, linefile, pdata
   linelist = *pdata[0]
   yanny_free, pdata

   vaclambda = linelist.lambda
   airtovac, vaclambda
   linelist.lambda = vaclambda

   ;----------
   ; Select the list of  plates

   if (keyword_set(plate)) then begin
      nplate = n_elements(plate)
      plist = replicate(create_struct('PLATE',0L,'MJD',0L), nplate)
      plist.plate = plate
      if (keyword_set(mjd1)) then begin
         if (n_elements(mjd1) NE nplate) then $
          message, 'Number of elements in PLATE and MJD do not agree'
         plist.mjd = mjd1
      endif
   endif else begin
      platelist, plist=plist
      plist = plist[where(strmatch(plist.status1d,'Done*') $
       AND strmatch(plist.platequality,'bad*') EQ 0, nplate)]
   endelse

   splog, 'Number of plates = ', nplate

   for iplate=0L, nplate-1 do begin
      t0 = systime(1)
      readspec, plist[iplate].plate, mjd=plist[iplate].mjd, $
       zans=zans, plug=plug, /silent
      if (keyword_set(zans[0])) then begin
         plist[iplate].mjd = zans[0].mjd ; Fill in if this was zero
         splog, 'Working on PLATE= ', plist[iplate].plate, $
          ' MJD= ', plist[iplate].mjd
         zans_trim = struct_selecttags(zans, select_tags=select_tags)
         qsky = (zans.zwarning AND 1) NE 0 AND (zans.zwarning AND 2^1+2^7) EQ 0
         qgalaxy = zans.zwarning EQ 0 AND strmatch(zans.class,'GALAXY*')
         qstar = zans.zwarning EQ 0 AND strmatch(zans.class,'STAR*')

         for ifiber=0, 639 do begin
            if (qsky[ifiber] OR qgalaxy[ifiber] OR qstar[ifiber]) then begin

               ; Read in all the individual exposures for this object
               ; (red cameras only)
               readonespec, plist[iplate].plate, mjd=plist[iplate].mjd, $
                zans[ifiber].fiberid, cameras='r', $
                flux=flux, invvar=invvar, mask=mask, loglam=loglam, $
                sky=sky, synflux=synflux, lineflux=lineflux, /silent
               invvar = invvar * ((mask AND pixelmask_bits('COMBINEREJ')) EQ 0)
               if (fitflux EQ 'lineflux') then synflux = lineflux

               if (keyword_set(npoly)) then begin
                  ; Generate background terms of orderr NPOLY, with a set of
                  ; background terms for each exposure
                  ndim = size(flux,/n_dimen)
                  dims = size(flux,/dimens)
                  npix = dims[0]
                  if (ndim EQ 1) then nobs = 1 else nobs = dims[1]
                  background = fltarr(npix,nobs,nobs,npoly)
                  for iobs=0, nobs-1 do $
                   background[*,iobs,iobs,*] = poly_array(npix,npoly)
                  background = reform(background,npix*nobs,nobs*npoly)
               endif else begin
                  npoly = 0
               endelse

               ; Set SYNFLUX=0 if this is a sky fiber
               synflux = synflux * (1-qsky[ifiber])

               ; Discard wavelengths near emission lines in the galaxy frame
               if (qgalaxy[ifiber]) then begin
                  for iline=0, n_elements(linelist)-1 do begin
                     thislam = alog10(linelist[iline].lambda*(1+zans[ifiber].z))
                     invvar = invvar * ((loglam LT thislam - 0.5*lwidth) $
                      OR (loglam GT thislam + 0.5*lwidth))
                  endfor
               endif

               ; Fit only within the fitting range
               invvar = invvar * (loglam GE alog10(fitrange[0]) $
                AND loglam LE alog10(fitrange[1]))

;               res1 = linebackfit(lambda, reform(loglam,npix*nobs), $
;                reform(flux+sky-synflux,npix*nobs), $
;                invvar=reform(invvar,npix*nobs), background=background, $
;                zindex=zindex, windex=windex, zguess=zguess, yfit=yfit)

               ii = where(invvar NE 0, ngpix)
               if (ngpix EQ 0) then ii = 0
               res1 = linebackfit(lambda, loglam[ii], $
                flux[ii]+sky[ii]-synflux[ii], $
                invvar=invvar[ii], background=background[ii,*], $
                zindex=zindex, windex=windex, zguess=zguess, $
                zlimits=zlimits, siglimits=siglimits, yfit=yfit1)
               yfit = fltarr(npix,nobs)
               yfit[ii] = yfit1

               if (NOT keyword_set(res_all)) then begin
                  res_blank = create_struct(zans_trim[0], res1[0])
                  struct_assign, {junk:0}, res_blank
                  res_all = replicate(res_blank, nline, 640)
               endif
               index_to = ifiber*nline+lindgen(nline)
               copy_struct_inx, replicate(zans_trim[ifiber],nline), res_all, $
                index_to=index_to
               copy_struct_inx, res1, res_all, index_to=index_to

               if (keyword_set(debug)) then begin
                  xplot = 10^loglam
                  yplot = flux+sky-synflux
                  if (ngpix GT 1) then begin
                     xrange = minmax(xplot[ii])
                     yrange = minmax(yplot[ii])
                  endif else begin
                     xrange = minmax(xplot)
                     yrange = minmax(yplot)
                  endelse
                  title = string(plist[iplate].plate, plist[iplate].mjd, $
                   zans[ifiber].fiberid, $
                   format='("Plate ", i4, " MJD ", i5, " Fiber ", i3)')
                  splot, [0], [0], /nodata, $
                   xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
                   charsize=csize, xtitle='Wavelength [Ang]', ytitle='Flux', $
                   title=title
                  for iobs=0, nobs-1 do begin
                     jj = where(invvar[*,iobs] NE 0, ct)
                     if (ct GT 1) then begin
                        soplot, xplot[jj,iobs], yplot[jj,iobs]
                        soplot, xplot[jj,iobs], yfit[jj,iobs], color='red'
                     endif
                  endfor
                  thisclass = qsky[ifiber] ? 'SKY' : zans[ifiber].class
                  xyouts, total([0.9,0.1]*!x.crange), $
                   total([0.1,0.9]*!y.crange), thisclass, charsize=csize
                  cc = get_kbrd(1) 
               endif
            endif
         endfor
         mwrfits_chunks, res_all, outfile, create=(iplate EQ 0), $
          append=(iplate NE 0), /silent
         splog, 'Elapsed time for plate ', iplate+1, ' of ', nplate, $
          ' = ', systime(1)-t0, ' sec'
      endif else begin
         splog, 'Skipping PLATE= ', plist[iplate].plate, $
          ' MJD= ', plist[iplate].mjd, ' (not found)'
      endelse
   endfor

   return
end
;------------------------------------------------------------------------------
