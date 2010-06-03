;+
; NAME:
;   nebularsky
;
; PURPOSE:
;   Fit nebular emission lines in SDSS sky + galaxy spectra
;
; CALLING SEQUENCE:
;   nebularsky, [ plate, mjd=, lambda=, skyfile=, fitrange=, $
;    zlimits=, siglimits=, npoly=, fitflux=, lwidth=, $
;    outfile=, logfile=logfile, /onlysky, /create, /debug ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   plate      - List of plates; default to using the PLATELIST procedure
;                to select all reduced plates without quality set to 'bad'
;   mjd        - MJD corresponding to each PLATE
;   lambda     - Wavelengths for sky emission lines [vacuum Ang]; default to
;                using the galaxy emission line list at > 3800 Ang, but where
;                the [O_I] lines are allowed to have a velocity and width
;                independent of the other lines
;   skyfile    - FITS file containing the PCA components for the sky,
;                with the spectra in HDU #0, and the log-wavelengths
;                in HDU #1; this file can be generated with NEBULAR_PCASKY.
;   fitrange   - Fitting region in vacuum Ang; default to using all wavelengths
;   zlimits    - Redshift limits for all emission lines; default to
;                [-200,200]/3e5
;   siglimits  - Velocity dispersion limits for all emission lines; default to
;                [60,150] km/s
;   npoly      - Number of polynomial terms for sky spectrum; default to 3
;                (quadratic)
;   fitflux    - Subtract out either the 'synflux' or 'lineflux' spectrum;
;                default to 'synflux'
;                broadened Elodie templates plus emission lines.
;   lwidth     - Full width for masking around possible galaxy emission lines;
;                default to 0.002 in log-wavelenghth (about 1382 km/s)
;   onlysky    - If set, then only fit to good sky fibers; otherwise, fits
;                are performed to good sky fibers, galaxies, and stars
;   outfile    - Output FITS file; default to 'nebular.fits'
;   logfile    - Output log file; default to a separate log file per plate,
;                named logfile-$PLATE-$MJD.log.  If specified, then a single
;                log file by that name is used for all plates.
;   create     - If set, then create a new output file; default to appending
;                to an existing file if it already exists
;   debug      - If set, then make debugging plots, and wait for keystroke
;                after each plot
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine creates the output file, and then appends to it one
;   plate at a time.  Several instances of this procedure can safely
;   write to the same file.
;
;   All wavelengths are in vacuum, and velocities are barycentric.
;
;   All lines are constrained to have the same velocity and width, except
;   for the [O_I] lines since those can be due to solar activity.
;
; EXAMPLES:
;   nebularsky, 231, mjd=51456
;
; BUGS:
;   Make use of the instrumental response at each line center ???
;   Test the FITFLUX='lineflux' option ???
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/emlines.par
;
; PROCEDURES CALLED:
;   copy_struct_inx
;   djs_lockfile()
;   djs_unlockfile
;   linebackfit()
;   mrdfits()
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
pro nebularsky, plate, mjd=mjd1, lambda=lambda1, skyfile=skyfile1, $
 fitrange=fitrange, $
 zlimits=zlimits1, siglimits=siglimits1, fitflux=fitflux1, lwidth=lwidth1, $
 npoly=npoly1, onlysky=onlysky, outfile=outfile1, logfile=logfile, $
 create=create1, debug=debug

   if (keyword_set(skyfile1)) then skyfile = skyfile1 $
    else skyfile = 'pcasky.fits'
   if (keyword_set(zlimits1)) then zlimits = zlimits1 $
    else zlimits = [-200.,200.]/3e5
   if (keyword_set(siglimits1)) then siglimits = siglimits1 $
    else siglimits = [60.,150.]
   if (n_elements(npoly1) GT 0) then npoly = npoly1 $
    else npoly = 3
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile = 'nebular.fits'
   if (keyword_set(fitflux1)) then fitflux = strlowcase(fitflux1) $
    else fitflux = 'synflux'
   if (keyword_set(lwidth1)) then lwidth = lwidth1 $
    else lwidth = 0.002
   if (fitflux NE 'synflux' AND fitflux NE 'lineflux') then $
    message, 'Invalid string for FITFLUX'

   cspeed = 2.99792458e5

   res_all = 0

   csize = 1.6
   select_tags = ['PLATE','MJD','FIBERID','OBJTYPE','PLUG_RA','PLUG_DEC', $
    'CLASS','SUBCLASS','Z','ZWARNING']
   res_append = create_struct('SFD_EBV', 0.)

   if (keyword_set(logfile)) then begin
      splog, filename=logfile
      splog, 'Log file opened ' + systime()
   endif

   ;----------
   ; Read the sky PCA components

   skypcaflux = mrdfits(skyfile)
   skyloglam = mrdfits(skyfile, 1)
   if (NOT keyword_set(skypcaflux)*keyword_set(skyloglam)) then $
    message, 'Unable to read sky file ' + skyfile
   if (keyword_set(npoly)) then $
    skypcaflux = [[skypcaflux], [poly_array(n_elements(skyloglam),npoly)]]

   ;----------
   ; Read line lists and convert to vacuum (used for masking lines
   ; in galaxy spectra)

   linefile = filepath('emlines.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
   yanny_read, linefile, pdata
   linelist = *pdata[0]
   yanny_free, pdata

   vaclambda = linelist.lambda
   airtovac, vaclambda
   linelist.lambda = vaclambda

   ;----------
   ; Generate the line list that will be used for fitting the sky nebular
   ; lines, where we force everything to have the same velocity, etc.,
   ; except for all the [O I] lines which can be atmospheric.

   if (keyword_set(lambda1)) then begin
      lambda = lambda1
      q_oxygen = 0
   endif else begin
      indx = where(linelist.lambda GT 3800.)
      lambda = linelist[indx].lambda
      linename = linelist[indx].name
      q_oxygen = strmid(linelist[indx].name,0,5) EQ '[O_I]'
   endelse
   nline = n_elements(lambda)
   zindex = lonarr(n_elements(lambda)) + q_oxygen
   windex = lonarr(n_elements(lambda)) + q_oxygen

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
      bterms_all = 0
      readspec, plist[iplate].plate, mjd=plist[iplate].mjd, $
       zans=zans, plug=plug, /silent

      if (NOT keyword_set(logfile)) then begin
         splog, filename='nebular-'+string(format='(i4.4,"-",i5.5)',plist[iplate].plate,zans[0].mjd)+'.log'
         splog, 'Log file opened ' + systime()
      endif
      splog, 'IDL version: ' + string(!version,format='(99(a," "))')
      spawn, 'uname -a', uname
      splog, 'UNAME: ' + uname[0]
      splog, 'DISPLAY=' + getenv('DISPLAY')
      splog, 'BOSS_SPECTRO_REDUX=' + getenv('BOSS_SPECTRO_REDUX')
      splog, 'idlspec2d version ' + idlspec2d_version()
      splog, 'idlutils version ' + idlutils_version()
      splog, ''
      splog, 'SKYFILE= ', skyfile
      splog, 'FITRANGE=', (keyword_set(fitrange) ? fitrange : '')
      splog, 'ZLIMITS=', (keyword_set(zlimits) ? zlimits : '')
      splog, 'SIGLIMITS=', (keyword_set(siglimits) ? siglimits : '')
      splog, 'NPOLY=', (keyword_set(siglimits) ? npoly : '')
      splog, 'NPOLY=', (n_elements(siglimits) NE 0 ? npoly : '')
      splog, 'FITFLUX= ', (keyword_set(fitflux) ? fitflux : '')
      splog, 'LWIDTH=', (keyword_set(lwidth) ? lwidth : '')
      splog, 'LWIDTH=', (keyword_set(lwidth) ? lwidth : '')
      splog, 'OUTFILE= ', outfile
      splog, 'ONLYSKY=', (keyword_set(onlysky) ? 1 : 0)

      t0 = systime(1)
      if (keyword_set(zans[0])) then begin
         plist[iplate].mjd = zans[0].mjd ; Fill in if this was zero
         zans_trim = struct_selecttags(zans, select_tags=select_tags)

         qsky = (zans.zwarning AND 1) NE 0 AND (zans.zwarning AND 2^1+2^7) EQ 0
         qgalaxy = zans.zwarning EQ 0 AND strmatch(zans.class,'GALAXY*')
         qstar = zans.zwarning EQ 0 AND strmatch(zans.class,'STAR*')
         if (keyword_set(onlysky)) then qgood = qsky $
          else qgood = qsky OR qgalaxy OR qstar

         nfiber = n_elements(zans)
         for ifiber=0L, nfiber-1L do begin
            if (qgood[ifiber]) then begin
               splog, 'Working on PLATE= ', zans[ifiber].plate, $
                ' MJD= ', zans[ifiber].mjd, ' FIBER=', zans[ifiber].fiberid

               ; Read in all the individual exposures for this object
               ; (blue and red cameras)
               readonespec, plist[iplate].plate, mjd=plist[iplate].mjd, $
                zans[ifiber].fiberid, $
                flux=flux, invvar=invvar, mask=mask, loglam=loglam, $
                sky=sky, synflux=synflux, lineflux=lineflux, expnum=expnum, $
                framehdr=framehdr, /silent
               invvar = invvar * ((mask AND pixelmask_bits('COMBINEREJ')) EQ 0)
               if (fitflux EQ 'lineflux') then synflux = lineflux

               ;----------
               ; Generate PCA sky background basis vectors,
               ; resampled to the wavelength mapping of the data

               ndim = size(flux,/n_dimen)
               dims = size(flux,/dimens)
               npix = dims[0]
               if (ndim EQ 1) then nobs = 1 else nobs = dims[1]

               ndim = size(skypcaflux,/n_dimen)
               dims = size(skypcaflux,/dimens)
               if (ndim EQ 1) then nsky = 1 else nsky = dims[1]

               explist = expnum[uniq(expnum,sort(expnum))]
               nexp = n_elements(explist)
               background = fltarr(npix, nobs, nsky, nexp)
               for iexp=0L, nexp-1L do begin
                  ithis = where(expnum EQ explist[iexp])
                  ; Convert the sky PCA spectra from Earth rest-frame
                  ; to heliocentric.
                  heliov = sxpar(*framehdr[ithis[0]], 'HELIO_RV')
                  for isky=0L, nsky-1L do begin
                     combine1fiber, $
                      skyloglam - alog10(1.d0 + heliov/2.99792458d5), $
                      skypcaflux[*,isky], $
                      newloglam=loglam[*,ithis], newflux=thisflux
                     background[*,ithis,isky,iexp] = thisflux
                  endfor
               endfor
               background = reform(background, npix*nobs, nsky*nexp)

               ; Discard wavelengths outside of the synthetic template fits
               if (qsky[ifiber] EQ 0) then $
                invvar = invvar * (synflux NE 0)

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
               if (keyword_set(fitrange)) then $
                invvar = invvar * (loglam GE alog10(fitrange[0]) $
                 AND loglam LE alog10(fitrange[1]))

               ; Initial guesses for the em. line Z and SIGMA from
               ; previous fibers on the same plate
               if (ifiber GT 0) then k = where(qgood[0:ifiber-1]) $
                else k = -1
               if (k[0] EQ -1) then begin
                  zguess = lonarr(n_elements(lambda))
                  sigguess = lonarr(n_elements(lambda)) + 1.5d-4
               endif else begin
                  zguess = djs_median( $
                   reform(res_all[*,k].linez,nline,n_elements(k)), 2)
                  sigguess = djs_median( $
                   reform(res_all[*,k].linesigma,nline,n_elements(k)), 2) $
                   / (alog(10.)*cspeed)
               endelse
               if (keyword_set(bterms_all)) then begin
                  if (size(bterms_all,/n_dimen) EQ 1) then $
                   backguess = bterms_all $
                  else $
                   backguess = djs_median(bterms_all, 2)
                  ; Replicate this guess for all exposures...
                  backguess = reform(rebin(backguess, nsky, nexp), nsky*nexp)
               endif else begin
                  backguess = fltarr(nsky*nexp)
               endelse

               ii = where(invvar NE 0, ngpix)
               if (ngpix EQ 0) then ii = 0
               res1 = linebackfit(lambda, loglam[ii], $
                flux[ii]+sky[ii]-synflux[ii], invvar=invvar[ii], $
                linename=linename, background=background[ii,*], $
                zindex=zindex, windex=windex, $
                zguess=zguess, sigguess=sigguess, backguess=backguess, $
                zlimits=zlimits, siglimits=siglimits, $
                yfit=yfit1, bfit=bfit1, bterms=bterms1)
               bterms1 = reform(bterms1, nsky, nexp)
               bterms_all = keyword_set(bterms_all) ? $
                [[bterms_all],[bterms1]] : bterms1
               yfit = fltarr(npix,nobs)
               yfit[ii] = yfit1
               bfit = fltarr(npix,nobs)
               bfit[ii] = bfit1

               if (NOT keyword_set(res_all)) then begin
                  res_blank = create_struct(zans_trim[0], res1[0], res_append)
                  struct_assign, {junk:0}, res_blank
                  res_all = replicate(res_blank, nline, nfiber)
               endif
               index_to = ifiber*nline+lindgen(nline)
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
;xrange = [6500,6800] & yrange=[-20,100] ; ???
                  title = string(plist[iplate].plate, plist[iplate].mjd, $
                   zans[ifiber].fiberid, $
                   format='("Plate ", i4, " MJD ", i5, " Fiber ", i3)')
                  splot, [0], [0], /nodata, $
                   xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
                   charsize=csize, xtitle='Wavelength [Ang]', ytitle='Flux', $
                   title=title
                  for iline=0, n_elements(lambda)-1 do $
                   soplot, lambda[iline]+[0,0], !y.crange, color='green'
                  for iobs=0, nobs-1 do begin
                     jj = where(invvar[*,iobs] NE 0, ct)
                     if (ct GT 1) then begin
                        soplot, xplot[jj,iobs], yplot[jj,iobs]
                        soplot, xplot[jj,iobs], yfit[jj,iobs], color='red'
                        soplot, xplot[jj,iobs], yfit[jj,iobs]-bfit[jj,iobs], $
                         color='green'
                     endif
                  endfor
                  thisclass = qsky[ifiber] ? 'SKY' : zans[ifiber].class
                  xyouts, total([0.9,0.1]*!x.crange), $
                   total([0.1,0.9]*!y.crange), thisclass, charsize=csize
                  cc = get_kbrd(1) 
               endif
            endif
         endfor

         ;----------
         ; Copy ZANS data into the output structure

         for ifiber=0L, nfiber-1L do begin
            index_to = ifiber*nline+lindgen(nline)
            copy_struct_inx, replicate(zans_trim[ifiber],nline), res_all, $
             index_to=index_to
            ; Copy over LINENAME,LINEWAVE even for fibers that were not fit
            res_all[index_to].linename = res1.linename
            res_all[index_to].linewave = res1.linewave
         endfor

         ;----------
         ; Add the SFD reddening values

         euler, res_all.plug_ra, res_all.plug_dec, ll, bb, 1
         res_all.sfd_ebv = dust_getval(ll, bb, /interp)

         ;----------
         ; Create/append to the output file

         while(djs_lockfile(outfile) EQ 0) do wait, 1
         create = ((iplate EQ 0) AND (file_test(outfile) EQ 0)) $
          OR keyword_set(create)
         splog, 'Writing PLATE=', plist[iplate].plate, $
          ' MJD= ', plist[iplate].mjd, ' to file ', outfile
         mwrfits_chunks, res_all, outfile, create=create, $
          append=(create EQ 0), /silent
         djs_unlockfile, outfile

         splog, 'Elapsed time for plate ', iplate+1, ' of ', nplate, $
          ' = ', systime(1)-t0, ' sec'
      endif else begin
         splog, 'Skipping PLATE= ', plist[iplate].plate, $
          ' MJD= ', plist[iplate].mjd, ' (not found)'
      endelse
      if (NOT keyword_set(logfile)) then splog, /close
   endfor

   if (keyword_set(logfile)) then splog, /close

   return
end
;------------------------------------------------------------------------------
