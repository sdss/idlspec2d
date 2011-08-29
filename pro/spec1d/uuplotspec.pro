;+
; NAME:
;   uuplotspec
;
; PURPOSE:
;   Based on D. Schlegel's plotspec routine for plotting spectra from Princeton-1D spectro outputs,
;   a floating X-Window was added in order to provide a functional interface for navigation and control and
;   to provide feedback to an online database at http://boss.astro.utah.edu
;
;
; CALLING SEQUENCE:
;   uuplotspec, plate, [ fiberid, mjd=, znum=, zmanual=, $
;    nsmooth=, /zline, /nosyn, /noerr, $
;    /sky, /ormask, /andmask, psfile=, /restframe, $
;    /zwarning, /allexp, topdir=, run1d=, run2d=, _EXTRA= ]
;
; INPUTS:
;   plate      - Plate number(s)
;
; OPTIONAL INPUTS:
;   fiberid    - Fiber number(s); if not set, then plot all fibers for
;                each plate specified.
;   mjd        - MJD number(s); if not set, then select the most recent
;                data for each plate (largest MJD).
;   znum       - If set, then return not the best-fit redshift, but the
;                ZUM-th best-fit; e.g., set ZNUM=2 for second-best fit.
;   zmanual    - If set, then do an on-the-fly fit to either a galaxy
;                template (if ZMANUAL[0]<1) or QSO (if ZMANUAL[0]>=1),
;                overriding ZNUM.  If this is a 2-element array, then fit
;                between all redshifts in the range [ZMANUAL[0],ZMANUAL[1]].
;                Do not include any polynomial terms with the PCA templates.
;   nsmooth    - If set, then boxcar smooth both the object and synthetic
;                spectra with a width equal to NSMOOTH.
;   zline      - If set, then overplot the emission line fits.
;   nosyn      - If set, then do not overplot the synthetic fit spectrum (blue).
;   noerr      - If set, then do not overplot the error vector (red).
;   sky        - If set, then overplot the sky spectrum (green).
;   ormask     - If set, then plot the OR-mask bits in yellow crosses.
;   andmask    - If set, then plot the AND-mask bits in red squares.
;   psfile     - If set, then send plot to a PostScript file instead of
;                to the SPLOT interactive widget.  The PostScript file name
;                can be set explicitly, e.g. with PSFILE='test.ps'.  Or if
;                you simply set this as a flag, e.g. with /PSFILE, then the
;                default file name is "spec-pppp-mmmmm-fff.ps",
;                where pppp=plate number, mmmmm=MJD, fff=fiber ID.
;                If FIBERID is specified, then put all plots in a single file
;                named "spec-pppp-mmmmm.ps".
;   restframe  - If set, then plot the wavelengths in the rest frame,
;                e.g. divide the wavelengths by (1+z).
;   zwarning   - If set, then only select those non-sky fibers where the
;                ZWARNING flag has been set; can be used with or without
;                specifying fiber numbers with FIBERID.
;   allexp     - If set, then plot all the individual exposure spectra,
;                rather than the co-added spectrum.
;   topdir     - TOPDIR; if not set, then default to BOSS_SPECTRO_REDUX
;   run1d      - RUN1D; if not set, then default to RUN1D
;   run2d      - RUN2D; if not set, then default to RUN2D
;   _EXTRA     - Keywords for SPLOT and XYOUTS, such as XRANGE, YRANGE, THICK,
;                or keywords for READSPEC and READONESPEC such as TOPDIR,
;                RUN2D, RUN1D.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The data are read with READSPEC.  See the documentation for that
;   routine to see how to set environment variables that describe where
;   the data files are.
;
; EXAMPLES:
;   Plot the spectrum of plate 401, fiber #100 using the SPLOT plotting tool:
;     IDL> uuplotspec, 401, 100
;
;   The spectrum is shown in white, the errors in red (except masked points
;   are set to zero), and the best-fit eigenspectrum in blue. The mouse
;   buttons will zoom in (left), recenter (center), or zoom out (right).
;   The frame can be saved as a PostScript file by selecting File->WriteEPS
;   from the left-hand corner.
;
;   Make the same plot, but boxcar-smooth the spectrum and limit the
;   wavelength range to [4000,5000] Angstroms:
;     IDL> uuplotspec, 401, 100, nsmooth=10, xrange=[5000,6000]
;
;   Some plates are observed on multiple nights. To select one of the two
;   observations of plate 306:
;     IDL> uuplotspec, 306, 20, mjd=51690
;
;   Loop through all the spectra for plate 401, interactively:
;     IDL> uuplotspec, 401
;
;   Plot all the spectra from plate 401 to a single PostScript file:
;     IDL> uuplotspec, 401, /psfile
;
;   Plot all the spectra from plate 401 to 640 individual PostScript files:
;     IDL> uuplotspec, 401, lindgen(640)+1, /psfile
;
;   Plot a list of 3 objects, each with its own plate, MJD, and fiberid:
;     IDL> plate = [400,400,401]
;     IDL> mjd = [51820,51820,51788]
;     IDL> fiberid = [10,11,20]
;     IDL> uuplotspec, plate, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, fiberid
;
; BUGS:
;   If the user interactively rescales in Y, then the labels for ORMASK
;   and ANDMASK are no longer lined up vertically with the bit mask plot.
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_icolor()
;   djs_oplot
;   djs_plot
;   djs_xyouts
;   plotspec_image
;   readspec
;   soplot
;   splot
;   sdss_flagname()
;   sxyouts
;   synthspec()
;   textoidl()
;
; INTERNAL SUPPORT ROUTINES (plotspec):
;   plotspec_mask
;   plotspec1
;
; REVISION HISTORY (plotspec):
;   01-Sep-2000  Written by D. Schlegel, Princeton
;
; INTERNAL SUPPORT ROUTINES (uuplotspec):
;   uuplotspec_mask
;   uuplotspec1
;   uuplotspec_init (modified from plotspec)
;  FLOATING WINDOW BASE ROUTINES
;   uuPlotspecBase
;   uuPlotspecBase_event
;   uuPlotspecBase_refresh
;   uuLogin
;   uuLogin_event
;  DATABASE BRIDGE
;   uuDatabase_webget
;   uuDatabase_download
;   uuDatabase_member
;   uuDatabase_comment
;   uuDatabase_recentcommentlist
;   uuDatabase_make_yanny
;   uuDatabase_post
;   uuDatabase_select
;  UTILITY FUNCTIONS
;  is_numeric
;  is_integer
;  undefine
;
;  REVISION HISTORY (uuplotspec):
;   14-May-2010   Written by Joel R. Brownstein, University of Utah
;-
;------------------------------------------------------------------------------
pro uuplotspec_mask, wave, thismask, psfile=psfile, nolabel=nolabel, $
    _EXTRA=Extra
  ;==============================================================================
  ; Plotspec procedure:  identical to plotspec_mask
  ;==============================================================================
    
  bitlabel = sdss_flagname('SPPIXMASK', 2UL^32-1, /silent)
  bitnum = where(bitlabel NE '', nlabel)
  bitlabel = bitlabel[bitnum]
  
  for ilabel=0, nlabel-1 do begin
    ypos = ( (ilabel+1) * !y.crange[1] + (nlabel-ilabel) * !y.crange[0] ) $
      / (nlabel+1)
    ynorm = ( (ilabel+0.8) * !y.window[1] + (nlabel-ilabel) * !y.window[0] ) $
      / float(nlabel+1)
    ipix = where((thismask AND 2L^bitnum[ilabel]) NE 0, npix)
    if (npix GT 0) then begin
      if (keyword_set(psfile)) then begin
        djs_oplot, [wave[ipix]], [replicate(ypos,npix)], $
          _EXTRA=Extra
      endif else begin
        soplot, [wave[ipix]], [replicate(ypos,npix)], $
          _EXTRA=Extra
      endelse
    endif
    if (NOT keyword_set(nolabel)) then begin
      if (keyword_set(psfile)) then begin
        djs_xyouts, !x.window[1], ynorm, bitlabel[ilabel]+' ', $
          charsize=1.4, align=1.0, _EXTRA=Extra, /normal
      endif else begin
        sxyouts, !x.window[1], ynorm, bitlabel[ilabel]+' ', $
          charsize=1.4, align=1.0, _EXTRA=Extra, /normal
      endelse
    endif
  endfor
  
  return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuplotspec1, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, $
    psfile=psfile, xrange=passxr, yrange=passyr, $
    allexp=allexp, _EXTRA=Extra
  ;==============================================================================
  ; Plotspec procedure:  modified from plotspec1 to accept keywords from the
  ; uuplotspec floating window event handler
  ;==============================================================================
  common plotspec_state, platelist, fiberidlist, mjdlist, topdirlist, run1dlist, run2dlist, ifiber, keyword, keywordset, uumessage
  common uuPlotspecBase_state, uuState, recentcommentlist
  
  cspeed = 2.99792458e5
  textcolor = 'green'
  linecolor = 'magenta'
  orcolor = 'yellow'
  andcolor = 'red'
 
  readspec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, znum=keyword.znum, flux=objflux, $
    wave=wave, plug=plug, zans=zans, _EXTRA=Extra, /silent
  if (NOT keyword_set(objflux)) then begin
    uumessage = "Spectrum not found for plate="+strtrim(plate,2)+", MJD="+strtrim(mjd,2)+", fiberID="+strtrim(fiberid,2) & print, uumessage
    return
  endif
  if (keywordset.allexp) then begin
    readonespec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, wave=allwave, flux=allflux, $
      _EXTRA=Extra, /silent
    if (n_elements(allflux) gt 1) then begin
      ndim = size(allflux,/n_dimen)
      if (ndim EQ 1) then nexp = 1 $
      else nexp = (size(allflux,/dimens))[1]
    endif else begin
      keywordset.allexp=0
      widget_control, uuState.uukeywordsid, get_value=uukeywordset
      uukeywordset[6] = keywordset.allexp
      widget_control, uuState.uukeywordsid, set_value=uukeywordset
      uumessage = 'File not found: spCFrame*.fits'
    endelse
  endif
  if (keywordset.restframe) then begin
    wave = wave / (1. + zans.z)
    if (keyword_set(allwave)) then allwave = allwave / (1. + zans.z)
  endif
  if (NOT keywordset.noerr) then $
    readspec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, flerr=objerr, _EXTRA=Extra, /silent
  if (keywordset.zmanual) then begin
    readspec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, invvar=objivar, loglam=loglam, $
      objhdr=hdr, _EXTRA=Extra, /silent
    if (keyword.zmanual[0] LT 1.) then eigenfile = 'spEigenGal-*.fits' $
    else eigenfile = 'spEigenQSO-*.fits'
    npoly = 0
    zrange=keyword.zmanual
    
    res_manual = zfind(objflux, objivar, hdr=hdr, $
      eigenfile=eigenfile, npoly=npoly, zmin=zrange[0], zmax=zrange[1], $
      pspace=1, nfind=1, width=1)
    keyword.manualz = res_manual.z
    synflux = synthspec(res_manual, loglam=loglam)
  endif else begin
    if (NOT keywordset.nosyn) then $
      readspec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, znum=keyword.znum, synflux=synflux, $
      _EXTRA=Extra, /silent
  endelse
  if (keywordset.sky) then $
    readspec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, sky=sky, _EXTRA=Extra, /silent
  if (keywordset.ormask) then $
    readspec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, ormask=ormask, _EXTRA=Extra, /silent
  if (keywordset.andmask) then $
    readspec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, andmask=andmask, _EXTRA=Extra, /silent
  if (keyword_set(zans) AND keywordset.zline) then $
    readspec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, zline=zline, lineflux=lineflux, $
    _EXTRA=Extra, /silent
    
  if (keywordset.nsmooth) then nsmooth = keyword.nsmooth else nsmooth = 1
  
  if (nsmooth GT 1) then begin
    objflux = smooth(objflux, nsmooth)
    if (keyword_set(allflux)) then begin
      for iexp=0, nexp-1 do $
        allflux[*,iexp] = smooth(allflux[*,iexp], nsmooth)
    endif
    if (keyword_set(synflux)) then $
      synflux = smooth(synflux, nsmooth)
    if (keyword_set(lineflux)) then $
      lineflux = smooth(lineflux, nsmooth)
  endif
  
  targstring = strmatch(plug.objtype,'SKY*') ? 'SKY ' : ''
  if (tag_exist(plug,'PRIMTARGET')) then $
    targstring += sdss_flagname('TARGET', plug.primtarget, /concat)+' '
  if (tag_exist(plug,'SECTARGET')) then $
    targstring += sdss_flagname('TTARGET', plug.sectarget, /concat)+' '
  if (tag_exist(plug,'BOSS_TARGET1')) then $
    targstring += sdss_flagname('BOSS_TARGET1', plug.boss_target1, /concat)+' '
  if (tag_exist(plug,'BOSS_TARGET2')) then $
    targstring += sdss_flagname('BOSS_TARGET2', plug.boss_target2, /concat)+' '
  if (tag_exist(plug,'ANCILLARY_TARGET1')) then targstring += $
    sdss_flagname('ANCILLARY_TARGET1', plug.ancillary_target1, /concat)+' '
  if (tag_exist(plug,'ANCILLARY_TARGET2')) then targstring += $
    sdss_flagname('ANCILLARY_TARGET2', plug.ancillary_target2, /concat)+' '
  targstring = strtrim(targstring) ; get rid of trailing spaces
  
  csize = 1.75
  if (keyword_set(passyr)) then begin
    yrange = passyr
    ymin = yrange[0]
    ymax = yrange[1]
  endif else begin
    if (keyword_set(synflux)) then $
      yrange = minmax(synflux) $
    else $
      yrange = minmax(objflux)
      
    if (yrange[0] EQ yrange[1]) then yrange = minmax(objflux)
    ymin = (1.3 * yrange[0] - 0.3 * yrange[1]) < 0
    ymax = -0.3 * yrange[0] + 1.3 * yrange[1]
    if (ymax EQ ymin) then ymax = ymin + 1
    yrange = [ymin, ymax]
  endelse
  if (keyword_set(passxr)) then xrange = passxr $
  else xrange = minmax(wave)
  if (keywordset.ormask OR keywordset.andmask) then $
    xrange[1] = 1.15 * xrange[1] - 0.15 * xrange[0]
    
  title = 'Plate ' + strtrim(string(plate),2) $
    + '  Fiber ' + strtrim(string(fiberid),2) $
    + '  MJD=' + strtrim(string(mjd),2)
  if (keywordset.restframe) then xtitle = 'Rest-Frame Wavelength [Ang]' $
  else xtitle = 'Observed Wavelength [Ang]'
  if (keyword_set(psfile)) then begin
    djs_plot, xrange, yrange, /nodata, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=TeXtoIDL('Flux [10^{-17} erg/s/cm^2/Ang]'), $
      title=title, charsize=csize, _EXTRA=KeywordsForSplot, /xstyle, /ystyle
    if (keywordset.allexp) then begin
      for iexp=0, nexp-1 do $
        djs_oplot, allwave[*,iexp], allflux[*,iexp], $
        _EXTRA=KeywordsForSplot
    endif else begin
      djs_oplot, wave, objflux, _EXTRA=KeywordsForSplot
    endelse
    if (NOT keywordset.noerr) then $
      djs_oplot, wave, objerr, color='red', _EXTRA=KeywordsForSplot
    if (keywordset.sky) then $
      djs_oplot, wave, sky, color='green', lw=2, _EXTRA=KeywordsForSplot
    if (keyword_set(synflux)) then $
      djs_oplot, wave, synflux, color='blue', lw=2, _EXTRA=KeywordsForSplot
  endif else begin
    if (NOT keyword_set(noerase)) then $
      splot, xrange, yrange, /nodata, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=TeXtoIDL('Flux [10^{-17} erg/s/cm^2/Ang]'), $
      title=title, charsize=csize, _EXTRA=KeywordsForSplot
    if (keywordset.allexp) then begin
      for iexp=0, nexp-1 do $
        soplot, allwave[*,iexp], allflux[*,iexp], _EXTRA=KeywordsForSplot
    endif else begin
      soplot, wave, objflux, _EXTRA=KeywordsForSplot
    endelse
    if (NOT keywordset.noerr) then $
      soplot, wave, objerr, color='red', _EXTRA=KeywordsForSplot
    if (keywordset.sky) then $
      soplot, wave, sky, color='green', lw=2, _EXTRA=KeywordsForSplot
    if (keyword_set(synflux)) then $
      soplot, wave, synflux, color='blue', lw=2, _EXTRA=KeywordsForSplot
  endelse
  
  xpos = 0.9 * !x.window[0] + 0.1 * !x.window[1]
  dypos = 0.05 * (!y.window[0] - !y.window[1])
  ypos = !y.window[1] + 1.5 * dypos
  
  if ((keywordset.zmanual) OR keyword_set(zans)) then begin
    keyword.manualclass = strtrim(zans.class,2)
    if (keywordset.zmanual) then begin
      zstring = res_manual.tfile $
        + '  z=' + string(res_manual.z,format='(f8.5)')
    endif else begin
      keyword.manualz=zans.z
      zstring = zans.class + ' ' + zans.subclass
      cz = zans.z * cspeed
      if (abs(cz) LT 3000) then $
        zstring += '  cz=' + string(cz,format='(f6.0)') + ' km/s' $
      else $
        zstring += '  z=' + string(zans.z,format='(f8.5)')
      if (zans.zwarning NE 0) then $
        zstring +=  ' (' $
        + sdss_flagname('ZWARNING', zans.zwarning, /concat) + ')'
      if (keywordset.znum) then $
        zstring += ' (fit #' + strtrim(string(keyword.znum),2) + ')'
    endelse
    
    if (keyword_set(psfile)) then $
      xyouts, xpos, ypos, zstring, $
      charsize=csize, color=djs_icolor(textcolor), /normal, $
      _EXTRA=KeywordsForSplot $
    else $
      sxyouts, xpos, ypos, zstring, $
      charsize=csize, color=textcolor, /normal, $
      _EXTRA=KeywordsForSplot
      
    ypos = ypos + dypos
    
    if (keyword_set(res_manual)) then thisrchi2 = res_manual.rchi2 $
    else thisrchi2 = zans.rchi2
    if (keyword_set(psfile)) then $
      xyouts, xpos, ypos, $
      TeXtoIDL('X^2_r =' + strtrim(string(thisrchi2, format='(f7.3)'),2)), $
      charsize=csize, color=djs_icolor(textcolor), /normal, $
      _EXTRA=KeywordsForSplot $
    else $
      sxyouts, xpos, ypos, $
      TeXtoIDL('X^2_r =' + strtrim(string(thisrchi2, format='(f7.3)'),2)), $
      charsize=csize, color=textcolor, /normal, $
      _EXTRA=KeywordsForSplot
  endif
  
  if (keyword_set(lineflux)) then begin
    if (keyword_set(psfile)) then $
      djs_oplot, wave, lineflux, color=linecolor, lw=2, _EXTRA=KeywordsForSplot $
    else $
      soplot, wave, lineflux, color=linecolor, lw=2, _EXTRA=KeywordsForSplot
      
    linewave = zline.linewave $
      * (1 + zline.linez * (keywordset.restframe EQ 0))
    ; Convert line sigma from km/sec to Angstroms
    linesigma = linewave * zline.linesigma / cspeed
    linepeak = zline.linecontlevel + zline.linearea / (sqrt(2*!pi) * linesigma)
    for iline=0, n_elements(zline)-1 do begin
      if (zline[iline].linearea_err GT 0) then begin
        if (keyword_set(psfile)) then $
          xyouts, linewave[iline], linepeak[iline], $
          '  '+zline[iline].linename, orient=90, $
          charsize=0.75*csize, color=djs_icolor(linecolor), $
          _EXTRA=KeywordsForSplot $
        else $
          sxyouts, linewave[iline], linepeak[iline], $
          '  '+zline[iline].linename, orient=90, $
          charsize=0.75*csize, color=linecolor, $
          _EXTRA=KeywordsForSplot
      endif
    endfor
  endif
  
  if (keywordset.ormask) then begin
    uuplotspec_mask, wave, ormask, psfile=psfile, $
      psym=1, symsize=0.6, color=orcolor, $
      nolabel=keywordset.andmask, _EXTRA=KeywordsForSplot
  endif
  
  if (keywordset.andmask) then begin
    uuplotspec_mask, wave, andmask, psfile=psfile, $
      psym=6, symsize=0.6, color=andcolor, _EXTRA=KeywordsForSplot
  endif
  
  if (keyword_set(targstring)) then begin
    ypos = ypos + dypos
    if (keyword_set(psfile)) then $
      xyouts, xpos, ypos, 'Target = '+targstring, $
      charsize=csize, color=djs_icolor(textcolor), /normal, $
      _EXTRA=KeywordsForSplot $
    else $
      sxyouts, xpos, ypos, 'Target = '+targstring, $
      charsize=csize, color=textcolor, /normal, $
      _EXTRA=KeywordsForSplot
  endif
  
  return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuplotspec_init, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, $
    psfile=psfile, xrange=xrange, yrange=yrange, $
    _EXTRA=Extra
  ;==============================================================================
  ; Plotspec procedure:  modified from plotspec to accept keywords from the
  ; uuplotspec floating window event handler
  ;==============================================================================
    
  common plotspec_state, platelist, fiberidlist, mjdlist, topdirlist, run1dlist, run2dlist, ifiber, keyword, keywordset, uumessage
  
  quiet = !quiet
  !quiet = 1
  ;----------
  ; If MJD is not set, then find the MJD for each plate
  
  nplate = n_elements(plate)
  if (NOT keyword_set(mjd)) then begin
    mjd = lonarr(nplate)
    for iplate=0, nplate-1 do begin
      mjd1 = 0
      readspec, plate[iplate], mjd=mjd1, topdir=topdir[iplate], run1d=run1d[iplate], run2d=run2d[iplate], _EXTRA=Extra, /silent
      if (NOT keyword_set(mjd1)) then begin
        uumessage = 'No MJD found for plate '+strtrim(plate[iplate],2) & print, uumessage
        !quiet = quiet
        return
      endif
      mjd[iplate] = mjd1
    endfor
  endif else begin
    if (n_elements(mjd) NE nplate) then begin
      uumessage = 'Number of elements in PLATE and MJD do not agree' & print, uumessage
      !quiet = quiet
      return
    endif
  endelse
  
  ;----------
  ; If /ZWARNING is set, then find the flagged fibers to plot.
  
  if (keywordset.zwarning) then begin
    if (keyword_set(fiberid)) then begin
      uumessage = 'FIBERID and /ZWARNING cannot both be set.' & print, uumessage
      !quiet = quiet
      return
    endif
    
    for iplate=0L, nplate-1L do begin
      readspec, plate[iplate], mjd=mjd[iplate], topdir=topdir[iplate], run1d=run1d[iplate], run2d=run2d[iplate], $
        _EXTRA=Extra, /silent, zans=zans
      if (NOT keyword_set(zans)) then begin
        uumessage = 'No spZ file found for selecting ZWARNING flags' & print, uumessage
        !quiet = quiet
        return
      endif
      indx = where((zans.zwarning AND 1) EQ 0 AND zans.zwarning NE 0, nthis)
      if (nthis GT 0) then begin
        if (NOT keyword_set(fiberid)) then begin
          platelist = replicate(plate[iplate], nthis)
          mjdlist = replicate(mjd[iplate], nthis)
          fiberid = zans[indx].fiberid
          topdirlist = replicate(topdir[iplate], nthis)
          run1dlist = replicate(run1d[iplate], nthis)
          run2dlist = replicate(run2d[iplate], nthis)
        endif else begin
          platelist = [platelist, replicate(plate[iplate], nthis)]
          mjdlist = [mjdlist, replicate(mjd[iplate], nthis)]
          fiberid = [fiberid, zans[indx].fiberid]
          topdirlist = [topdirlist, replicate(topdir[iplate], nthis)]
          run1dlist = [run1dlist, replicate(run1d[iplate], nthis)]
          run2dlist = [run2dlist, replicate(run2d[iplate], nthis)]
        endelse
      endif
    endfor
    if (NOT keyword_set(fiberid)) then begin
      uumessage = 'No non-sky fibers with ZWARNING flag set' & print, uumessage
      !quiet = quiet
      return
    endif
    nfiber = n_elements(fiberid)
    uumessage = 'Selecting '+strtrim(nfiber,2)+' non-sky fibers with ZWARNING flag set' & print, uumessage
  endif
  
  ;----------
  ; Set FIBERID to [1,...,NFIBER] (for each plate) if not set.
  ;
  ; If writing to a PostScript file, then all plots are in the same file
  ; either if PSFILE is that file name, or if FIBERID is not specified
  ; (and then all spectra are being plotted).
  
  if (NOT keyword_set(fiberid)) then begin
    readspec, plate, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, 0*plate+1, nfiber=nfiber_tmp, $
      _EXTRA=Extra, /silent
    nfiber_tot = long(total(nfiber_tmp))
    if (nfiber_tot EQ 0) then begin
      uumessage = 'No fibers found' & print, uumessage
      !quiet = quiet
      return
    endif
    platelist = lonarr(nfiber_tot)
    mjdlist = lonarr(nfiber_tot)
    fiberid = lonarr(nfiber_tot)
    topdirlist = strarr(nfiber_tot)
    run1dlist = strarr(nfiber_tot)
    run2dlist = strarr(nfiber_tot)
    j = 0L
    for iplate=0L, nplate-1L do begin
      platelist[j:j+nfiber_tmp[iplate]-1] = plate[iplate]
      mjdlist[j:j+nfiber_tmp[iplate]-1] = mjd[iplate]
      fiberid[j:j+nfiber_tmp[iplate]-1] = lindgen(nfiber_tmp[iplate]) + 1
      topdirlist[j:j+nfiber_tmp[iplate]-1] = topdir[iplate]
      run1dlist[j:j+nfiber_tmp[iplate]-1] = run1d[iplate]
      run2dlist[j:j+nfiber_tmp[iplate]-1] = run2d[iplate]
      j += nfiber_tmp[iplate]
    endfor
    if (keyword_set(psfile)) then begin
      q_onefile = 1
      psfilename = string(plate[0], mjd[0], $
        format='("spec-",i4.4,"-",i5.5,".ps")')
    endif
  endif else begin
    if (n_elements(platelist) ne 0) then undefine, platelist
    if (n_elements(mjdlist) ne 0) then undefine, mjdlist
    if (n_elements(fiberidlist) ne 0) then undefine, fiberidlist
    if (n_elements(topdirlist) ne 0) then undefine, topdirlist
    if (n_elements(run1dlist) ne 0) then undefine, run1dlist
    if (n_elements(run2dlist) ne 0) then undefine, run2dlist
  endelse
  
  ;   if (min(fiberid) LT 1 OR max(fiberid) GT 640) then begin
  ;      uumessage = 'Invalid FIBERID (must be between 1 and 640)' & print, uumessage
  ;      return
  ;   endif
  
  ;----------
  ; If FIBERID is specified, and writing to a PostScript file,
  ; then open only one PS file for all plots.
  
  nfiber = n_elements(fiberid)
  if (size(psfile,/tname) EQ 'STRING' AND nfiber GT 1) then begin
    psfilename = psfile
    q_onefile = 1
  endif
  
  ;----------
  ; If /ZWARNING is not set, then construct the PLATELIST,MJDLIST.
  
  if (n_elements(platelist) EQ 0) then begin
    if (nplate EQ 1) then begin
      platelist = replicate(plate, nfiber)
      mjdlist = replicate(mjd, nfiber)
      topdirlist = replicate(topdir, nfiber)
      run1dlist = replicate(run1d, nfiber)
      run2dlist = replicate(run2d, nfiber)
    endif else begin
      platelist = plate
      mjdlist = mjd
      topdirlist = topdir
      run1dlist = run1d
      run2dlist = run2d
    endelse
  endif else if (n_elements(platelist) NE n_elements(fiberid)) then begin
    uumessage = 'Number of elements in PLATE and FIBERID do not agree.' & print, uumessage
    !quiet = quiet
    return
  endif else if (n_elements(platelist) NE n_elements(topdirlist)) then begin
    uumessage = 'Number of elements in PLATE and TOPDIR do not agree.' & print, uumessage
    !quiet = quiet
    return
  endif else if (n_elements(platelist) NE n_elements(run1dlist)) then begin
    uumessage = 'Number of elements in PLATE and RUN1D do not agree.' & print, uumessage
    !quiet = quiet
    return
  endif else if (n_elements(platelist) NE n_elements(run2dlist)) then begin
    uumessage = 'Number of elements in PLATE and RUN2D do not agree.' & print, uumessage
    !quiet = quiet
    return
  endif
  
  ;----------
  ; Loop over each plot if creating postscript files.
  ; otherwise navigate to each plate-mjd-fiberid identified by ifiber,
  ; starting with first in the platelist (ifiber=0)
  
  
  if (keyword_set(psfile)) then begin
  
    ifiber = 0L
    while (ifiber LT nfiber) do begin
    
      ;----------
      ; Open the PostScript file if appropriate
    
      if (NOT keyword_set(q_onefile)) then begin
        if (size(psfile,/tname) EQ 'STRING') then $
          psfilename = psfile $
        else $
          psfilename = string(platelist[ifiber], mjdlist[ifiber], $
          fiberid[ifiber], format='("spec-",i4.4,"-",i5.5,"-",i3.3,".ps")')
      endif
      
      if (NOT keyword_set(q_onefile) OR ifiber EQ 0) then begin
        dfpsplot, psfilename, /color, /square
      endif
      
      uuplotspec1, platelist[ifiber], fiberid[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber], $
        psfile=psfile, $
        xrange=xrange, yrange=yrange, $
        allexp=allexp, _EXTRA=Extra
        
      if (NOT keyword_set(q_onefile) OR ifiber EQ nfiber-1) then dfpsclose
      ifiber = ifiber + 1
    endwhile
  endif else begin
    ifiber = 0
    fiberidlist = fiberid
    uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber], $
      psfile=psfile, $
      xrange=xrange, yrange=yrange, $
      allexp=allexp, _EXTRA=Extra
      
  endelse
  
  !quiet = quiet
  return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function uuDatabase_webget, url,POST=post
  ;==============================================================================
  ; wrapper for webget that catches errors in case user is offline/site is down.
  ;==============================================================================
  servers = n_elements(url)
  server = 0
  catch, webget_error
  if (webget_error ne 0) then begin
    print, "unable to connect to host at http://boss.astro.utah.edu.  Switching to backup URL"
    server = server+1
    if (server eq servers) then begin
      catch,/cancel
      response =  {Text:'NULL'}
    endif else response = webget(url[server],POST=post,/silent)
  endif
  if (server lt servers) then begin
    response = webget(url[server],POST=post,/silent)
    ;help, response, /struc
    ;print, response.Text
    if (response.Text eq '') then begin
      if (server lt servers-1) then begin
        print, "unable to get response from host at http://boss.astro.utah.edu.  Switching to backup URL"
        response = webget(url[server+1],POST=post,/silent)
      endif
    endif
  endif else return, {Text:'NULL'}
  return, response
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function is_numeric,input
  ;==============================================================================
  ; utility function:  check input is a numeric value
  ;==============================================================================
  on_ioerror, false
  test = (double(input) eq double(input))
  return, test
  false: return, 0
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function is_integer,input
  ;==============================================================================
  ; utility function:  check input is a integer value
  ;==============================================================================
  on_ioerror, false
  test = (long(input) eq double(input))
  return, test
  false: return, 0
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro undefine, key
  ;==============================================================================
  ; utility function:  refresh common block variables when called.
  ;==============================================================================
  on_error, 1
  IF n_params() eq 0 then message, 'One argument required in call to UNDEFINE'
  tempvar = size(temporary(key))
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME:
;  SETINTERSECTION
;
; PURPOSE:
;
;   This function is used to find the intersection between two sets of integers.
;
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;
;   Utilities
;
; CALLING SEQUENCE:
;
;   intersection = SetIntersection(set_a, set_b)
;
; RETURN VALUE:
;
;   intersection:  A vector of values that are found in both set_a and set_b.
;
; ARGUMENTS:
;
;   set_a:         A vector of integers.
;
;   set_b:         A vector of integers.
;
; KEYWORDRS:
;
;  NORESULT:       Set this keyword to a value that will be returned from the function
;                  if no intersection between the two sets of numbers is found. By default, -1.
;
;  SUCCESS:        An output keyword that is set to 1 if an intersection was found, and to 0 otherwise.
;
; EXAMPLE:
;
;  IDL> set_a = [1,2,3,4,5]
;  IDL> set_b = [4,5,6,7,8,9,10,11]
;  IDL> Print, SetIntersection(set_a, set_b)
;          4   5
;
;  See http://www.dfanning.com/tips/set_operations.html for other types of set operations.
;
; NOTES:
;
;  If you read the Set Operations article pointed to above, you will see quite a lot of
;  discussion about what kinds of algorithms are faster than others. The Histogram
;  algorithms implemented here are sometimes NOT the fastest algorithms, especially
;  for sparse arrays. If this is a concern in your application, please be sure to read
;  that article.
;
; MODIFICATION HISTORY:
;
;  Written by: David W. Fanning, October 31, 2009, from code originally supplied to the IDL
;     newsgroup by Research Systems software engineers.
;  Yikes, bug in original code only allowed positive integers. Fixed now. 2 Nov 2009. DWF.
;  Fixed a problem when one or both of the sets was a scalar value. 18 Nov 2009. DWF.
;
;******************************************************************************************;
;  Copyright (c) 2009, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
FUNCTION SetIntersection, set_a, set_b, $
    NORESULT=noresult, $
    SUCCESS=success
    
  Compile_Opt StrictArr, DefInt32
  
  ; Set up noresult value.
  IF N_Elements(noresult) EQ 0 THEN noresult = -1
  
  ; Error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
    void = Error_Message()
    success = 0
    RETURN, noresult
  ENDIF
  
  ; Check parameters.
  IF N_Params() NE 2 THEN Message, 'Two input parameters or sets are required.'
  
  ; The input sets must be integers.
  IF (Size(set_a, /TYPE) GT 3) AND (Size(set_a, /TYPE) LT 12) THEN $
    Message, 'Set A must be an integer array.'
  IF (Size(set_b, /TYPE) GT 3) AND (Size(set_b, /TYPE) LT 12) THEN $
    Message, 'Set B must be an integer array.'
    
  ; If either of the sets is a scalar, make it a vector.
  IF N_Elements(set_a) EQ 1 && (Size(set_a))[0] EQ 0 THEN set_a = [set_a]
  IF N_Elements(set_b) EQ 1 && (Size(set_b))[0] EQ 0 THEN set_b = [set_b]
  
  ; Assume success.
  success = 1
  
  ; Find the intersection of the ranges.
  mina = Min(set_a, Max=maxa)
  minb = Min(set_b, Max=maxb)
  minab = mina > minb
  maxab = maxa < maxb
  
  ; If the set ranges don't intersect, leave now.
  IF ((maxa LT minab) AND (minb GT maxab)) OR ((maxb LT minab) AND (mina GT maxab)) THEN BEGIN
    success = 0
    RETURN, noresult
  ENDIF
  
  ; Find the intersection.
  r = Where((Histogram(set_a, Min=minab, Max=maxab) NE 0) AND  $
    (Histogram(set_b, Min=minab, Max=maxab) NE 0), count)
    
  ; Was there an intersection? If not, leave now.
  IF count EQ 0 THEN BEGIN
    success = 0
    RETURN, noresult
  ENDIF
  
  ; Here is the result.
  result = Temporary(r) + minab
  
  ; Return the result. Make sure to return scalar if only a single element.
  IF N_Elements(result) EQ 1 THEN RETURN, result[0] ELSE RETURN, result
  
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
Function uuLogin, Group_Leader=group_leader
  ;==============================================================================
  ; Login Function: create a modal window to accept username/password
  ;==============================================================================
  common uuPlotspecBase_state, uuState, recentcommentlist
  
  if n_elements(group_leader) eq 0 then  uuLogin = Widget_Base(Title='Login', Row=3) ELSE  uuLogin = Widget_Base(Title='Login', Row=4, /Modal, Group_Leader=group_leader)
  uuLoginrow1 = widget_base(uuLogin,/row)
  uuLoginrow2 = widget_base(uuLogin,/row)
  uuLoginrow3 = widget_base(uuLogin,Col=2)
  uuLoginrow3col1 = widget_base(uuLoginrow3)
  uuLoginrow3col2 = widget_base(uuLoginrow3,/nonexclusive)
  uuLoginrow4 = widget_base(uuLogin,Col=2)
  uuLoginrow4col1 = widget_base(uuLoginrow4)
  uuLoginrow4col2 = widget_base(uuLoginrow4)
  
  uuState.usernameid = cw_field(uuLoginrow1,/row,title = 'Username:',value = '', /string,/return_events,uvalue='username')
  label = widget_label(uuLoginrow2, Value='Password:',uvalue='label')
  password = widget_text(uuLoginrow2, Scr_XSize=142, All_Events=1, Editable=0,uvalue='password')
  void = widget_label(uuLoginrow3col1, Value='        ')
  uuState.remoteid = widget_button(uuLoginrow3col2,value='Remember Login',uvalue='remote')
  void = widget_label(uuLoginrow4col1, Value='        ')
  submitloginbutton = widget_button(uuLoginrow4col2,value='Login',uvalue='submitlogin',xsize=80,ysize=36)
  device, get_screen_size=screenSize
  xCenter = screenSize[0] / 2
  yCenter = screenSize[1] / 2
  geom = widget_info(uuLogin, /Geometry)
  xHalfSize = geom.Scr_XSize / 2
  yHalfSize = geom.Scr_YSize / 2
  widget_control, uuLogin, XOffset = xCenter-xHalfSize, YOffset = yCenter-yHalfSize
  widget_control, uuLogin, /Realize
  ptr = ptr_New("")
  widget_control, uuLogin, set_uvalue=ptr
  xmanager, 'uuLogin', uuLogin
  uuState.password = *ptr
  ptr_Free, ptr
  return, ((uuState.username ne '') and (uuState.password ne ''))
  
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO uuLogin_event, event
  ;==============================================================================
  ; Login Function: login window event handler to validate username/password
  ;==============================================================================
  common uuPlotspecBase_state, uuState, recentcommentlist
  
  widget_control, uuState.usernameid, get_value=theusername
  uuState.username=theusername
  widget_control, event.id, get_uvalue=uvalue
  if (uvalue eq 'remote') then begin
    uuState.remote = strtrim(event.select,2)
  endif else if (uvalue eq 'submitlogin') then begin
    uuState.remote = strtrim(event.select,2)
    widget_control, event.top, Get_UValue=ptr
    if (uuState.username ne '') and (*ptr ne '') then begin
      widget_control, event.top, /Destroy
      return
    endif
  endif else if (uvalue eq 'password') and (event.type eq 0) then begin
    if (event.ch eq 10) then begin
      widget_control, event.top, /Destroy
      return
    endif else begin
      widget_control, event.id, Get_Value=text
      text = text[0]
      length = StrLen(text)
      selection = Widget_Info(event.id, /Text_Select)
      widget_control, event.id, /Use_Text_Select, Set_Value='*'
      widget_control, event.id, Set_Text_Select=event.offset + 1
      widget_control, event.top, Get_UValue=ptr
      if *ptr eq "" THEN *ptr=String(event.ch) ELSE  *ptr = *ptr + String(event.ch)
      widget_control, event.top, set_uvalue=ptr
    endelse
  endif else if (uvalue eq 'password') and (event.type eq 2) then begin
    widget_control, event.top, Get_UValue=ptr
    text = *ptr
    length = StrLen(text)
    *ptr = StrMid(text, 0, length-event.length)
    passwordLen = StrLen(*ptr)
    if (passwordLen gt 0) then widget_control, event.id, Set_Value=Replicate('*', passwordLen) else widget_control, event.id, Set_Value=''
    widget_control, event.id, Set_Text_Select=event.offset
    widget_control, event.top, set_uvalue=ptr
  endif else if (uvalue eq 'username') then begin
    widget_control, event.top, /Destroy
    return
  endif
  
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuPlotspecBase
  ;==============================================================================
  ; Plotspec Base Function: create a floating window to accomodate plotspec
  ; arguments (promoted to textfields) and plotspec keywords (promoted to
  ; nonexclusive buttons).  The state of the window is maintained by the common
  ; uuState structure.
  ;==============================================================================
  common splot_state, state, graphkeys
  common plotspec_state, platelist, fiberidlist, mjdlist, topdirlist, run1dlist, run2dlist, ifiber, keyword, keywordset, uumessage
  common uuPlotspecBase_state, uuState, recentcommentlist
  if (NOT xregistered('uuplotspecbase')) then begin
    issues = ['Reduction/Calibration','Redshift/Class','Sky Subtraction','Non-masked Artifacts','Little/No Data','Other/Unknown']
    uuState = {uuplotspecbase:0L,commentheader:0L,commentbase:0L,yannybase:0L,run1d:0L,run2d:0L,loginbuttonid:0L,plateid:0L,mjdid:0L,fiberid:0L,ifiberid:0L,nfiberid:0L,usernameid:0L,username:'',remoteid:0L,remote:'0',password:'',loggedin:0,fullname:'',sid:'',messageid:0L,recentcommentid:0L,commentid:0L,comment:'',issueid:0L,issues:issues,issue:issues[0],zid:0L,zmanual0id:0L,zmanual1id:0L,uukeywordsid:0L,z:'',znumid:0L,nsmoothid:0L,classid:0L,class:'',zconfid:0,zconf:'',yannyid:0L,yanny:'spinspect',yannygroupid:0L,yannygroup:0,valid:0,action:''}
    
    recentcommentlist = [{comment:'Paste from recent comments                     ',commentid:0L}]
    if (xregistered('splot')) then begin
      uuState.uuplotspecbase = widget_base(/base_align_left,/floating,xsize=630, group_leader = state.base_id,  /column,  title='uuplotspec', uvalue = 'uuplotspecbase')
      uurow0 = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'uurow0')
      uurow1 = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'uurow1')
      uurow2 = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'uurow2')
      uurow3 = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'uurow3')
      uuState.plateid = cw_field(uurow0,/row,title = 'Plate-MJD-FiberID:',uvalue='uunav',value = '', /string,/return_events,xsize=4)
      uuState.mjdid = cw_field(uurow0,/row,title = '',uvalue='uunav',value = '', /string,/return_events,xsize=5)
      uuState.fiberid = cw_field(uurow0,/row,title = '',uvalue='uunav',value = '', /string,/return_events,xsize=4)
      uunavButton = cw_bgroup(uurow0, [' Get ', ' < ',' > '],uvalue='uunav', /ROW, SPACE=2,ysize=31)
      uuState.ifiberid = cw_field(uurow0,/row,title = '',uvalue='uuifiberid',value = '0', /string,/return_events,xsize=4)
      uuState.nfiberid = widget_label(uurow0,value = '0 of 0')
      void = widget_label(uurow0,value =' ')
      uuState.loginbuttonid = widget_button(uurow0,value = 'Logout', uvalue = 'uuLoginbutton',xsize=80,ysize=36)
      void = widget_label(uurow0,value = ' ')
      uuState.zmanual0id = cw_field(uurow1,/row,title = 'zmanual:',uvalue='uuzmanual0field',value ='', /string,/return_events,xsize=6)
      uuState.zmanual1id = cw_field(uurow1,/row,title = '',uvalue='uuzmanual1field',value ='', /string,/return_events,xsize=6)
      uuzmanualbutton = cw_bgroup(uurow1, [' z Manual '],uvalue='uuzmanualbutton',/ROW, SPACE=2,ysize=31)
      void = widget_label(uurow1,value =' ')
      uuState.znumid = cw_field(uurow1,/row,title = 'Num:',uvalue='uuznumfield',value ='', /string,/return_events,xsize=2)
      uuznumbutton = cw_bgroup(uurow1, [' z Num '],uvalue='uuznumbutton',/ROW,ysize=31)
      znumchangeid= cw_bgroup(uurow1, ['+','-'],/row,uvalue='uuznumchangebutton', space=0)
      void = widget_label(uurow1,value = ' ')
      uuState.nsmoothid = cw_field(uurow1,/row,title = 'N:',uvalue='uunsmoothfield',value ='', /string,/return_events,xsize=2)
      uunsmoothbutton = cw_bgroup(uurow1, [' N Smooth '],uvalue='uunsmoothbutton',/ROW,ysize=31)
      void = widget_label(uurow1,value = ' ')
      uukeywords = ['zline','nosyn','noerr','sky','ormask','andmask','allexp','restframe','zwarning']
      uukeywordset = [keywordset.zline,keywordset.nosyn,keywordset.noerr,keywordset.sky,keywordset.ormask,keywordset.andmask,keywordset.allexp,keywordset.restframe,keywordset.zwarning]
      uuState.uukeywordsid = cw_bgroup(uurow2,uukeywords,/row,/nonexclusive,set_value=uukeywordset,uvalue='uukeywords')
      void = widget_label(uurow2,value = ' ')
      uuState.messageid = widget_label(uurow3,uvalue='uumessagefield',value = 'Please Login, or register for an account at http://boss.astro.utah.edu', xsize=600, frame=1)
      void = widget_label(uurow3,value = ' ')
      uuState.commentheader = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'commentheader')
      uuState.commentbase = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'commentbase')
      uuState.yannybase = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'yannybase')
      uuDatabase_member, 'login'
      uuPlotspecBase_refresh
      widget_control, uuState.uuplotspecbase, /realize
      xmanager, 'uuplotspecbase', uuState.uuplotspecbase, /no_block
    endif
  endif
  return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuPlotspecBase_event, event
  ;==============================================================================
  ; Plotspec Base Function: event handler for the uuplotspec floating window.
  ;==============================================================================
  common splot_state, state, graphkeys
  common plotspec_state, platelist, fiberidlist, mjdlist, topdirlist, run1dlist, run2dlist, ifiber, keyword, keywordset, uumessage
  common uuPlotspecBase_state, uuState, recentcommentlist
  
  if (uuState.loggedin) then uumessage='' else uumessage = 'Please Login to provide spectrum feedback to http://boss.astro.utah.edu'
  widget_control, event.id, get_uvalue=uvalue
  case uvalue of
  
    'uuLoginbutton': begin
      if (uuState.loggedin) then begin
        uuState.action = 'logout'
        uuDatabase_member, uuState.action
        uumessage = 'Please Login to provide spectrum feedback to http://boss.astro.utah.edu'
        uuPlotspecBase_refresh
      endif else if uuLogin(Group_Leader=event.top) then begin
        uuState.action = 'login'
        uuDatabase_member, uuState.action
        if (uuState.loggedin) then uumessage = 'Logged in on '+SYSTIME() else uumessage = "Invalid Login.  Please try again."
        uuPlotspecBase_refresh
      endif
    end
    
    'uunav': begin
      keywordset.zmanual=0
      keyword.zmanual=[0D,0D]
      keyword.manualz = 0D
      widget_control, uuState.zmanual0id, set_value=''
      widget_control, uuState.zmanual1id, set_value=''
      keywordset.znum=0
      keyword.znum=0D
      widget_control, uuState.znumid, set_value=''
      if (uuState.loggedin) then begin
        widget_control, uuState.zid, set_value=''
        widget_control, uuState.commentid, set_value=''
      endif
      widget_control, uuState.plateid, get_value=plate
      widget_control, uuState.mjdid, get_value=mjd
      widget_control, uuState.fiberid, get_value=fiberid
      if (plate ne '') and is_integer(plate) and (n_elements(fiberidlist) eq 1) then begin
        if (mjd ne '') and is_integer(mjd) then begin
          uumessage = "Getting all fibers for this plate/mjd" & print, uumessage
          uuplotspec_init, plate, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d
        endif else begin
          uumessage = "Getting all fibers for this plate" & print, uumessage
          uuplotspec_init, plate, topdir=topdir, run1d=run1d, run2d=run2d
        endelse
        if (fiberid ne '') and is_integer(fiberid) then ifiber = long(fiberid[0])-1
      endif
      case event.value of
        1: begin ; < Button
          if (ifiber gt 0) then ifiber = ifiber - 1 else begin
            ifiber =  n_elements(fiberidlist)-1
            uumessage = "Wrapping around to show final fiber" & print, uumessage
          endelse
          uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
        end
        2: begin ; > Button
          if (ifiber lt n_elements(fiberidlist)-1) then ifiber = ifiber + 1 else begin
            ifiber = 0
            uumessage = "Wrapping around to show initial fiber" & print, uumessage
          endelse
          uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
        end
        else: begin ; Get Button or plate-mjd-fiberid field
          if (plate ne '') and is_integer(plate) and is_integer(mjd) and is_integer(fiberid) then begin
            platefibers = where(platelist eq plate[0], ncomplement=nonplatefibers)
            if (nonplatefibers ne n_elements(fiberidlist)) then begin
              if (mjd eq '' and fiberid eq '') then begin
                ifiber = platefibers[0]
              endif else if (mjd ne '' and fiberid eq '') then begin
                mjdfibers = where(mjdlist eq mjd[0], ncomplement=nonmjdfibers)
                if (nonmjdfibers ne n_elements(fiberidlist)) then begin
                  plateandmjdfibers = setIntersection(platefibers,mjdfibers)
                  ifiber = plateandmjdfibers[0]
                endif else ifiber = -1
              endif else if (mjd eq '' and fiberid ne '') then begin
                fiberidfibers = where(fiberidlist eq fiberid[0], ncomplement=nonfiberidfibers)
                if (nonfiberidfibers ne n_elements(fiberidlist)) then begin
                  plateandfiberidfibers = setIntersection(platefibers,fiberidfibers)
                  ifiber = plateandfiberidfibers[0]
                endif else ifiber = -1
              endif else if (mjd ne '' and fiberid ne '') then begin
                mjdfibers = where(mjdlist eq mjd[0], ncomplement=nonmjdfibers)
                fiberidfibers = where(fiberidlist eq fiberid[0], ncomplement=nonfiberidfibers)
                if (nonmjdfibers ne n_elements(fiberidlist)) and (nonfiberidfibers ne n_elements(fiberidlist)) then begin
                  plateandmjdandfiberidfibers = setIntersection(setIntersection(platefibers,mjdfibers),fiberidfibers)
                  ifiber = plateandmjdandfiberidfibers[0]
                endif else ifiber = -1
              endif
            endif else ifiber = -1
            if (ifiber ne -1) then begin
              uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
            endif else if (mjd ne '') then begin
              uumessage = "Getting new plate/mjd" & print, uumessage
              uuplotspec_init, plate, mjd=mjd, topdir=topdirlist[0], run1d=run1dlist[0], run2d=run2dlist[0]
            endif else begin
              uumessage = "Getting new plate" & print, uumessage
              uuplotspec_init, plate, topdir=topdirlist[0], run1d=run1dlist[0], run2d=run2dlist[0]
            endelse
          endif
          if (ifiber eq -1) then ifiber=0
        end
      endcase
      widget_control, uuState.plateid, set_value=strtrim(platelist[ifiber],2)
      widget_control, uuState.mjdid, set_value=strtrim(mjdlist[ifiber],2)
      widget_control, uuState.fiberid, set_value=strtrim(fiberidlist[ifiber],2)
      nfibertext = strtrim(ifiber+1,2)+' of '+strtrim(n_elements(fiberidlist),2)+'              '
      nfibertext = 'of '+strtrim(n_elements(fiberidlist),2)+'    '
      widget_control, uuState.ifiberid, set_value = strtrim(ifiber+1,2)
      widget_control, uuState.nfiberid, set_value = strmid(nfibertext,0,16)
      if (uuState.loggedin) then begin
        widget_control, uuState.run1d, set_value=run1dlist[ifiber]
        widget_control, uuState.run2d, set_value=run2dlist[ifiber]
      endif
    end
    
    'uuifiberid': begin
      widget_control, uuState.ifiberid, get_value=ifiber0
      if (ifiber0 ne '') and (is_integer(ifiber0)) then begin
        if (ifiber0 ge 1) and (ifiber0 le n_elements(fiberidlist)) then begin
          if (ifiber ne long(ifiber0)-1) then begin
            ifiber = long(ifiber0)-1
            keywordset.zmanual=0
            keyword.zmanual=[0D,0D]
            keyword.manualz = 0D
            widget_control, uuState.zmanual0id, set_value=''
            widget_control, uuState.zmanual1id, set_value=''
            keywordset.znum=0
            keyword.znum=0D
            widget_control, uuState.znumid, set_value=''
            if (uuState.loggedin) then begin
              widget_control, uuState.zid, set_value=''
              widget_control, uuState.commentid, set_value=''
            endif
            widget_control, uuState.plateid, set_value=strtrim(platelist[ifiber],2)
            widget_control, uuState.mjdid, set_value=strtrim(mjdlist[ifiber],2)
            widget_control, uuState.fiberid, set_value=strtrim(fiberidlist[ifiber],2)
            if (uuState.loggedin) then begin
              widget_control, uuState.run1d, set_value=run1dlist[ifiber]
              widget_control, uuState.run2d, set_value=run2dlist[ifiber]
            endif
            uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
          endif
        endif else ifiber0 = -1
      endif else ifiber0 = -1
      if (ifiber0 eq -1) then begin
        widget_control, uuState.ifiberid, set_value=strtrim(ifiber+1,2)
        uumessage = "Please enter a valid fiber index from 1 to "+strtrim(n_elements(fiberidlist),2) & print, uumessage
      endif
    end
    
    'uuzmanualbutton': begin
      keywordset.znum=0
      keyword.znum=0D
      widget_control, uuState.znumid, set_value=''
      widget_control, uuState.zmanual0id, get_value=zmanual0
      widget_control, uuState.zmanual1id, get_value=zmanual1
      zmanual0=double(zmanual0)
      zmanual1=double(zmanual1)
      if (zmanual0 ne '' and zmanual1 ne '') and (zmanual1 lt zmanual0) then begin
        zmanual01 = zmanual1
        zmanual1 = zmanual0
        zmanual0 = zmanual01
        widget_control, uuState.zmanual0id, set_value=zmanual0
        widget_control, uuState.zmanual1id, set_value=zmanual1
      endif else if (zmanual0 eq '' and zmanual1 ne '') then begin
        zmanual0=double(zmanual1)
        zmanual1=''
        widget_control, uuState.zmanual0id, set_value=zmanual0
        widget_control, uuState.zmanual1id, set_value=''
      endif
      if (zmanual0 ne '') then begin
        keywordset.zmanual=1
        if (zmanual1 ne '') then keyword.zmanual=[zmanual0,zmanual1] else keyword.zmanual=[zmanual0,zmanual0]
      endif else begin
        keywordset.zmanual=0
        keyword.zmanual=[0D,0D]
      endelse
      uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
      if (uuState.loggedin) and (zmanual0 ne '') then begin
        widget_control, uuState.zid, set_value=keyword.manualz
        if (keyword.manualclass eq 'GALAXY') then widget_control, uuState.classid, set_droplist_select=1 else if (keyword.manualclass eq 'QSO') then widget_control, uuState.classid, set_droplist_select=2 else if (keyword.manualclass eq 'STAR') then widget_control, uuState.classid, set_droplist_select=3 else widget_control, uuState.classid, set_droplist_select=0
        widget_control, uuState.issueid, set_droplist_select=1
      endif else if (uuState.loggedin) then begin
        widget_control, uuState.zid, set_value=''
        widget_control, uuState.classid, set_droplist_select=0
      endif
    end
    
    'uuzmanual0field': begin
      keywordset.znum=0
      keyword.znum=0D
      widget_control, uuState.znumid, set_value=''
      widget_control, uuState.zmanual0id, get_value=zmanual0
      widget_control, uuState.zmanual1id, get_value=zmanual1
      zmanual0=double(zmanual0)
      zmanual1=double(zmanual1)
      if (zmanual0 ne '' and zmanual1 ne '') and (zmanual1 lt zmanual0) then begin
        zmanual01 = zmanual1
        zmanual1 = zmanual0
        zmanual0 = zmanual01
        widget_control, uuState.zmanual0id, set_value=zmanual0
        widget_control, uuState.zmanual1id, set_value=zmanual1
      endif else if (zmanual0 eq '' and zmanual1 ne '') then begin
        zmanual0=double(zmanual1)
        zmanual1=''
        widget_control, uuState.zmanual0id, set_value=zmanual0
        widget_control, uuState.zmanual1id, set_value=''
      endif
      if (zmanual0 ne '') then begin
        keywordset.zmanual=1
        if (zmanual1 ne '') then keyword.zmanual=[zmanual0,zmanual1] else keyword.zmanual=[zmanual0,zmanual0]
      endif else begin
        keywordset.zmanual=0
        keyword.zmanual=[0D,0D]
      endelse
      uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
      if (uuState.loggedin) and (zmanual0 ne '') then begin
        widget_control, uuState.zid, set_value=keyword.manualz
        if (keyword.manualclass eq 'GALAXY') then widget_control, uuState.classid, set_droplist_select=1 else if (keyword.manualclass eq 'QSO') then widget_control, uuState.classid, set_droplist_select=2 else if (keyword.manualclass eq 'STAR') then widget_control, uuState.classid, set_droplist_select=3 else widget_control, uuState.classid, set_droplist_select=0
        widget_control, uuState.issueid, set_droplist_select=1
      endif else if (uuState.loggedin) then begin
        widget_control, uuState.zid, set_value=''
        widget_control, uuState.classid, set_droplist_select=0
      endif
    end
    
    'uuzmanual1field': begin
      keywordset.znum=0
      keyword.znum=0D
      widget_control, uuState.znumid, set_value=''
      widget_control, uuState.zmanual0id, get_value=zmanual0
      widget_control, uuState.zmanual1id, get_value=zmanual1
      zmanual0=double(zmanual0)
      zmanual1=double(zmanual1)
      if (zmanual0 ne '' and zmanual1 ne '') and (zmanual1 lt zmanual0) then begin
        zmanual01 = zmanual1
        zmanual1 = zmanual0
        zmanual0 = zmanual01
        widget_control, uuState.zmanual0id, set_value=zmanual0
        widget_control, uuState.zmanual1id, set_value=zmanual1
      endif else if (zmanual0 eq '' and zmanual1 ne '') then begin
        zmanual0=double(zmanual1)
        zmanual1=''
        widget_control, uuState.zmanual0id, set_value=zmanual0
        widget_control, uuState.zmanual1id, set_value=''
      endif
      if (zmanual0 ne '') then begin
        keywordset.zmanual=1
        if (zmanual1 ne '') then keyword.zmanual=[zmanual0,zmanual1] else keyword.zmanual=[zmanual0,zmanual0]
      endif else begin
        keywordset.zmanual=0
        keyword.zmanual=[0D,0D]
      endelse
      uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
      if (uuState.loggedin) and (zmanual0 ne '') then begin
        widget_control, uuState.zid, set_value=keyword.manualz
        if (keyword.manualclass eq 'GALAXY') then widget_control, uuState.classid, set_droplist_select=1 else if (keyword.manualclass eq 'QSO') then widget_control, uuState.classid, set_droplist_select=2 else if (keyword.manualclass eq 'STAR') then widget_control, uuState.classid, set_droplist_select=3 else widget_control, uuState.classid, set_droplist_select=0
        widget_control, uuState.issueid, set_droplist_select=1
      endif else if (uuState.loggedin) then begin
        widget_control, uuState.zid, set_value=''
        widget_control, uuState.classid, set_droplist_select=0
      endif
    end
    
    'uuznumbutton': begin
      keywordset.zmanual=0
      keyword.zmanual=[0D,0D]
      keyword.manualz = 0D
      widget_control, uuState.zmanual0id, set_value=''
      widget_control, uuState.zmanual1id, set_value=''
      widget_control, uuState.znumid, get_value=znum
      if (znum ne '') then begin
        keywordset.znum=1
        keyword.znum=znum
      endif else begin
        keywordset.znum=0
        keyword.znum=0L
      endelse
      uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
      if (uuState.loggedin) and (znum ne '') then begin
        widget_control, uuState.zid, set_value=keyword.manualz
        if (keyword.manualclass eq 'GALAXY') then widget_control, uuState.classid, set_droplist_select=1 else if (keyword.manualclass eq 'QSO') then widget_control, uuState.classid, set_droplist_select=2 else if (keyword.manualclass eq 'STAR') then widget_control, uuState.classid, set_droplist_select=3 else widget_control, uuState.classid, set_droplist_select=0
        widget_control, uuState.issueid, set_droplist_select=1
      endif else if (uuState.loggedin) then begin
        widget_control, uuState.zid, set_value=''
        widget_control, uuState.classid, set_droplist_select=0
      endif
    end
    
    'uuznumfield': begin
      keywordset.zmanual=0
      keyword.zmanual=[0D,0D]
      keyword.manualz = 0D
      widget_control, uuState.zmanual0id, set_value=''
      widget_control, uuState.zmanual1id, set_value=''
      widget_control, uuState.znumid, get_value=znum
      if (znum ne '') then begin
        keywordset.znum=1
        keyword.znum=znum
      endif else begin
        keywordset.znum=0
        keyword.znum=0L
      endelse
      uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
      if (uuState.loggedin) and (znum ne '') then begin
        widget_control, uuState.zid, set_value=keyword.manualz
        if (keyword.manualclass eq 'GALAXY') then widget_control, uuState.classid, set_droplist_select=1 else if (keyword.manualclass eq 'QSO') then widget_control, uuState.classid, set_droplist_select=2 else if (keyword.manualclass eq 'STAR') then widget_control, uuState.classid, set_droplist_select=3 else widget_control, uuState.classid, set_droplist_select=0
        widget_control, uuState.issueid, set_droplist_select=1
      endif else if (uuState.loggedin) then begin
        widget_control, uuState.zid, set_value=''
        widget_control, uuState.classid, set_droplist_select=0
      endif
    end
    
    'uuznumchangebutton': begin
      keywordset.zmanual=0
      keyword.zmanual=[0D,0D]
      keyword.manualz = 0D
      widget_control, uuState.zmanual0id, set_value=''
      widget_control, uuState.zmanual1id, set_value=''
      widget_control, uuState.znumid, get_value=znum
      if (znum ne '' and is_integer(znum)) then begin
        case event.value of
          0: begin
            keyword.znum=znum+1
            keywordset.znum=1
            widget_control, uuState.znumid, set_value=keyword.znum
          end
          1: begin
            if (znum gt 1) then begin
              keyword.znum=znum-1
              keywordset.znum=1
              widget_control, uuState.znumid, set_value=keyword.znum
            endif else begin
              keyword.znum=0
              keywordset.znum=0
              widget_control, uuState.znumid, set_value=''
            endelse
          end
        endcase
      endif else begin
         case event.value of
          0: begin
            keyword.znum=1
            keywordset.znum=1
            widget_control, uuState.znumid, set_value=keyword.znum
          end
          1: begin 
            keyword.znum=0
            keywordset.znum=0
            widget_control, uuState.znumid, set_value=''
          end
        endcase
      endelse
      uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
      if (uuState.loggedin) and (keyword.znum ne '') then begin
        widget_control, uuState.zid, set_value=keyword.manualz
        if (keyword.manualclass eq 'GALAXY') then widget_control, uuState.classid, set_droplist_select=1 else if (keyword.manualclass eq 'QSO') then widget_control, uuState.classid, set_droplist_select=2 else if (keyword.manualclass eq 'STAR') then widget_control, uuState.classid, set_droplist_select=3 else widget_control, uuState.classid, set_droplist_select=0
        widget_control, uuState.issueid, set_droplist_select=1
      endif else if (uuState.loggedin) then begin
        widget_control, uuState.zid, set_value=''
        widget_control, uuState.classid, set_droplist_select=0
      endif
    end

    'uunsmoothbutton': begin
      widget_control, uuState.nsmoothid, get_value=nsmooth
      if (nsmooth ne '') and (nsmooth ne '1') then begin
        keywordset.nsmooth=1
        keyword.nsmooth=nsmooth
      endif else begin
        keywordset.nsmooth=0
        keyword.nsmooth=1L
      endelse
      uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
    end
    
    'uunsmoothfield': begin
      widget_control, uuState.nsmoothid, get_value=nsmooth
      if (nsmooth ne '') and (nsmooth ne '1') then begin
        keywordset.nsmooth=1
        keyword.nsmooth=nsmooth
      endif else begin
        keywordset.nsmooth=0
        keyword.nsmooth=1L
      endelse
      uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
    end
    
    'uukeywords': begin      
      widget_control, uuState.uukeywordsid, get_value=uukeyword
      i=0
      keywordset.zline=uukeyword[i++]
      keywordset.nosyn=uukeyword[i++]
      keywordset.noerr=uukeyword[i++]
      keywordset.sky=uukeyword[i++]
      keywordset.ormask=uukeyword[i++]
      keywordset.andmask=uukeyword[i++]
      keywordset.allexp=uukeyword[i++]
      keywordset.restframe=uukeyword[i++]
      keywordset.zwarning=uukeyword[i++]
      case event.value of
        8: begin ; zwarning
          widget_control, uuState.plateid, get_value=plate
          widget_control, uuState.mjdid, get_value=mjd
          widget_control, uuState.fiberid, get_value=fiberid
          if (plate ne '') and is_integer(plate) and is_integer(mjd) and is_integer(fiberid) then begin
            if (mjd eq '' and fiberid eq '') then begin
              uuplotspec_init, plate, topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
            endif else if (mjd ne '' and fiberid eq '') then begin
              uuplotspec_init, plate, mjd=mjd, topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
            endif else if (mjd eq '' and fiberid ne '') then begin
              uuplotspec_init, plate, topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
              platefibers = where(platelist eq plate[0], ncomplement=nonplatefibers)
              if (nonplatefibers ne n_elements(fiberidlist)) then begin
                fiberidfibers = where(fiberidlist eq fiberid[0], ncomplement=nonfiberidfibers)
                if (nonfiberidfibers ne n_elements(fiberidlist)) then begin
                  plateandfiberidfibers = setIntersection(platefibers,fiberidfibers)
                  ifiber = plateandfiberidfibers[0]
                  uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
                endif else begin
                  uumessage = "Getting new plate" & print, uumessage
                  uuplotspec_init, plate, topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
                endelse
              endif else begin
                uumessage = "Getting new plate" & print, uumessage
                uuplotspec_init, plate, topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
              endelse
            endif else if (mjd ne '' and fiberid ne '') then begin
              uuplotspec_init, plate, topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
              platefibers = where(platelist eq plate[0], ncomplement=nonplatefibers)
              if (nonplatefibers ne n_elements(fiberidlist)) then begin
                mjdfibers = where(mjdlist eq mjd[0], ncomplement=nonmjdfibers)
                fiberidfibers = where(fiberidlist eq fiberid[0], ncomplement=nonfiberidfibers)
                if (nonmjdfibers ne n_elements(fiberidlist)) and (nonfiberidfibers ne n_elements(fiberidlist)) then begin
                  plateandmjdandfiberidfibers = setIntersection(setIntersection(platefibers,mjdfibers),fiberidfibers)
                  ifiber = plateandmjdandfiberidfibers[0]
                  uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
                endif else begin
                  uumessage = "Getting new plate/mjd" & print, uumessage
                  uuplotspec_init, plate, mjd=mjd, topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
                endelse
              endif else begin
                uumessage = "Getting new plate/mjd" & print, uumessage
                uuplotspec_init, plate, mjd=mjd, topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
              endelse
            endif
          endif
          widget_control, uuState.plateid, set_value=strtrim(platelist[ifiber],2)
          widget_control, uuState.mjdid, set_value=strtrim(mjdlist[ifiber],2)
          widget_control, uuState.fiberid, set_value=strtrim(fiberidlist[ifiber],2)
          nfibertext = strtrim(ifiber+1,2)+' of '+strtrim(n_elements(fiberidlist),2)+'              '
          widget_control, uuState.nfiberid, set_value = strmid(nfibertext,0,16)
        end
        else: begin
          uuplotspec1, platelist[ifiber], fiberidlist[ifiber], mjd=mjdlist[ifiber], topdir=topdirlist[ifiber], run1d=run1dlist[ifiber], run2d=run2dlist[ifiber]
        end
      endcase
      
    end
    
    'uuclasslist': begin
      widget_control, uuState.classid, get_value=class
      if (event.index gt 0) then uuState.class = class[event.index] else uuState.class=''
    end
    
    'uuzconflist': begin
      widget_control, uuState.zconfid, get_value=zconf
      if (event.index gt 0) then uuState.zconf = zconf[event.index] else uuState.zconf=''
    end
    
    'uuissuelist': begin
      uuState.issue = uuState.issues[event.index]
    end
    'uurecentcommentlist': begin
      uuDatabase_select_comment, strtrim(recentcommentlist[event.index].commentid,2)
    end
    
    'uuyannybutton': begin
      uuDatabase_make_yanny
    end
    
    'uuyannygroup': begin
      widget_control, uuState.yannyid, get_value=yanny
      if (yanny eq '') then begin
        yanny = uuState.yanny
        widget_control, uuState.yannyid, set_value=yanny
      endif
      if (strlen(yanny) gt 128) then begin
        yanny = strmid(yanny,0,128)
        widget_control, uuState.yannyid, set_value=yanny
      endif
      uuState.yanny = yanny
      widget_control, uuState.yannygroupid, get_value=yannygroup
      if (yanny eq 'spinspect' and yannygroup ne 0) then begin
        widget_control, uuState.yannygroupid, set_value=0
        uumessage = 'Choose a custom name for the Yanny file to group all fibers together'
      endif
      uuState.yannygroup = yannygroup
    end
    
    else: begin
      uuState.valid = 1
      widget_control, uuState.zid, get_value=z
      widget_control, uuState.commentid, get_value=comment
      widget_control, uuState.yannyid, get_value=yanny
      if (yanny eq '') then begin
        yanny = uuState.yanny
        widget_control, uuState.yannyid, set_value=yanny
      endif
      if (strlen(yanny) gt 128) then begin
        yanny = strmid(yanny,0,128)
        widget_control, uuState.yannyid, set_value=yanny
      endif
      uuState.yanny = yanny
      if (z eq '') and (comment eq '') and (uuState.class eq '') then begin
        uumessage = 'Enter a comment, or manual z or class inspection'
        uuState.valid=0
      endif else begin
        if (strlen(comment) le 256) then begin
          uuState.comment = comment
        endif else begin
          uumessage = 'Comment cannot exceed 256 characters'
          uuState.valid = 0
        endelse
        if (is_numeric(z)) then begin
          uuState.z = z
        endif else begin
          if (uumessage ne '') then uumessage = uumessage + ', and z must be numeric' else uumessage = 'z must be numeric'
          uuState.valid = 0
        endelse
      endelse
      
      if (uuState.valid) then begin
        uuDatabase_comment
        uuDatabase_recentcommentlist
        widget_control, uuState.recentcommentid, set_value=recentcommentlist.comment
      endif
    end
    
  endcase
  
  widget_control, uuState.messageid, set_value=uumessage[0]
  return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuPlotspecBase_refresh
  ;==============================================================================
  ; Plotspec Base Function: refreshes the uuplotspec floating window to keep
  ; it up to date during navigation functions and database functions.
  ;==============================================================================
  common splot_state, state, graphkeys
  common plotspec_state, platelist, fiberidlist, mjdlist, topdirlist, run1dlist, run2dlist, ifiber, keyword, keywordset, uumessage
  common uuPlotspecBase_state, uuState, recentcommentlist
  if (uuState.action eq 'login') or (uuState.action eq 'logout') then begin
    widget_control, uuState.commentheader, /DESTROY
    widget_control, uuState.commentbase, /DESTROY
    widget_control, uuState.yannybase, /DESTROY
  endif
  widget_control, uuState.plateid, set_value = strtrim(platelist[ifiber],2)
  widget_control, uuState.mjdid, set_value = strtrim(mjdlist[ifiber],2)
  widget_control, uuState.fiberid, set_value=strtrim(fiberidlist[ifiber],2)
  widget_control, uuState.ifiberid, set_value = strtrim(ifiber+1,2)
  nfibertext = 'of '+strtrim(n_elements(fiberidlist),2)+'    '
  widget_control, uuState.nfiberid, set_value = strmid(nfibertext,0,16)
  if (keywordset.zmanual) then begin
    if (keyword.zmanual[0] ne keyword.zmanual[1]) then zmanual = strtrim(keyword.zmanual,2) else zmanual = [strtrim(keyword.zmanual[0],2),'']
  endif else zmanual = ['','']
  widget_control, uuState.zmanual0id, set_value = zmanual[0]
  widget_control, uuState.zmanual1id, set_value = zmanual[1]
  if (keywordset.znum) then znum = strtrim(keyword.znum,2) else znum = ''
  widget_control, uuState.znumid, set_value = znum
  if (keyword.nsmooth eq 1) then keywordset.nsmooth = 0
  if (keywordset.nsmooth) then nsmooth = strtrim(keyword.nsmooth,2) else nsmooth = '1'
  widget_control, uuState.nsmoothid, set_value = nsmooth
  ;==============================================================================
  ; LOGGED IN
  ;==============================================================================
  if (uuState.loggedin) then begin
    widget_control, uuState.loginbuttonid, set_value = 'Logout'
    uuState.commentheader = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'commentheader')
    uuState.commentbase = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'commentbase')
    uuState.yannybase = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'yannybase')
    commentheader_col0 = widget_base(uuState.commentheader,/align_right,  group_leader = state.base_id,  /col, uvalue = 'commentheader_col0')
    commentheader_col1 = widget_base(uuState.commentheader,/align_right,  group_leader = state.base_id,  /col, uvalue = 'commentheader_col1');
    commentheader_col2 = widget_base(uuState.commentheader,/align_right,  group_leader = state.base_id,  /col, uvalue = 'commentheader_col2');
    commentheader_col3 = widget_base(uuState.commentheader,/align_right,  group_leader = state.base_id,  /col, uvalue = 'commentheader_col3');
    commentheader_col4 = widget_base(uuState.commentheader,/align_right,  group_leader = state.base_id,  /col, uvalue = 'commentheader_col4');
    commentbase_col0 = widget_base(uuState.commentbase,/align_right,  group_leader = state.base_id,  /col, uvalue = 'commentbase_col0')
    commentbase_col1 = widget_base(uuState.commentbase,/align_right,  group_leader = state.base_id,  /col, uvalue = 'commentbase_col1');
    commentbase_col2 = widget_base(uuState.commentbase,/align_right,  group_leader = state.base_id,  /col, uvalue = 'commentbase_col2');
    yannybase_col0 = widget_base(uuState.yannybase,/align_right,  group_leader = state.base_id,  /col, uvalue = 'yannybase_col0')
    yannybase_col1 = widget_base(uuState.yannybase,/align_right,  group_leader = state.base_id,  /col, uvalue = 'yannybase_col1');
    yannybase_col2 = widget_base(uuState.yannybase,/align_right,  group_leader = state.base_id,  /col, uvalue = 'yannybase_col2');
    yannybase_col3 = widget_base(uuState.yannybase,/align_right,  group_leader = state.base_id,  /col, uvalue = 'yannybase_col3');
    widget_control, uuState.messageid, set_value = ' '
    void = widget_label(commentheader_col0,value ='Provide Feedback via Manual Inspection           ')
    uuState.recentcommentid = widget_droplist(commentbase_col0,title='        ', uvalue='uurecentcommentlist',value = recentcommentlist.comment)
    uuState.issueid = widget_droplist(commentbase_col0,title='  Issue:', uvalue='uuissuelist',value = uuState.issues)
    uuState.commentid = cw_field(commentbase_col0,/row,title = 'Comment:',uvalue='uucommentfield',value = '', /string,/return_events,xsize=50)
    uuState.yannyid = cw_field(yannybase_col0,/row,title = 'Yanny:',uvalue='uuyannyfield',value = 'spinspect', /string,/return_events,xsize=10)
    yannygroup = ['Group by plate-mjd','All fibers']
    uuState.yannygroupid = cw_bgroup(yannybase_col1,yannygroup,/row,/exclusive,set_value=0,uvalue='uuyannygroup')
    void = cw_field(yannybase_col2,/row,title = '   Dir:',uvalue='uuusername',value = uuState.username, /string,/noedit, xsize=10)
    void = widget_label(commentheader_col1,value = 'RUN1D:')
    uustate.run1d = widget_label(commentheader_col2,uvalue='uurun1d',value = run1dlist[ifiber], xsize=60, frame=1)
    void = widget_label(commentheader_col3,value = 'RUN2D:')
    uustate.run2d = widget_label(commentheader_col4,uvalue='uurun2d',value = run2dlist[ifiber], xsize=60, frame=1)
    void = widget_label(uuState.commentheader,value = ' ')
    uuState.zid = cw_field(commentbase_col1,/row,title = '     z:',uvalue='uuzfield',value = '', /string,/return_events,xsize=9)
    uuState.classid = widget_droplist(commentbase_col1,title=' Class:', uvalue='uuclasslist',value = [' ','Galaxy','QSO','Star'])
    uuState.zconfid = widget_droplist(commentbase_col1,title='z conf:', uvalue='uuzconflist',value = ['      ','3','2','1','0'])
    submit = widget_button(commentbase_col2,value = 'Submit'+string(10B)+'Feedback',uvalue='uusubmitbutton', xsize=80,ysize=110)
    yanny = widget_button(yannybase_col3,value = 'Make Yanny',uvalue='uuyannybutton', xsize=80,ysize=36)
    void = widget_label(uuState.commentbase,value = ' ')
    void = widget_label(uuState.yannybase,value = ' ')
    
  ;==============================================================================
  ; NOT LOGGED IN
  ;==============================================================================
  endif else begin
    widget_control, uuState.loginbuttonid, set_value = 'Login'
    widget_control, uuState.messageid, set_value = 'Please Login to provide spectrum feedback to http://boss.astro.utah.edu'
    uuState.commentheader = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'commentheader')
    uuState.commentbase = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'commentbase')
    uuState.yannybase = widget_base(uuState.uuplotspecbase,/align_right,  group_leader = state.base_id,  /row,  uvalue = 'yannybase')
  endelse
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuDatabase_member, action
  ;==============================================================================
  ; Database Function: authenticate member (login) and get recent comments
  ;==============================================================================
  common uuPlotspecBase_state, uuState, recentcommentlist
  post = {func:'member',mysql:'boss',ver:'040620111405',siteID:'1',action:action,username:uuState.username,password:uuState.password,remote:uuState.remote}
  item = {loggedin:'',username:'',sid:''}
  uuDatabase_select, post, uuState.sid, item, select
  if (n_elements(select) eq 1) then begin
    uuState.loggedin=fix(select[0].loggedin)
    uuState.sid=select[0].sid
    uuState.username=select[0].username
  endif else begin
    uuState.loggedin=0
    uuState.sid=''
    uuState.username=''
  endelse
  if (uuState.loggedin) then uuDatabase_recentcommentlist
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuDatabase_comment
  ;==============================================================================
  ; Database Function: insert comment and receive message response from database
  ;==============================================================================
  common plotspec_state, platelist, fiberidlist, mjdlist, topdirlist, run1dlist, run2dlist, ifiber, keyword, keywordset, uumessage
  common uuPlotspecBase_state, uuState, recentcommentlist
  
  escapedcomment = uuState.comment
  if (strpos(escapedcomment,'%') ge 0) then escapedcomment = STRJOIN(STRSPLIT(escapedcomment,'%',/EXTRACT),'%25')
  if (strpos(escapedcomment,'&') ge 0) then escapedcomment = STRJOIN(STRSPLIT(escapedcomment,'&',/EXTRACT),'%26')
  post = {func:'comment',mysql:'boss',ver:'040620111405',action:'insert',username:uuState.username,password:uuState.password,plate:platelist[ifiber],mjd:mjdlist[ifiber],fiberid:fiberidlist[ifiber],run1d:run1dlist[ifiber],run2d:run2dlist[ifiber],comment:escapedcomment,z:uuState.z,issue:uuState.issue,class:uuState.class,zconf:uuState.zconf,yanny:uuState.yanny,yannygroup:uuState.yannygroup}
  item = {message:'',confirm:'',loggedin:'',sid:''}
  uuDatabase_select, post, uuState.sid, item, select
  if (n_elements(select) eq 1) then begin
    uuState.loggedin=fix(select[0].loggedin)
    uuState.sid=select[0].sid
    if (select[0].confirm eq '1') then begin
      widget_control, uuState.zid, set_value=''
      widget_control, uuState.commentid, set_value=''
      widget_control, uuState.classid, set_droplist_select=0
      widget_control, uuState.zconfid, set_droplist_select=0
    endif
    uumessage=select[0].message
  endif else begin
    uuState.loggedin=0
    uuState.sid=''
  endelse
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function uuDatabase_download, url
  ;==============================================================================
  ; Database Function: download a dynamic resource from the database using a
  ; Java method to launch the OS default browser.
  ;==============================================================================
  catch = 1
  catch, awt_error
  if (awt_error ne 0) then begin
    catch,/cancel
    return, 1
  endif
  oJavaDesktop = OBJ_NEW('IDLJavaObject$Static$JAVA_AWT_DESKTOP', 'java.awt.Desktop')
  oJavaURI = OBJ_NEW('IDLJavaObject$JAVA_NET_URI', 'java.net.URI', url)
  if (oJavaDesktop->isDesktopSupported()) then begin
    oBrowser = oJavaDesktop->getDesktop();
    oBrowser->browse,oJavaURI
    OBJ_DESTROY, oBrowser
    catch = 0
  endif
  OBJ_DESTROY, oJavaURI, oJavaDesktop
  return, catch
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuDatabase_make_yanny
  ;==============================================================================
  ; Database Function: generate a yanny file from feedback in the database
  ; and email the files to the member, or do a local save to dir (needs widget)
  ;==============================================================================
  common uuPlotspecBase_state, uuState, recentcommentlist
  common plotspec_state, platelist, fiberidlist, mjdlist, topdirlist, run1dlist, run2dlist, ifiber, keyword, keywordset, uumessage
  
  post = {func:'comment',mysql:'boss',ver:'040620111405',action:'yanny',username:uuState.username,password:uuState.password}
  item = {count:'',inserted:'',updated:'',unchanged:'',url:''}
  uuDatabase_select, post, uuState.sid, item, select
  if (n_elements(select) eq 1) then begin
    count = select[0].count
    if (count gt '0') then begin
      inserted = select[0].inserted + " new"
      inserted = select[0].inserted & if (inserted eq '1') then inserted = "Created 1 new yanny file" else inserted = "Created " + inserted + " new yanny files"
      updated = select[0].updated & if (updated ne '0') then updated = ", updated "+updated else updated = ''
      unchanged = select[0].unchanged & if (unchanged ne '0') then unchanged = ", "+ unchanged + " are up-to-date" else unchanged = ''
      url = select[0].url
      uumessage = inserted+updated+unchanged+", for a total of "+count+"."
      catch = 1
      if (url ne '') then catch = uuDatabase_download(url)
      if (catch eq 1) then uuMessage = "Please login to http://boss.astro.utah.edu/spinspect/data using the sdss3 username/password."
    endif else if (count eq 0) then begin
      uumessage = "No yanny files were generated.  Please provide feedback first."
    endif else uumessage = "No yanny files were generated.  Please contact administrator."
  endif else begin
    uumessage = "No yanny files were generated."
  endelse
  print, uumessage
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuDatabase_recentcommentlist
  ;==============================================================================
  ; Database Function: select list of recent comments from database
  ;==============================================================================
  common uuPlotspecBase_state, uuState, recentcommentlist
    
  post = {func:'comment',mysql:'boss',ver:'040620111405',action:'recent',username:uuState.username,password:uuState.password}
  item = {comment:'',commentid:0L}
  uuDatabase_select, post, uuState.sid, item, select
  if (n_elements(select) eq 0) then begin
    recentcommentlist = [{comment:'Paste from recent comments                     ',commentid:0L}]
  endif else begin
    recentcommentlist = select
    recentcommentlist[0].comment='Paste from recent comments                     '
    recentcommentlist[0].commentid=0L
  endelse
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuDatabase_select_comment, commentid
  ;==============================================================================
  ; Database Function: select comment from database
  ;==============================================================================
  common uuPlotspecBase_state, uuState, recentcommentlist
  
  if (commentid ne 0) then begin
    post = {func:'comment',mysql:'boss',ver:'040620111405',action:'select',commentid:commentid}
    item = {comment:'',commentid:0L,issueid:0L}
    uuDatabase_select, post, uuState.sid, item, select
    if (n_elements(select) eq 1) then begin
      if (select[0].commentid ne '0') then begin
        widget_control, uuState.commentid, set_value=select[0].comment
        widget_control, uuState.issueid, set_droplist_select=long(select[0].issueid-1)
      endif else commentid = 0
    endif else commentid = 0
  endif
  if (commentid eq 0) then widget_control, uuState.commentid, set_value=''
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuDatabase_post, cmd, post, sid, response=response
  ;==============================================================================
  ; Database Function: communicate with database URL by HTTP POST
  ;==============================================================================
  url = ['http://boss.astro.utah.edu/','http://cosmo.astro.utah.edu/boss/']
  ;url = ['http://cosmo.astro.utah.edu/boss/','http://cosmo.astro.utah.edu/boss/']
  if (sid ne '') then cmdurl = url + cmd + '.php?PHPSESSID='+sid else cmdurl = url + cmd + '.php'
  response = uuDatabase_webget(cmdurl,POST=post)
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuDatabase_select, post, sid, item, select
  ;==============================================================================
  ; Database Function: parse response from database on selected item
  ;==============================================================================
  cmd = 'mysql'
  uuDatabase_post, cmd, post, sid, response=response
  if (response.Text ne 'NULL') and (response.Text ne '') then begin
    escape = '`'
    responses = STRSPLIT(response.Text, ESCAPE=escape, /EXTRACT,';')
    nresponses = N_ELEMENTS(responses)
    select = replicate(item,nresponses)
    if (nresponses gt 0) then begin
      select = replicate(item,nresponses)
      for i = 0,nresponses-1 do begin
        pairs = STRSPLIT(responses[i], ESCAPE=escape, /EXTRACT,',')
        npairs = N_ELEMENTS(pairs)
        itemtag = TAG_NAMES(item)
        if (npairs eq N_ELEMENTS(itemtag)) then begin
          for j=0,npairs-1 do begin
            pair = STRSPLIT(pairs[j], ESCAPE=escape, /EXTRACT,'=')
            if (N_ELEMENTS(pair) eq 2) then begin
              key = pair[0]
              value = pair[1]
              if (key = itemtag[j]) then select[i].(j)=value
            endif
          endfor
        endif else begin
          uumessage = "bad response from database, please contact administrator." & print, uumessage 
        endelse
      endfor
    endif else begin
      uumessage = "null response from database, please contact administrator." & print, uumessage
    endelse
  endif else begin
  endelse
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uuplotspec, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, znum=znum, nsmooth=nsmooth, $
    zline=zline, nosyn=nosyn, noerr=noerr, sky=sky, $
    ormask=ormask, andmask=andmask, $
    psfile=psfile, xrange=xrange, yrange=yrange, allexp=allexp, $
    restframe=restframe, zwarning=zwarning, $
    _EXTRA=Extra
  ;==============================================================================
  ; Main program: reorganize keywords for plotspec (uuplotspec_init) and
  ; setup the uuPlotspecBase floating window.
  ;==============================================================================
    
  common splot_state, state, graphkeys
  common plotspec_state, platelist, fiberidlist, mjdlist, topdirlist, run1dlist, run2dlist, ifiber, keyword, keywordset, uumessage
  
  if (getenv('IDLUTILS_DIR') eq '') then begin
    print, 'Please set your IDLUTILS_DIR environment variable, and start again'
    return
  endif
  
  if (n_params() LT 1) then begin
    doc_library, 'uuplotspec'
    return
  endif
   
  if (xregistered ('splot')) then widget_control, state.base_id, /destroy
  keyword = {zmanual:[0D,0D],znum:0L,nsmooth:1L,manualz:0D,manualclass:''}
  keywordset = {zmanual:0L,znum:keyword_set(znum),nsmooth:keyword_set(nsmooth),zline:keyword_set(zline),nosyn:keyword_set(nosyn),noerr:keyword_set(noerr),sky:keyword_set(sky),ormask:keyword_set(ormask),andmask:keyword_set(andmask),psfile:keyword_set(psfile),xrange:keyword_set(xrange),yrange:keyword_set(yrange),allexp:keyword_set(allexp),restframe:keyword_set(restframe),zwarning:keyword_set(zwarning)}
  if (not keyword_set(topdir)) or (n_elements(topdir) lt n_elements(plate)) then topdir = replicate(getenv('BOSS_SPECTRO_REDUX'),n_elements(plate))
  if (not keyword_set(run1d)) or (n_elements(run1d) lt n_elements(plate)) then run1d = replicate(getenv('RUN1D'),n_elements(plate)) 
  if (not keyword_set(run2d)) or (n_elements(run2d) lt n_elements(plate)) then run2d = replicate(getenv('RUN2D'),n_elements(plate)) 
  if (keyword_set(zmanual)) then begin
    if (n_elements(zmanual) eq 1) then begin
      keywordset.zmanual = 1
      keyword.zmanual = [zmanual[0],zmanual[0]]
    endif else if (n_elements(zmanual) gt 1) then begin
      keywordset.zmanual = 1
      keyword.zmanual = zmanual[0:1]
    endif
  endif
  if (keywordset.znum) then keyword.znum = znum
  if (keywordset.nsmooth) then keyword.nsmooth = nsmooth
  uumessage =''
  uuplotspec_init, plate, fiberid, mjd=mjd, topdir=topdir, run1d=run1d, run2d=run2d, $
    psfile=psfile, xrange=xrange, yrange=yrange, $
    _EXTRA=Extra
    
  if (uumessage eq '') then begin
    uuPlotspecBase
    splot_set_minmax
  endif
  
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
