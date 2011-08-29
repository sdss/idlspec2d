;+
; NAME:
;   plotspec
;
; PURPOSE:
;   Routine for plotting spectra from Princeton-1D spectro outputs.
;
; CALLING SEQUENCE:
;   plotspec, plate, [ fiberid, mjd=, znum=, zmanual=, $
;    nsmooth=, /zline, /nosyn, /noerr, $
;    /sky, /ormask, /andmask, psfile=, /restframe, /netimage, $
;    /zwarning, /allexp, _EXTRA= ]
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
;   netimage   - If set, then launch a Netscape browser with the object
;                image from Steve Kent's web site.  This only works if
;                Netscape is running and has permissions at the site
;                "http://sdssmosaic.fnal.gov:8015". ???
;                This is disabled if PSFILE is set.
;   zwarning   - If set, then only select those non-sky fibers where the
;                ZWARNING flag has been set; can be used with or without
;                specifying fiber numbers with FIBERID.
;   allexp     - If set, then plot all the individual exposure spectra,
;                rather than the co-added spectrum.
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
;     IDL> plotspec, 401, 100
;
;   The spectrum is shown in white, the errors in red (except masked points
;   are set to zero), and the best-fit eigenspectrum in blue. The mouse
;   buttons will zoom in (left), recenter (center), or zoom out (right).
;   The frame can be saved as a PostScript file by selecting File->WriteEPS
;   from the left-hand corner. 
;
;   Make the same plot, but boxcar-smooth the spectrum and limit the
;   wavelength range to [4000,5000] Angstroms:
;     IDL> plotspec, 401, 100, nsmooth=10, xrange=[5000,6000]
;
;   Some plates are observed on multiple nights. To select one of the two
;   observations of plate 306: 
;     IDL> plotspec, 306, 20, mjd=51690
;
;   Loop through all the spectra for plate 401, interactively:
;     IDL> plotspec, 401
;
;   Plot all the spectra from plate 401 to a single PostScript file:
;     IDL> plotspec, 401, /psfile
;
;   Plot all the spectra from plate 401 to 640 individual PostScript files:
;     IDL> plotspec, 401, lindgen(640)+1, /psfile
;
;   Plot a list of 3 objects, each with its own plate, MJD, and fiberid:
;     IDL> plate = [400,400,401]
;     IDL> mjd = [51820,51820,51788]
;     IDL> fiberid = [10,11,20]
;     IDL> plotspec, plate, mjd=mjd, fiberid
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
; INTERNAL SUPPORT ROUTINES:
;   plotspec_mask
;   plotspec1
;
; REVISION HISTORY:
;   01-Sep-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro plotspec_mask, wave, thismask, psfile=psfile, nolabel=nolabel, $
 _EXTRA=Extra

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
;------------------------------------------------------------------------------
pro plotspec1, plate, fiberid, mjd=mjd, znum=znum, zmanual=zmanual, $
 nsmooth=nsmooth1, zline=q_zline, nosyn=nosyn, noerr=noerr, sky=sky, $
 ormask=ormask, andmask=andmask, $
 psfile=psfile, xrange=passxr, yrange=passyr, noerase=noerase, $
 restframe=restframe, netimage=netimage, allexp=allexp, _EXTRA=Extra

   cspeed = 2.99792458e5
   textcolor = 'green'
   linecolor = 'magenta'
   orcolor = 'yellow'
   andcolor = 'red'

   readspec, plate, fiberid, mjd=mjd, znum=znum, flux=objflux, $
    wave=wave, plug=plug, zans=zans, _EXTRA=Extra, /silent
   if (NOT keyword_set(objflux)) then begin
      print, plate, mjd, fiberid, $
       format='("Spectrum not found for plate=", i4, " MJD=", i5, " fiber=", i3)'
      return
   endif
   if (keyword_set(allexp)) then begin
      readonespec, plate, fiberid, mjd=mjd, wave=allwave, flux=allflux, $
                   _EXTRA=Extra, /silent
      ndim = size(allflux,/n_dimen)
      if (ndim EQ 1) then nexp = 1 $
       else nexp = (size(allflux,/dimens))[1]
   endif
   if (keyword_set(restframe)) then begin
      wave = wave / (1. + zans.z)
      if (keyword_set(allwave)) then allwave = allwave / (1. + zans.z)
   endif
   if (NOT keyword_set(noerr)) then $
    readspec, plate, fiberid, mjd=mjd, flerr=objerr, _EXTRA=Extra, /silent
   if (keyword_set(zmanual)) then begin
      readspec, plate, fiberid, mjd=mjd, invvar=objivar, loglam=loglam, $
       objhdr=hdr, _EXTRA=Extra, /silent
      if (zmanual[0] LT 1.) then eigenfile = 'spEigenGal-*.fits' $
       else eigenfile = 'spEigenQSO-*.fits'
      npoly = 0
      if (n_elements(zmanual) EQ 1) then zrange=[zmanual,zmanual] $
       else zrange=zmanual[0:1]
      res_manual = zfind(objflux, objivar, hdr=hdr, $
       eigenfile=eigenfile, npoly=npoly, zmin=zrange[0], zmax=zrange[1], $
       pspace=1, nfind=1, width=1)
      synflux = synthspec(res_manual, loglam=loglam)
   endif else begin
      if (NOT keyword_set(nosyn)) then $
       readspec, plate, fiberid, mjd=mjd, znum=znum, synflux=synflux, $
        _EXTRA=Extra, /silent
   endelse
   if (keyword_set(sky)) then $
    readspec, plate, fiberid, mjd=mjd, sky=sky, _EXTRA=Extra, /silent
   if (keyword_set(ormask)) then $
    readspec, plate, fiberid, mjd=mjd, ormask=ormask, _EXTRA=Extra, /silent
   if (keyword_set(andmask)) then $
    readspec, plate, fiberid, mjd=mjd, andmask=andmask, _EXTRA=Extra, /silent
   if (keyword_set(zans) AND keyword_set(q_zline)) then $
    readspec, plate, fiberid, mjd=mjd, zline=zline, lineflux=lineflux, $
     _EXTRA=Extra, /silent

   if (keyword_set(nsmooth1)) then nsmooth = nsmooth1 $
    else nsmooth = 1

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
   if (keyword_set(ormask) OR keyword_set(andmask)) then $
    xrange[1] = 1.15 * xrange[1] - 0.15 * xrange[0]

   title = 'Plate ' + strtrim(string(plate),2) $
    + '  Fiber ' + strtrim(string(fiberid),2) $
    + '  MJD=' + strtrim(string(mjd),2)
   if (keyword_set(restframe)) then xtitle = 'Rest-Frame Wavelength [Ang]' $
    else xtitle = 'Observed Wavelength [Ang]'
   if (keyword_set(psfile)) then begin
      djs_plot, xrange, yrange, /nodata, xrange=xrange, yrange=yrange, $
       xtitle=xtitle, ytitle=TeXtoIDL('Flux [10^{-17} erg/s/cm^2/Ang]'), $
       title=title, charsize=csize, _EXTRA=KeywordsForSplot, /xstyle, /ystyle
      if (keyword_set(allexp)) then begin
         for iexp=0, nexp-1 do $
          djs_oplot, allwave[*,iexp], allflux[*,iexp], $
           psym=(nsmooth GT 1) ? 0 : 3, _EXTRA=KeywordsForSplot
      endif else begin
         djs_oplot, wave, objflux, _EXTRA=KeywordsForSplot
      endelse
      if (NOT keyword_set(noerr)) then $
       djs_oplot, wave, objerr, color='red', _EXTRA=KeywordsForSplot
      if (keyword_set(sky)) then $
       djs_oplot, wave, sky, color='green', lw=2, _EXTRA=KeywordsForSplot
      if (keyword_set(synflux)) then $
       djs_oplot, wave, synflux, color='blue', lw=2, _EXTRA=KeywordsForSplot
   endif else begin
      if (NOT keyword_set(noerase)) then $
       splot, xrange, yrange, /nodata, xrange=xrange, yrange=yrange, $
        xtitle=xtitle, ytitle=TeXtoIDL('Flux [10^{-17} erg/s/cm^2/Ang]'), $
        title=title, charsize=csize, _EXTRA=KeywordsForSplot
      if (keyword_set(allexp)) then begin
         for iexp=0, nexp-1 do $
          soplot, allwave[*,iexp], allflux[*,iexp], $
           _EXTRA=KeywordsForSplot
      endif else begin
         soplot, wave, objflux, _EXTRA=KeywordsForSplot
      endelse
      if (NOT keyword_set(noerr)) then $
       soplot, wave, objerr, color='red', _EXTRA=KeywordsForSplot
      if (keyword_set(sky)) then $
       soplot, wave, sky, color='green', lw=2, _EXTRA=KeywordsForSplot
      if (keyword_set(synflux)) then $
       soplot, wave, synflux, color='blue', lw=2, _EXTRA=KeywordsForSplot
   endelse

   xpos = 0.9 * !x.window[0] + 0.1 * !x.window[1]
   dypos = 0.05 * (!y.window[0] - !y.window[1])
   ypos = !y.window[1] + 1.5 * dypos

   if (keyword_set(zmanual) OR keyword_set(zans)) then begin
      if (keyword_set(zmanual)) then begin
         zstring = res_manual.tfile $
          + '  z=' + string(res_manual.z,format='(f8.5)')
      endif else begin
         zstring = zans.class + ' ' + zans.subclass
         cz = zans.z * cspeed
         if (abs(cz) LT 3000) then $
           zstring += '  cz=' + string(cz,format='(f6.0)') + ' km/s' $
          else $
           zstring += '  z=' + string(zans.z,format='(f8.5)')
         if (zans.zwarning NE 0) then $
          zstring +=  ' (' $
           + sdss_flagname('ZWARNING', zans.zwarning, /concat) + ')'
         if (keyword_set(znum)) then $
          zstring += ' (fit #' + strtrim(string(znum),2) + ')'
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
       * (1 + zline.linez * (keyword_set(restframe) EQ 0))
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

   if (keyword_set(ormask)) then begin
      plotspec_mask, wave, ormask, psfile=psfile, $
       psym=1, symsize=0.6, color=orcolor, $
       nolabel=keyword_set(andmask), _EXTRA=KeywordsForSplot
   endif

   if (keyword_set(andmask)) then begin
      plotspec_mask, wave, andmask, psfile=psfile, $
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

   if (keyword_set(netimage) AND NOT keyword_set(psfile)) then begin
      netstring = 'http://sdssmosaic.fnal.gov:8015/template/tsSingle.tml?run=' $
       + strtrim(string(plug.objid[0]),2) $
       + '&camcol=' + strtrim(string(plug.objid[2]),2) $
       + '&field=' + strtrim(string(plug.objid[3]),2) $
       + '&ra=' + strtrim(string(plug.ra),2) $
       + '&dec=' + strtrim(string(plug.dec),2)
      spawn, '\netscape -remote "openURL(' + netstring + ')"'
print,netstring
   endif

   return
end
;------------------------------------------------------------------------------
pro plotspec, plate, fiberid, mjd=mjd, znum=znum, nsmooth=nsmooth, $
 zline=zline, nosyn=nosyn, noerr=noerr, sky=sky, $
 ormask=ormask, andmask=andmask, $
 psfile=psfile, xrange=xrange, yrange=yrange, noerase=noerase, $
 restframe=restframe, netimage=netimage, zwarning=zwarning, allspec=allspec, $
 _EXTRA=Extra

   if (n_params() LT 1) then begin
      doc_library, 'plotspec'
      return
   endif

   quiet = !quiet
   !quiet = 1

   ;----------
   ; If MJD is not set, then find the MJD for each plate

   nplate = n_elements(plate)
   if (NOT keyword_set(mjd)) then begin
      mjd = lonarr(nplate)
      for iplate=0, nplate-1 do begin
         mjd1 = 0
         readspec, plate[iplate], mjd=mjd1, _EXTRA=Extra, /silent
         if (NOT keyword_set(mjd1)) then begin
            print, 'No MJD found for plate ', plate[iplate]
            !quiet = quiet
            return
         endif
         mjd[iplate] = mjd1
      endfor
   endif else begin
      if (n_elements(mjd) NE nplate) then begin
         print, 'Number of elements in PLATE and MJD do not agree'
         !quiet = quiet
         return
      endif
   endelse
   
   ;----------
   ; If /ZWARNING is set, then find the flagged fibers to plot.

   if (keyword_set(zwarning)) then begin
      if (keyword_set(fiberid)) then begin
         print, 'FIBERID and /ZWARNING cannot both be set.'
         !quiet = quiet
         return
      endif

      for iplate=0L, nplate-1L do begin
         readspec, plate[iplate], mjd=mjd[iplate], $
          _EXTRA=Extra, /silent, zans=zans
         if (NOT keyword_set(zans)) then begin
            print, 'No spZ file found for selecting ZWARNING flags'
            !quiet = quiet
            return
         endif
         indx = where((zans.zwarning AND 1) EQ 0 AND zans.zwarning NE 0, nthis)
         if (nthis GT 0) then begin
            if (NOT keyword_set(fiberid)) then begin
               platelist = replicate(plate[iplate], nthis)
               mjdlist = replicate(mjd[iplate], nthis)
               fiberid = zans[indx].fiberid
            endif else begin
               platelist = [platelist, replicate(plate[iplate], nthis)]
               mjdlist = [mjdlist, replicate(mjd[iplate], nthis)]
               fiberid = [fiberid, zans[indx].fiberid]
            endelse
         endif
      endfor
      if (NOT keyword_set(fiberid)) then begin
         print, 'No non-sky fibers with ZWARNING flag set'
         !quiet = quiet
         return
      endif
      nfiber = n_elements(fiberid)
      print, 'Selecting ', nfiber, ' non-sky fibers with ZWARNING flag set'
   endif

   ;----------
   ; Set FIBERID to [1,...,NFIBER] (for each plate) if not set.
   ;
   ; If writing to a PostScript file, then all plots are in the same file
   ; either if PSFILE is that file name, or if FIBERID is not specified
   ; (and then all spectra are being plotted).

   if (NOT keyword_set(fiberid)) then begin
      readspec, plate, mjd=mjd, 0*plate+1, nfiber=nfiber_tmp, $
       _EXTRA=Extra, /silent
      nfiber_tot = long(total(nfiber_tmp))
      if (nfiber_tot EQ 0) then begin
         print, 'No fibers found'
         !quiet = quiet
         return
      endif
      platelist = lonarr(nfiber_tot)
      mjdlist = lonarr(nfiber_tot)
      fiberid = lonarr(nfiber_tot)
      j = 0L
      for iplate=0L, nplate-1L do begin
         platelist[j:j+nfiber_tmp[iplate]-1] = plate[iplate]
         mjdlist[j:j+nfiber_tmp[iplate]-1] = mjd[iplate]
         fiberid[j:j+nfiber_tmp[iplate]-1] = lindgen(nfiber_tmp[iplate]) + 1
         j += nfiber_tmp[iplate]
      endfor
      if (keyword_set(psfile)) then begin
         q_onefile = 1
         psfilename = string(plate[0], mjd[0], $
          format='("spec-",i4.4,"-",i5.5,".ps")')
      endif
   endif

;   if (min(fiberid) LT 1 OR max(fiberid) GT 640) then begin
;      print, 'Invalid FIBERID (must be between 1 and 640)'
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
      endif else begin
         platelist = plate
         mjdlist = mjd
      endelse
   endif else if (n_elements(platelist) NE n_elements(fiberid)) then begin
      print, 'Number of elements in PLATE and FIBERID do not agree.'
      !quiet = quiet
      return
   endif

   ;----------
   ; Loop over each plot

   ifiber = 0L
   while (ifiber LT nfiber) do begin

      ;----------
      ; Open the PostScript file if appropriate

      if (keyword_set(psfile)) then begin
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
      endif

      plotspec1, platelist[ifiber], fiberid[ifiber], mjd=mjdlist[ifiber], $
       znum=znum, zmanual=zmanual, nsmooth=nsmooth, zline=zline, $
       nosyn=nosyn, noerr=noerr, $
       sky=sky, ormask=ormask, andmask=andmask, psfile=psfile, $
       xrange=xrange, yrange=yrange, noerase=noerase, netimage=netimage, $
       restframe=restframe, allexp=allexp, _EXTRA=Extra

      if (keyword_set(psfile)) then begin
         if (NOT keyword_set(q_onefile) OR ifiber EQ nfiber-1) then dfpsclose
         ifiber = ifiber + 1
      endif else begin
         if (ifiber LT nfiber-1) then begin
            if (keyword_set(nsmooth)) then $ 
             sstring = ' (currently=' + strtrim(string(nsmooth),2) + ')' $
            else $
             sstring = ''

            print, 'Press b=back one fiber'
            print, '      p=select new plate'
            print, '      f=select new fiber number'
            print, '      n=change which PCA-fit to plot'
            print, '      q=quit (and enter interactive mode for this plot)'
            print, '      s=change boxcar smoothing' + sstring
            print, '      x=change X plotting range'
            print, '      y=change Y plotting range'
            print, '      z=manual z'
            print, '      v=view reconstructed frame'
            print, '      any other key=forward'

            cc = strupcase(get_kbrd(1))
            print, cc
            case cc of
            'V': begin
                    if (keyword_set(getenv('PHOTOOP_DIR'))) then begin
                       readspec, platelist[ifiber], fiberid[ifiber], $
                        mjd=mjdlist[ifiber], zans=zans, _EXTRA=Extra, /silent
                       plotspec_image, ra=zans.plug_ra, dec=zans.plug_dec, $
                        cutout=300, /calibrate, /register, /allid
                    endif else begin
                       print, 'Need to set up photoop product to display images'
                    endelse
                 end
            'B': begin
                    ifiber = (ifiber - 1) > 0
                    zmanual = 0.
                 end
            'P': begin
                    read, plate, mjd, prompt='Enter new plate and MJD (enter 0 for unknown MJD): '
                    if (NOT keyword_set(mjd)) then $
                     readspec, plate, mjd=mjd, _EXTRA=Extra, /silent
                    if (NOT keyword_set(mjd)) then begin
                       print, 'MJD not found for plate ', plate
                       !quiet = quiet
                       return
                    endif
                    readspec, plate, mjd=mjd, nfiber=nfiber, $
                     _EXTRA=Extra, /silent
                    platelist = replicate(plate,nfiber)
                    mjdlist = replicate(mjd,nfiber)
                    zmanual = 0.
                 end
            'N': begin
                    read, znum, prompt='Enter 1=best redshift, 2=2nd best, ...: '
                    znum = long(znum) > 0
                    zmanual = 0.
                 end
            'F': begin
                    read, newfiber, prompt='Enter new fiber number: '
                    ifiber = ((long(newfiber)-1) > 0) < (nfiber-1)
                    zmanual = 0.
                 end
            'Q': ifiber = nfiber
            'S': begin
                    read, nsmooth, prompt='Enter boxcar smoothing width (0=none): '
                    nsmooth = long(nsmooth) > 0
                 end
            'X': begin
                    read, xmin, xmax, prompt='Enter new X range values (0 0=full range): '
                    if (xmin EQ 0 AND xmax EQ 0) then xrange = 0 $
                     else xrange = [xmin, xmax]
                 end
            'Y': begin
                    read, ymin, ymax, prompt='Enter new Y range values (0 0=full range): '
                    if (ymin EQ 0 AND ymax EQ 0) then yrange = 0 $
                     else yrange = [ymin, ymax]
                 end
            'Z': begin
                    read, zmanual, prompt='Enter redshift guess (zmin zmax)...: '
                 end
            else: begin
                     ifiber = ifiber + 1
                     zmanual = 0.
                 end
            endcase
         endif else begin
            ifiber = nfiber
         endelse
      endelse
   endwhile

   !quiet = quiet
   return
end
;------------------------------------------------------------------------------
