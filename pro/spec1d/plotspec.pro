;+
; NAME:
;   plotspec
;
; PURPOSE:
;   Routine for plotting spectra from Princeton-1D spectro outputs.
;
; CALLING SEQUENCE:
;   plotspec, plate, [ fiberid, mjd=, znum=, nsmooth=, /zline, /nosyn, /noerr, $
;    /ormask, /andmask, psfile=, /netimage, /zwarning, topdir=, _EXTRA= ]
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
;   nsmooth    - If set, then boxcar smooth both the object and synthetic
;                spectra with a width equal to NSMOOTH.
;   zline      - If set, then overplot the emission line fits.
;   nosyn      - If set, then do not overplot the synthetic fit spectrum.
;   noerr      - If set, then do not overplot the error vector.
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
;   netimage   - If set, then launch a Netscape browser with the object
;                image from Steve Kent's web site.  This only works if
;                Netscape is running and has permissions at the site
;                "http://sdssmosaic.fnal.gov:8015".
;                This is disabled if PSFILE is set.
;   zwarning   - If set, then only select those non-sky fibers where the
;                ZWARNING flag has been set; can be used with or without
;                specifying fiber numbers with FIBERID.
;   topdir     - Top-level directory for data; default to the environment
;                variable $SPECTRO_DATA.
;   _EXTRA     - Kewords for SPLOT, such as XRANGE, YRANGE, THICK.
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
 _EXTRA=KeywordsForPlot

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
             _EXTRA=KeywordsForPlot
         endif else begin
            soplot, [wave[ipix]], [replicate(ypos,npix)], $
             _EXTRA=KeywordsForPlot
         endelse
      endif
      if (NOT keyword_set(nolabel)) then begin
         if (keyword_set(psfile)) then begin
            djs_xyouts, !x.window[1], ynorm, bitlabel[ilabel]+' ', $
             charsize=1.4, align=1.0, _EXTRA=KeywordsForPlot, /normal
         endif else begin
            sxyouts, !x.window[1], ynorm, bitlabel[ilabel]+' ', $
             charsize=1.4, align=1.0, _EXTRA=KeywordsForPlot, /normal
         endelse
      endif
   endfor

   return
end
;------------------------------------------------------------------------------
pro plotspec1, plate, fiberid, mjd=mjd, znum=znum, nsmooth=nsmooth, $
 zline=q_zline, nosyn=nosyn, noerr=noerr, ormask=ormask, andmask=andmask, $
 psfile=psfile, xrange=passxr, yrange=passyr, noerase=noerase, $
 netimage=netimage, topdir=topdir, EXTRA=KeywordsForSplot

   cspeed = 2.99792458e5
   textcolor = 'green'
   linecolor = 'magenta'
   orcolor = 'yellow'
   andcolor = 'red'

   readspec, plate, fiberid, mjd=mjd, znum=znum, flux=objflux, $
    wave=wave, plug=plug, zans=zans, topdir=topdir, /silent
   if (NOT keyword_set(objflux)) then begin
      print, plate, mjd, fiberid, $
       format='("Spectrum not found for plate=", i4, " MJD=", i5, " fiber=", i3)'
      return
   endif
   if (NOT keyword_set(noerr)) then $
    readspec, plate, fiberid, mjd=mjd, flerr=objerr, topdir=topdir, /silent
   if (NOT keyword_set(nosyn)) then $
    readspec, plate, fiberid, mjd=mjd, znum=znum, synflux=synflux, $
     topdir=topdir, /silent
   if (keyword_set(ormask)) then $
    readspec, plate, fiberid, mjd=mjd, ormask=ormask, topdir=topdir, /silent
   if (keyword_set(andmask)) then $
    readspec, plate, fiberid, mjd=mjd, andmask=andmask, topdir=topdir, /silent
   if (keyword_set(zans) AND keyword_set(q_zline)) then $
    readspec, plate, fiberid, mjd=mjd, zline=zline, lineflux=lineflux, $
     topdir=topdir, /silent

   if (keyword_set(nsmooth)) then begin
      if (nsmooth GT 1) then begin
         objflux = smooth(objflux, nsmooth)
         if (keyword_set(synflux)) then $
          synflux = smooth(synflux, nsmooth)
         if (keyword_set(lineflux)) then $
          lineflux = smooth(lineflux, nsmooth)
      endif
   endif

   primtarget = sdss_flagname('TARGET', plug.primtarget, /concat)
   sectarget = sdss_flagname('TTARGET', plug.sectarget, /concat)

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
   if (keyword_set(psfile)) then begin
      djs_plot, wave, objflux, xrange=xrange, yrange=yrange, $
       xtitle='Wavelength [Ang]', $
       ytitle=TeXtoIDL('Flux [10^{-17} erg/cm/s/Ang]'), $
       title=title, charsize=csize, _EXTRA=KeywordsForSplot, /xstyle, /ystyle
      if (NOT keyword_set(noerr)) then $
       djs_oplot, wave, objerr, color='red', _EXTRA=KeywordsForSplot
      if (keyword_set(synflux)) then $
       djs_oplot, wave, synflux, color='blue', _EXTRA=KeywordsForSplot, lw=2
   endif else begin
      if (NOT keyword_set(noerase)) then $
       splot, wave, objflux, xrange=xrange, yrange=yrange, $
        xtitle='Wavelength [Ang]', $
        ytitle=TeXtoIDL('Flux [10^{-17} erg/cm/s/Ang]'), $
        title=title, charsize=csize, _EXTRA=KeywordsForSplot $
      else $
       soplot, wave, objflux, _EXTRA=KeywordsForSplot
      if (NOT keyword_set(noerr)) then $
       soplot, wave, objerr, color='red', _EXTRA=KeywordsForSplot
      if (keyword_set(synflux)) then $
       soplot, wave, synflux, color='blue', _EXTRA=KeywordsForSplot, lw=2
   endelse

   xpos = 0.9 * !x.window[0] + 0.1 * !x.window[1]
   dypos = 0.05 * (!y.window[0] - !y.window[1])
   ypos = !y.window[1] + 1.5 * dypos

   if (keyword_set(zans)) then begin
      cz = zans.z * cspeed
      if (abs(cz) LT 3000) then $
        zstring = '  cz=' + string(cz,format='(f6.0)') + ' km/s' $
       else $
        zstring = '  z=' + string(zans.z,format='(f8.5)')
      if (zans.zwarning NE 0) then $
;       zstring = zstring + '  ZWARNING=' + strtrim(string(zans.zwarning),2) + ''
       zstring = zstring + ' (' $
        + sdss_flagname('ZWARNING', zans.zwarning, /concat) + ')'
      if (keyword_set(znum)) then $
       zstring = zstring + ' (fit #' + strtrim(string(znum),2) + ')'

      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, zans.class + ' ' + zans.subclass + zstring, $
        charsize=csize, color=djs_icolor(textcolor), /normal $
      else $
       sxyouts, xpos, ypos, zans.class + ' ' + zans.subclass + zstring, $
        charsize=csize, color=textcolor, /normal

      ypos = ypos + dypos

      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, $
        TeXtoIDL('X^2_r =' + strtrim(string(zans.rchi2, format='(f6.2)'),2)), $
        charsize=csize, color=djs_icolor(textcolor), /normal $
      else $
       sxyouts, xpos, ypos, $
        TeXtoIDL('X^2_r =' + strtrim(string(zans.rchi2, format='(f6.2)'),2)), $
        charsize=csize, color=textcolor, /normal
   endif

   if (keyword_set(lineflux)) then begin
      if (keyword_set(psfile)) then $
       djs_oplot, wave, lineflux, color=linecolor, _EXTRA=KeywordsForSplot, lw=2 $
      else $
       soplot, wave, lineflux, color=linecolor, _EXTRA=KeywordsForSplot, lw=2

      linewave = zline.linewave * (1 + zline.linez)
      ; Convert line sigma from km/sec to Angstroms
      linesigma = linewave * zline.linesigma / cspeed
      linepeak = zline.linecontlevel + zline.linearea / (sqrt(2*!pi) * linesigma)
      for iline=0, n_elements(zline)-1 do begin
         if (zline[iline].linearea_err GT 0) then begin
            if (keyword_set(psfile)) then $
             xyouts, linewave[iline], linepeak[iline], $
              '  '+zline[iline].linename, orient=90, $
              charsize=0.75*csize, color=djs_icolor(linecolor) $
            else $
             sxyouts, linewave[iline], linepeak[iline], $
              '  '+zline[iline].linename, orient=90, $
              charsize=0.75*csize, color=linecolor
         endif
      endfor
   endif

   if (keyword_set(ormask)) then begin
      plotspec_mask, wave, ormask, psfile=psfile, $
       psym=1, symsize=0.6, color=orcolor, nolabel=keyword_set(andmask)
   endif

   if (keyword_set(andmask)) then begin
      plotspec_mask, wave, andmask, psfile=psfile, $
       psym=6, symsize=0.6, color=andcolor
   endif

   if (keyword_set(primtarget)) then begin
      ypos = ypos + dypos
      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, $
        charsize=csize, color=djs_icolor(textcolor), /normal $
      else $
       sxyouts, xpos, ypos, 'PRIMTARGET = ' + primtarget, $
        charsize=csize, color=textcolor, /normal
   endif

   if (keyword_set(sectarget)) then begin
      ypos = ypos + dypos
      if (keyword_set(psfile)) then $
       xyouts, xpos, ypos, 'SECTARGET = ' + sectarget, $
        charsize=csize, color=djs_icolor(textcolor), /normal $
      else $
       sxyouts, xpos, ypos, 'SECTARGET = ' + sectarget, $
        charsize=csize, color=textcolor, /normal
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
 zline=zline, nosyn=nosyn, noerr=noerr, ormask=ormask, andmask=andmask, $
 psfile=psfile, xrange=xrange, yrange=yrange, noerase=noerase, $
 netimage=netimage, zwarning=zwarning, topdir=topdir, _EXTRA=KeywordsForSplot

   if (n_params() LT 1) then begin
      print, 'Syntax - plotspec, plate, [ fiberid, mjd=, znum=, nsmooth=, $'
      print, '         /zline, /nosyn, /noerr, /ormask, /andmask, $'
      print, '         psfile=, xrange=, yrange=, /noerase, $'
      print, '         /netimage, /zwarning, _EXTRA=KeywordsForSplot'
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
         readspec, plate[iplate], mjd=mjd1, topdir=topdir, /silent
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

      for iplate=0, nplate-1 do begin
         readspec, plate[iplate], mjd=mjd[iplate], $
          topdir=topdir, /silent, zans=zans
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
   ; Set FIBERID to [1,...,640] (for each plate) if not set.
   ;
   ; If writing to a PostScript file, then all plots are in the same file
   ; either if PSFILE is that file name, or if FIBERID is not specified
   ; (and then all 640 spectra are being plotted).

   if (NOT keyword_set(fiberid)) then begin
      fiberid = reform((lindgen(640L)+1L) # replicate(1L,nplate), 640L*nplate)
      platelist = reform(replicate(1L,640L) # plate, 640L*nplate)
      mjdlist = reform(replicate(1L,640L) # mjd, 640L*nplate)
      if (keyword_set(psfile)) then begin
         q_onefile = 1
         psfilename = string(plate[0], mjd[0], $
          format='("spec-",i4.4,"-",i5.5,".ps")')
      endif
   endif

   if (min(fiberid) LT 1 OR max(fiberid) GT 640) then begin
      print, 'Invalid FIBERID (must be between 1 and 640)'
      return
   endif

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

   ifiber = 0
   while (ifiber LT nfiber) do begin

      ;----------
      ; Open the PostScript file if appropriate

      if (keyword_set(psfile)) then begin
         if (NOT keyword_set(q_onefile)) then $
          psfilename = string(platelist[ifiber], mjdlist[ifiber], $
           fiberid[ifiber], format='("spec-",i4.4,"-",i5.5,"-",i3.3,".ps")')

         if (NOT keyword_set(q_onefile) OR ifiber EQ 0) then begin
            dfpsplot, psfilename, /color, /square
         endif
      endif

      plotspec1, platelist[ifiber], fiberid[ifiber], mjd=mjdlist[ifiber], $
       znum=znum, nsmooth=nsmooth, zline=zline, nosyn=nosyn, noerr=noerr, $
       ormask=ormask, andmask=andmask, psfile=psfile, $
       xrange=xrange, yrange=yrange, noerase=noerase, netimage=netimage, $
       topdir=topdir, _EXTRA=KeywordsForSplot

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
            print, '      q=quit (and enter interactive mode for this plot)'
            print, '      s=change boxcar smoothing' + sstring
            print, '      x=change X plotting range'
            print, '      y=change Y plotting range'
            print, '      z=change which PCA-fit to overplot'
            print, '      any other key=forward'

            cc = strupcase(get_kbrd(1))
            case cc of
            'B': ifiber = (ifiber - 1) > 0
            'P': begin
                    read, plate, mjd, prompt='Enter new plate and MJD (enter 0 for unknown MJD): '
                    if (NOT keyword_set(mjd)) then $
                     readspec, plate, mjd=mjd, topdir=topdir, /silent
                    if (NOT keyword_set(mjd)) then begin
                       print, 'MJD not found for plate ', plate
                       !quiet = quiet
                       return
                    endif
                    platelist = replicate(plate,640)
                    mjdlist = replicate(mjd,640)
                 end
            'F': begin
                    read, newfiber, prompt='Enter new fiber number: '
                    nfiber = 640
                    fiberid = lindgen(nfiber) + 1L
                    ifiber = ((long(newfiber)-1) > 0) < (nfiber-1)
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
                    read, znum, prompt='Enter 1=best redshift, 2=2nd best, ...: '
                    znum = long(znum) > 0
                 end
            else: ifiber = ifiber + 1
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
