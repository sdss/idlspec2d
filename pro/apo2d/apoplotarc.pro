;+
; NAME:
;   apoplotarc
;
; PURPOSE:
;   Routine for plotting arc spectra from the Son-of-Spectro outputs at APO.
;
; CALLING SEQUENCE:
;   apoplotarc, expnum, [ camname=, mjd=, everyn=, psfile=, _EXTRA= ]
;
; INPUTS:
;   expnum     - Exposure number
;
; OPTIONAL INPUTS:
;   camname    - Camera name; default to 'r1'
;   mjd        - MJD; must be set if this exposure is not in the most
;                recent MJD directory
;   everyn     - Plot every EVERYN-th spectrum; default to 40, which will
;                plot fiber numbers 1, 41, 81, ...281.
;   psfile     - If set, then send plot to a PostScript file instead of
;                to the SPLOT interactive widget.  The PostScript file name
;                can be set explicitly, e.g. with PSFILE='test.ps'.  Or if
;                you simply set this as a flag, e.g. with /PSFILE, then the
;                default file name is arc-pppp-cc-eeeeeeee.ps,
;                where pppp=plate number, cc=camname, eeeeeeee=exposure number.
;   _EXTRA     - Kewords for SPLOT, such as XRANGE, YRANGE, THICK.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The arc spectrum is plotted as a red line.  The individual
;   arc values at each pixel in the arc fibers are plotted as points.
;
;   The Son-of-Spectro outputs are first read from the file:
;     $SPECTROLOG_DIR/$MJD/logsheet-$MJD.fits
;   This then points to the other files that are read:
;     $SPECTROLOG_DIR/$MJD/fflat-$MJD-$PLATE-$EXPNUM-$CAMERA.fits
;     $SPECTROLOG_DIR/$MJD/wset-$MJD-$PLATE-$EXPNUM-$CAMERA.fits
;     $SPECTROLOG_DIR/$MJD/sci-$PLATE-$CAMERA-$EXPNUM.fits
;
;   If $SPECTROLOG_DIR is not set, then it is assumed to be
;     /data/spectro/spectrologs
;
; EXAMPLES:
;   Plot the arc spectrum for exposure number 9437 taken in the most
;   recent night's data in the r1-camera, using the SPLOT plotting tool:
;     IDL> apoplotarc, 9437
;   The mouse buttons will zoom in (left), recenter (center), or zoom out
;   (right).  The frame can be saved as a PostScript file by selecting
;   File->WriteEPS from the left-hand corner. 
;
;   Make the same plot, but change the Y limits and write a PostScript file:
;     IDL> apoplotarc, 9437, yrange=[0,1000], /psfile
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_oplot
;   djs_plot
;   get_mjd_dir()
;   soplot
;   splot
;
; REVISION HISTORY:
;   20-Nov-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro apoplotarc, expnum, camname=camname, mjd=mjd, everyn=everyn, $
 psfile=psfile, xrange=xrange, yrange=yrange, $
 _EXTRA=KeywordsForSplot

   if (n_params() LT 1) then begin
      print, 'Syntax - apoplotarc, expnum, [ camname=, mjd=, everyn=, $'
      print, '         psfile=, xrange=, yrange=, $'
      print, '         _EXTRA=KeywordsForSplot'
      return
   endif

   quiet = !quiet
   !quiet = 1

   ;----------
   ; Set defaults

   if (n_elements(expnum) NE 1) then $
    message, 'EXPNUM must be a scalar'
   if (NOT keyword_set(camname)) then camname = 'r1'
   if (camname NE 'b1' AND camname NE 'b2' $
    AND camname NE 'r1' AND camname NE 'r2') then begin
      print, 'Invalid CAMNAME ', camname
      !quiet = quiet
      return
   endif
   if (NOT keyword_set(everyn)) then everyn = 40
   nplot = long(320 / everyn)
   fibernum = lindgen(nplot) * everyn

   ;----------
   ; If MJD is not specified, then find the most recent MJD for output files

   spectrolog_dir = getenv('SPECTROLOG_DIR')
   if (NOT keyword_set(spectrolog_dir)) then $
    spectrolog_dir = '/data/spectro/spectrologs'

   if (NOT keyword_set(mjd)) then begin
      mjdlist = get_mjd_dir(spectrolog_dir, mjstart=1, mjend=99999, mjd='?????')
      mjd = (reverse(mjdlist[sort(mjdlist)]))[0]
      splog, 'Selecting MJD=', mjd, ' (override this with MJD keyword)'
   endif

   ;----------
   ; Find the wset file, from which we'll read the wavelengths
   ; and the arc spectra

   expstr = string(expnum, format='(i8.8)')
   mjdstr = string(mjd, format='(i5.5)')
   wsetfile = filepath('wset-'+mjdstr+'-????-'+expstr+'-'+camname+'.fits', $
    root_dir=spectrolog_dir, subdir=mjdstr)
   thisfile = (findfile(wsetfile, count=ct))[0]
   if (ct EQ 0) then begin
      splog, 'File not found ' + wsetfile
      !quiet = quiet
      return
   endif
   platestr = strmid(fileandpath(thisfile), 11, 4)

   wset = mrdfits(thisfile, 1)
   flux = mrdfits(thisfile, 2)
   if (NOT keyword_set(wset) OR NOT keyword_set(flux)) then begin
      splog, 'Missing data in file ' + thisfile
      !quiet = quiet
      return
   endif

   traceset2xy, wset, xx, loglam

   ;----------
   ; Make plots

   title = 'Arc for Plate ' + platestr + '  MJD=' + mjdstr $
    + '  ' + camname+'-'+expstr
   xtitle = 'Vacuum Wavelength [Ang]'
   ytitle = 'Flux [electrons]'
   if (NOT keyword_set(xrange)) then xrange = 10.d0^minmax(loglam)
   if (NOT keyword_set(yrange)) then yrange = minmax(flux)
   colorvec = ['default','red','green','blue','magenta','cyan', $
    'yellow']
   ncolor = n_elements(colorvec)

   if (keyword_set(psfile)) then begin
      if (size(psfile,/tname) EQ 'STRING') then psfilename = psfile $
       else psfilename = string(platestr, camname, expnum, $
        format='("arc-",a4,"-",a2,"-",i8.8,".ps")')
      splog, 'Writing to file ' + psfilename
      dfpsplot, psfilename, /color, /square
      djs_plot, 10.d0^loglam[*,fibernum[0]], flux[*,fibernum[0]], $
       xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
       xtitle=xtitle, ytitle=ytitle, title=title, charsize=1.3, $
       color=colorvec[0], _EXTRA=KeywordsForSplot
      for i=0, nplot-1 do $
       djs_oplot, 10.d0^loglam[*,fibernum[i]], flux[*,fibernum[i]], $
        color=colorvec[i MOD ncolor], _EXTRA=KeywordsForSplot
      dfpsclose
   endif else begin
      splot, 10.d0^loglam[*,fibernum[0]], flux[*,fibernum[0]], $
       xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
       xtitle=xtitle, ytitle=ytitle, title=title, charsize=1.3, $
       color=colorvec[0], _EXTRA=KeywordsForSplot
      for i=0, nplot-1 do $
       soplot, 10.d0^loglam[*,fibernum[i]], flux[*,fibernum[i]], $
        color=colorvec[i MOD ncolor], _EXTRA=KeywordsForSplot
   endelse

   !quiet = quiet
   return
end
;------------------------------------------------------------------------------
