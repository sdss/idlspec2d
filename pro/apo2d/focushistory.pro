;+
; NAME:
;   focushistory
;
; PURPOSE:
;   Plot the history of wavelength focus based upon the reduced spPlate files.
;
; CALLING SEQUENCE:
;   focushistory, [ mjdrange= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjdrange   - 2-element vector of plotting range in MJD; default to
;                using all MJDs.
;
; OUTPUT:
;
; COMMENTS:
;   A single PostScript file is created "Focushistory-$MJDSTART-$MJDEND.ps".
;   The plate list must exist for the PLATELIST procedure, and the spPlate
;   files should be in $SPECTRO_DATA.
;
; EXAMPLES:
;   Make plots of history of wavelength focus for all time:
;     IDL> focushistory
;
;   Make plots of history of wavelength focus between MJD 52000 and MJD 52100:
;     IDL> focusistory, mjdrange=[5200,52100]
;
; BUGS:
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   fileandpath()
;   yanhy_free
;   yanny_par()
;   yanny_read
;
; REVISION HISTORY:
;   21-Mar-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro focushistory, mjdrange=mjdrange

   platelist, plist=plist
   if (NOT keyword_set(plist)) then begin
      print, 'No reduced plates found'
      return
   endif

   idone = where(strtrim(plist.status1d,2) EQ 'Done')
   if (idone[0] EQ -1) then begin
      print, 'No reduced plates found'
      return
   endif
   plist = plist[idone]

   if (keyword_set(mjdrange)) then begin
      if (n_elements(mjdrange) NE 2) then begin
         splog, 'MJDRANGE must be 2-element vector, e.g. MJDRANGE=[52000,52100]'
         return
      endif
      indx = where(plist.mjd GE mjdrange[0] AND plist.mjd LE mjdrange[1], nfile)
      if (nfile EQ 0) then begin
         splog, 'No reduced plates within specified MJD range'
         return
      endif
      splog, 'Trimming to ', nfile, ' plates within specified MJD range'
      plist = plist[indx]
   endif

   plist = plist[sort(plist.mjd)] ; Sort by MJD
   nfile = n_elements(plist)
   disparr = fltarr(4,nfile)
   dispmax = fltarr(nfile)

   for ifile=0, nfile-1 do begin
      splog, 'Reading file ', ifile+1, ' of ', nfile, ': ', $
       plist[ifile].plate, '/', plist[ifile].mjd
      readspec, plist[ifile].plate, mjd=plist[ifile].mjd, disp=dispimg
      dims = size(dispimg,/dimens)
      disparr[0,ifile] = median(dispimg[0:dims[0]/2-1,0:dims[1]/2-1])
      disparr[1,ifile] = median(dispimg[dims[0]/2:dims[0]-1,0:dims[1]/2-1])
      disparr[2,ifile] = median(dispimg[0:dims[0]/2-1,dims[1]/2-1:dims[1]-1])
      disparr[3,ifile] = median(dispimg[dims[0]/2:dims[0]-1,dims[1]/2-1:dims[1]-1])
      dispmax[ifile] = max(disparr[*,ifile])
   endfor

   plotfile = string(min(plist.mjd), max(plist.mjd), $
    format='("Focushistory-",i5.5,"-",i5.5,".ps")')
   dfpsplot, plotfile

   xrange = minmax(plist.mjd) + [-30,30]
   yrange = [0.0,2.0]

   plot, plist.mjd, dispmax, psym=4, $
    xrange=xrange, /xstyle, yrange=yrange, /ystyle, $
    xtickformat='(i5)', $
    xtitle='MJD', ytitle='Focus [pix]', $
    title='Focus History'
   oplot, !x.crange, [1,1]

   dfpsclose
stop

   return
end
;------------------------------------------------------------------------------
