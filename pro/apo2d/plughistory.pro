;+
; NAME:
;   plughistory
;
; PURPOSE:
;   Plot the history of unplugged fibers based upon plPlugMapM files.
;
; CALLING SEQUENCE:
;   plughistory, [ mjdrange= ]
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
;   A single PostScript file is created "Plughistory-$MJDSTART-$MJDEND.ps",
;   with one page per cartridge.  The top panel shows the number of fibers
;   plugged as a function of MJD.  The bottom panel shows the fraction
;   of times each fiber is plugged, with < 90%-plugged fibers labelled
;   with their fiber number.
;
; EXAMPLES:
;   Make plots of history of unplugged fibers for all time:
;     IDL> plughistory
;
;   Make plots of history of unplugged fibers between MJD 52000 and MJD 52100:
;     IDL> plughistory, mjdrange=[52000,52100]
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
;   27-Feb-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro plughistory, mjdrange=mjdrange

   astrolog = getenv('ASTROLOG_DIR')
   if (NOT keyword_set(astrolog)) then astrolog = '/astrolog'
   splog, 'Searching astrolog directory ASTROLOG_DIR=' + astrolog

   searchname = filepath('plPlugMapM-????-?????-??.par', $
    root_dir=astrolog, subdir='?????')
   plugfiles = findfile(searchname, count=nfile)
   splog, 'Found ', nfile, ' plPlugMapM files'
   if (nfile EQ 0) then begin
      splog, 'No matches to ' + searchname
      return
   endif

   mjdvec = long( strmid(fileandpath(plugfiles),16,5) )
   if (keyword_set(mjdrange)) then begin
      if (n_elements(mjdrange) NE 2) then begin
         splog, 'MJDRANGE must be 2-element vector, e.g. MJDRANGE=[52000,52100]'
         return
      endif
      indx = where(mjdvec GE mjdrange[0] AND mjdvec LE mjdrange[1], nfile)
      if (nfile EQ 0) then begin
         splog, 'No plPlugMapM files within specified MJD range'
         return
      endif
      splog, 'Trimming to ', nfile, ' files within specified MJD range'
      plugfiles = plugfiles[indx]
   endif

   cartvec = lonarr(nfile)
;   mjdvec = lonarr(nfile)
   plugarr = lonarr(nfile,641)

   plotfile = string(min(mjdvec), max(mjdvec), $
    format='("Plughistory-",i5.5,"-",i5.5,".ps")')
   dfpsplot, plotfile

   for ifile=0, nfile-1 do begin
      splog, 'Reading file ', ifile+1, ' of ', nfile, ': ', $
       fileandpath(plugfiles[ifile])
      yanny_read, plugfiles[ifile], pp, hdr=hdr, /anonymous
      cartvec[ifile] = yanny_par(hdr,'cartridgeId')
;      mjdvec[ifile] = yanny_par(hdr,'fscanMJD')
      plugmap = *pp[0]
      yanny_free, pp
      indx = where(strtrim(plugmap.holetype,2) EQ 'OBJECT', ct)
      if (ct NE 640) then $
       splog, 'Warning: Fewer than 640 objects in ' $
        + fileandpath(plugfiles[ifile])
      fiberid = plugmap[indx].fiberid
      plugarr[ifile,fiberid] = plugarr[ifile,fiberid] + 1
   endfor

   !p.multi = [0,1,2]
   xrange = minmax(mjdvec[where(mjdvec NE 0)]) + [-30,30]
   yrange = [620,645]

   for cartid=1, 9 do begin
      ifile = where(cartvec EQ cartid, nmatch)
      if (nmatch GT 0) then begin
         splog, 'Generating plots for cartridge #', cartid

         plot, mjdvec[ifile], total(plugarr[ifile,1:640], 2), $
          psym=2, symsize=0.5, $
          xrange=xrange, /xstyle, yrange=yrange, /ystyle, $
          xtickformat='(i5)', $
          xtitle='MJD', ytitle='Number plugged fibers', $
          title='Plugging History Cartridge #' + strtrim(string(cartid),2)
         oplot, !x.crange, [640,640]
         xyouts, !x.crange[0], 642, $
          '  Number of pluggings = ' + strtrim(string(nmatch),2)

         meannum = total(plugarr[ifile,1:640], 1) / nmatch
         plot, lindgen(640)+1, meannum, psym=10, $
          xrange=[-20,660], /xstyle, yrange=[-0.1,1.2], /ystyle, $
          xtitle='Fiber number', ytitle='Fraction times plugged'
         ; Overplot wherever this fraction is less than 90%
         ibad = where(meannum LT 0.90, nbad)
         for ii=0, nbad-1 do $
          xyouts, ibad[ii]+1, 1.0, string(ibad[ii]+1, format='(i4)'), $
           orient=90, charsize=0.75
      endif
   endfor

   dfpsclose
   !p.multi = 0

   return
end
;------------------------------------------------------------------------------
