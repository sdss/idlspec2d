;+
; NAME:
;   slithistory
;
; PURPOSE:
;   Plot the history of spectro slit-head positions based upon plSlitpos files.
;
; CALLING SEQUENCE:
;   slithistory, [ mjdrange= ]
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
;   A single PostScript file is created "Slithistory-$MJDSTART-$MJDEND.ps",
;   with one page per cartridge.  The fiber bundle spacing is plotted
;   in units of the fiber-mapper stepper motor steps.
;
; EXAMPLES:
;   Make plots of history of slit-head positions for all time:
;     IDL> slithistory
;
;   Make plots of history of slit-head positions between MJD 52000 and MJD 52020:
;     IDL> setenv, 'SPECLOG_DIR=/scr/spectro1/spectro/scan'
;     IDL> slithistory, mjdrange=[52000,52020]
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
;   19-Nov-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro slithistory, mjdrange=mjdrange

   astrolog = getenv('SPECLOG_DIR')
   if (NOT keyword_set(astrolog)) then astrolog = '/astrolog'
   splog, 'Searching astrolog directory SPECLOG_DIR=' + astrolog

   searchname = filepath('plSlitpos-????-?????-??.par', $
    root_dir=astrolog, subdir='?????')
   slitfiles = findfile(searchname, count=nfile)
   splog, 'Found ', nfile, ' plSlitpos files'
   if (nfile EQ 0) then begin
      splog, 'No matches to ' + searchname
      return
   endif

   ; Sort these by MJD
   mjdvec = long( strmid(fileandpath(slitfiles),15,5) )
   isort = sort(mjdvec)
   slitfiles = slitfiles[isort]
   mjdvec = mjdvec[isort]

   if (keyword_set(mjdrange)) then begin
      if (n_elements(mjdrange) NE 2) then begin
         splog, 'MJDRANGE must be 2-element vector, e.g. MJDRANGE=[52000,52100]'
         return
      endif
      indx = where(mjdvec GE mjdrange[0] AND mjdvec LE mjdrange[1], nfile)
      if (nfile EQ 0) then begin
         splog, 'No plSlitpos files within specified MJD range'
         return
      endif
      splog, 'Trimming to ', nfile, ' files within specified MJD range'
      slitfiles = slitfiles[indx]
   endif

   cartvec = lonarr(nfile)
   slitarr = fltarr(nfile,640)

   plotfile = string(min(mjdvec), max(mjdvec), $
    format='("Slithistory-",i5.5,"-",i5.5,".ps")')
   dfpsplot, plotfile, /color

   for ifile=0, nfile-1 do begin
      splog, 'Reading file ', ifile+1, ' of ', nfile, ': ', $
       fileandpath(slitfiles[ifile])
      yanny_read, slitfiles[ifile], pp, hdr=hdr, /anonymous
      cartvec[ifile] = yanny_par(hdr,'cartridgeId')
;      mjdvec[ifile] = yanny_par(hdr,'fscanMJD')
      slitpos = *pp[0]
      yanny_free, pp
      slitarr[ifile,*] = slitpos.motorpos
   endfor

   !p.multi = [0,1,2]
   yrange = [140, 260]
   psym = -4

   for cartid=1, 9 do begin
      ifile = where(cartvec EQ cartid, nmatch)
      if (nmatch GT 0) then begin
         for specid=1, 2 do begin
            splog, 'Generating plots for cartridge #', cartid, $
             ' spectro-', specid

            ; Select the fiber numbers at which there is a fiber bundle space
            fibernum = (lindgen(15)+1)*20 + (specid-1)*320
            nbund = n_elements(fibernum)

            bundspace = fltarr(nmatch,nbund)
            bundmean = fltarr(nbund)
            for ibund=0, nbund-1 do begin
               pos1 = slitarr[ifile,fibernum[ibund]-1]
               pos2 = slitarr[ifile,fibernum[ibund]]
               bundspace[*,ibund] = (pos2 - pos1) * (pos1 NE 0) * (pos2 NE 0)
               igood = where(bundspace[*,ibund] NE 0, ngood)
               if (ngood GT 0) then bundmean = mean(bundspace[igood,ibund])
            endfor

            xrange = minmax(mjdvec[ifile]) + [-30, 30]
            title='Slithead Gap History Cartridge #' $
             + string(cartid,format='(I1)') + ' Spectro-' $
             + string(specid,format='(I1)')
            colorvec = ['default','red','green','blue','magenta','cyan']
            for ibund=0, nbund-1 do begin
               thiscolor = colorvec[ibund MOD n_elements(colorvec)]
               if (ibund EQ 0) then $
                djs_plot, [0,1], [0,1], /nodata, $
                 xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
                 xtickformat='(i5)', charsize=1.2, $
                 xtitle='MJD', ytitle='Gap size [motor steps]', title=title
               ii = where(bundspace[*,ibund] NE 0) > 0
               djs_oplot, [mjdvec[ifile[ii]]], [bundspace[ii,ibund]], $
                psym=psym, color=thiscolor
               ilast = ii[n_elements(ii)-1]
               djs_xyouts, mjdvec[ifile[ilast]], bundspace[ilast,ibund], $
                ' '+strtrim(string(fibernum[ibund]),2), color=thiscolor
            endfor
         endfor
      endif
   endfor

   dfpsclose
   !p.multi = 0

   return
end
;------------------------------------------------------------------------------
