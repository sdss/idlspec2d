;+
; NAME:
;   bundlethru
;
; PURPOSE:
;   Measure relative throughput of fiber bundles from flat-fields
;
; CALLING SEQUENCE:
;   bundlethru, [ wrange=, mjdrange=, /first ]
;
; INPUTS:
;   wrange     - Wavelength range; default to [5000,5100] Ang
;   mjdrange   - If set, then limit to this MJD range
;   first      - If set, then use the first rather than last good observation
;                of each plate
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine relies upon the platelist file to find good plates
;   for each cartridge.  Then the spPlate, spCFrame, spFlat, spArc files
;   are read.
;
;   A plot file and a data file are output with names bundlethru-$MJD.ps
;   and bundlethru-$MJD.dat .  The MJD in the name is the lowest MJD used.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   splog
;
; REVISION HISTORY:
;   26-Aug-2011  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro bundlethru, wrange=wrange, mjdrange=mjdrange, first=first

   if (NOT keyword_set(wrange)) then wrange = [5000,5100]

   ; Select the last good observation of each cartridge
   platelist, plist=plist
   if (NOT keyword_set(plist)) then begin
      print, 'No platelist found'
      return
   end
   indx = where(strmatch(plist.platequality,'good*'), ct)
   if (ct EQ 0) then begin
      print, 'No good plates found'
      return
   end
   plist = plist[indx]
   if (keyword_set(mjdrange)) then begin
      indx = where(plist.mjd GE mjdrange[0] AND plist.mjd LE mjdrange[1], ct)
      if (ct EQ 0) then begin
         print, 'No good plates found'
         return
      end
   endif
   if (keyword_set(first)) then $
    plist = plist[sort((plist.cartid * max(plist.mjd)+1) - plist.mjd)] $
   else $
    plist = plist[sort((plist.cartid * max(plist.mjd)+1) + plist.mjd)]
   i = uniq(plist.cartid)
   plist = plist[i]
   ncart = n_elements(plist)

   mjdstr = string(min(plist.mjd),format='(i5.5)')
   logfile = 'bundlethru-'+mjdstr+'.dat'
   plotfile = 'bundlethru-'+mjdstr+'.ps'

   nbundle = 25
   fbundle = fltarr(ncart, 2, nbundle)

   splog, file=logfile, /noname
   dfpsplot, plotfile, /color, /square

   for icart=0, ncart-1 do begin
      for ispec=0, 1 do begin ; loop over spectrographs
         fiberid = 1 + 500*ispec
         cameras = (wrange[0] LT 6000) ? 'b' : 'r'
         readonespec, plist[icart].plate, fiberid, mjd=plist[icart].mjd, $
          cameras=cameras,framehdr=phdr
         hdr1 = *phdr[0]
         ptr_free, phdr
         flatfile = 'spFlat-'+strmid(sxpar(hdr1,'FLATFILE'),4,11)+'.fits.gz'
         arcfile = 'spArc-'+strmid(sxpar(hdr1,'ARCFILE'),4,11)+'.fits.gz'
         flatfile = filepath(flatfile, root_dir=getenv('BOSS_SPECTRO_REDUX'), $
          subdir=[getenv('RUN2D'),string(plist[icart].plate,format='(i4.4)')])
         arcfile = filepath(arcfile, root_dir=getenv('BOSS_SPECTRO_REDUX'), $
          subdir=[getenv('RUN2D'),string(plist[icart].plate,format='(i4.4)')])
         flatimg = mrdfits(flatfile, 0, /silent)
         wset = mrdfits(arcfile, 2, /silent)
         traceset2xy, wset, xx, loglam
         mask = loglam GT alog10(wrange[0]) AND loglam LT alog10(wrange[1])
         ffiber = total(flatimg * mask, 1) / total(mask, 1)
         for ib=0, nbundle-1 do $
          fbundle[icart,ispec,ib] = median(ffiber[ib*20:(ib+1)*20-1])

         ; Renormalize to be relative to the maximum-throughput bundle
         ffiber /= max(fbundle[icart,ispec,*])
         fbundle[icart,ispec,*] /= max(fbundle[icart,ispec,*])

         for ib=0, nbundle-1 do $
          splog, plist[icart].cartid, ispec+1, ib, fbundle[icart,ispec,ib], $
           /noname

         title = 'Cartridge '+strtrim(plist[icart].cartid,2) $
          +'  Slithead '+strtrim(ispec+1,2)
         djs_plot, lindgen(500)+0.5, ffiber, yrange=[0,1.3], psym=10, $
          xtitle='Fiber number', ytitle='Relative throughput', $
          title=title, charsize=1.5, color='red'
         djs_oplot, lindgen(nbundle)*20+10.5, fbundle[icart,ispec,*], $
          psym=6
         for ib=0, nbundle-1 do $
          djs_oplot, ib*20+[1,20], fbundle[icart,ispec,ib]+[0,0]
         for ib=0, nbundle-1 do $
          xyouts, ib*20+10.5, fbundle[icart,ispec,ib]+0.1, $
           strtrim(ib+1,2), charsize=1.4
         xyouts, 50, 0.1, 'From plate '+strtrim(plist[icart].plate,2) $
          +' MJD '+strtrim(plist[icart].mjd,2) $
          +' wave='+strtrim(fix(wrange[0]),2)+'-'+strtrim(fix(wrange[1]),2)
      endfor
   endfor

   plothist, fbundle, bin=0.02, xtitle='Relative throughput', $
    ytitle='Number of bundles', charsize=1.5

   splog, /close
   dfpsclose
end

