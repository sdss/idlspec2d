;+
; NAME:
;   plugmap_tycho
;
; PURPOSE:
;   Plot bright Tycho stars on a plate
;
; CALLING SEQUENCE:
;   plugmap_tycho, plugfile, [ matchdist=, mlimit= ]
;
; INPUTS:
;   plugfile   - Name of plugmap file (either plPlugMapP or plPlugMapM)
;
; OPTIONAL INPUTS:
;   matchdist  - Match distance for Tycho stars; default to 5./3600 degrees
;   mlimit     - Magnitude limit for Tycho stars in VTMAG; default to 15.0
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Plot bright Tycho stars on plates, for the purpose of having these
;   fibers not be plugged at APO.  In particular, this procedure was
;   written to deal with the two SEGUE plates 2333 and 2338, as per
;   discussion in the mailing lists on 31 Oct 2005.
;
; EXAMPLES:
;   IDL> plugmap_tycho, 'plPlugMapP-2333.par', mlimit=10.5
;   IDL> plugmap_tycho, 'plPlugMapP-2338.par', mlimit=11.5
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   djs_oplot
;   djs_plot
;   splog
;   tycho_read()
;   yanny_par()
;   yanny_readone()
;
; INTERNAL SUPPORT ROUTINES:
;   tycho_plot_unplugged
;
; REVISION HISTORY:
;   31-Oct-2005  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
; This internal routine is a modified version of PLOT_UNPLUGGED from
; EVILMAP in the idlmapper product.
pro tycho_plot_unplugged, plugdat, plotfile=plotfile, title=title, $
 ibright=ibright

   iobject = where(plugdat.holetype EQ 'OBJECT', nobject)
   iguide = where(plugdat.holetype EQ 'GUIDE' $
    OR plugdat.holetype EQ 'COHERENT_SKY', nguide)
   iquality = where(plugdat.holetype EQ 'QUALITY', nquality)
   itrap = where(plugdat.holetype EQ 'LIGHT_TRAP', ntrap)

   pobj = plugdat[iobject]

   if (keyword_set(plotfile)) then $
    dfpsplot, plotfile, /square, /color

   plot, pobj.xfocal, pobj.yfocal, psym=1, symsize=0.5, $
     xrange=[-350,350], yrange=[-350,350], xstyle=1, ystyle=1, $
     title=title, xtitle='xFocal [mm]', ytitle='yFocal [mm]'
   oplot, [-320], [320], psym=1, symsize=0.5
   xyouts, [-320], [320], '  Plugged fiber'
   if (ibright[0] NE -1) then begin
      djs_oplot, [plugdat[ibright].xfocal], [plugdat[ibright].yfocal], $
       psym=2, symsize=2, color='red'
      nbright = n_elements(ibright)
   endif else nbright = 0
   djs_oplot, [-320], [300], psym=2, symsize=2, color='red'
   xyouts, [-320], [300], '  Bright Tycho stars (' + strtrim(nbright,2) + ')'

   if (nguide GT 0) then begin
      usersym, cos(findgen(21)*!pi/10), sin(findgen(21)*!pi/10) ; open circ
      djs_oplot, [plugdat[iguide].xfocal], [plugdat[iguide].yfocal], $
       psym=8, symsize=3
      xyouts, [plugdat[iguide].xfocal], [plugdat[iguide].yfocal-4], $
       strtrim(string(plugdat[iguide].fiberid),2), alignment=0.5, charsize=1
   endif
   oplot, [-320], [280], psym=8, symsize=3
   xyouts, [-320], [280], '  GUIDE'

   if (ntrap GT 0) then begin
      djs_oplot, [plugdat[itrap].xfocal], [plugdat[itrap].yfocal], $
       psym=6, symsize=1.5
   endif
   oplot, [-320], [260], psym=6, symsize=1.5
   xyouts, [-320], [260], '  LIGHT TRAP'

   if (nquality GT 0) then begin
      usersym, cos(findgen(21)*!pi/10), sin(findgen(21)*!pi/10), /fill ; circ
      djs_oplot, [plugdat[iquality].xfocal], [plugdat[iquality].yfocal], $
       psym=8, symsize=1
   endif
   oplot, [-320], [240], psym=8, symsize=1
   xyouts, [-320], [240], '  QUALITY'

   if (keyword_set(plotfile)) then dfpsclose
   return
end
;------------------------------------------------------------------------------
pro plugmap_tycho, plugfile, matchdist=matchdist1, mlimit=mlimit1

   if (keyword_set(matchdist1)) then matchdist = matchdist1 $
    else matchdist = 5./3600
   if (keyword_set(mlimit1)) then mlimit = mlimit1 $
    else mlimit = 15.

   plugdat = yanny_readone(plugfile, hdr=hdr)
   if (NOT keyword_set(plugdat)) then begin
      splog, 'No plugmap file found ', plugfile
      return
   endif
   iobj = where(strmatch(plugdat.holetype,'OBJECT'))

   racen = yanny_par(hdr, 'raCen')
   deccen = yanny_par(hdr, 'decCen')
   plotfile = 'tycho-' + repstr(fileandpath(plugfile), '.par', '.ps')
   title = 'TYCHO STARS ON PLATE ' + strtrim(yanny_par(hdr, 'plateId'),2) $
    + ' V < ' + string(mlimit,format='(f5.2)')

   tyc = tycho_read(racen=racen, deccen=deccen, radius=1.6)
   spherematch, plugdat[iobj].ra, plugdat[iobj].dec, tyc.ramdeg, tyc.demdeg, $
    matchdist, i1, i2, d12
   ibright = -1
   if (i2[0] NE -1) then begin
      ii = where(tyc[i2].vtmag LE mlimit, ct)
      if (ct GT 0) then begin
         ibright = iobj[i1[ii]]
         for i=0L, ct-1L do $
          splog, 'Bright' $
           + ' RA= ' + string(plugdat[ibright[i]].ra,format='(f8.4)') $
           + ' Dec= ' + string(plugdat[ibright[i]].dec,format='(f8.4)') $
           + ' X= ' + string(plugdat[ibright[i]].xfocal,format='(f5.0)') $
           + ' Y= ' + string(plugdat[ibright[i]].xfocal,format='(f5.0)') $
           + ' Vmag= ' + string(tyc[i2[ii[i]]].vtmag,format='(f4.1)')
      endif
      splog, 'Number of bright stars = ', ct
   endif

   tycho_plot_unplugged, plugdat, plotfile=plotfile, title=title, $
    ibright=ibright

   return
end
;------------------------------------------------------------------------------
