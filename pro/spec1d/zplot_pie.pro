;+
; NAME:
;   zplot_pie
;
; PURPOSE:
;   Plot pie diagram
;
; CALLING SEQUENCE:
;   zplot_pie, [ rarange=, decrange=, zmax=, deltaz=, plotfile ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   rarange    - RA range; default to [130.,250] deg
;   decrange   - Declination range; default to [-0.5,0.5]
;   zmax       - Maximum redshift; default to 0.15
;   deltaz     - Spacing for labeling the redshift; default to 0.05
;   plotfile   - Plot file; default to 'zplot.ps'
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   Plot the pie diagram for the Southern stripe (stripe 82):
;     IDL> zplot, rarange=[-45,50]
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   platelist
;   readspec
;
; REVISION HISTORY:
;   14-May-2007  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro zplot_pie, rarange=rarange, decrange=decrange, zmax=zmax, deltaz=deltaz, $
 plotfile=plotfile

   if (n_elements(rarange) NE 2) then rarange = [130.,250]
   if (n_elements(decrange) NE 2) then decrange = [-0.5,0.5]
   if (NOT keyword_set(zmax)) then zmax = 0.15
   if (NOT keyword_set(deltaz)) then deltaz = 0.05
   if (NOT keyword_set(plotfile)) then plotfile = 'zplot.ps'

   psym = 6
   symsize = 0.25

   ramid = 0.5*(rarange[0]+rarange[1])
   ra0 = ramid - 90.

   platelist, plist=plist
   ii = where(plist.dec GT decrange[0]-3. AND plist.dec LT decrange[1]+3. $
    AND plist.qsurvey AND strtrim(plist.status1d,2) EQ 'Done')
   readspec, plist[ii].plate, mjd=plist[ii].mjd, plug=plug, zans=zans

   xplot = zans.z * cos((plug.ra - ra0) / !radeg)
   yplot = zans.z * sin((plug.ra - ra0) / !radeg)

   ; Select all galaxies in this pie slice
   spec_gal = strtrim(zans.class,2) EQ 'GALAXY' AND zans.zwarning EQ 0
   if (rarange[0] GE 0) then $
    qgoodra = zans.plug_ra GE rarange[0] AND zans.plug_ra LE rarange[1] $
   else $
    qgoodra = zans.plug_ra GE rarange[0]+360 OR zans.plug_ra LE rarange[1]
   iplot = where(spec_gal AND zans.z LT zmax AND qgoodra $
    AND zans.plug_dec GT decrange[0] AND zans.plug_dec LT decrange[1])

   !x.margin = [4,4]
   !y.margin = [2,2]
   thick = 2
   csize = 1.2

   dfpsplot, plotfile, /square
   plot, [0], [0], /nodata, xrange=[-zmax,zmax], yrange=[-zmax,zmax], $
    xstyle=5, ystyle=5
   oplot, xplot[iplot], yplot[iplot], psym=psym, symsize=symsize, $
    thick=thick, charthick=thick

   ravec = rarange[0]-ra0 + 0.1*findgen(10.*(rarange[1]-rarange[0])+1)
   xtmp = [0,zmax * cos(ravec/!radeg),0]
   ytmp = [0,zmax * sin(ravec/!radeg),0]
   oplot, xtmp, ytmp, thick=thick

   dra = 20.
   for thisra = 0., 359.99, dra do $
    if (thisra GE rarange[0] AND thisra LE rarange[1]) then $
     xyouts, 1.02*zmax*cos((thisra-ra0)/!radeg), $
      1.02*zmax*sin((thisra-ra0)/!radeg), $
      textoidl(strtrim(string(thisra,format='(i3)'),2)+'^o'), $
      align=0.5, orient=thisra-ramid, charthick=thick, charsize=csize
   for thisra = 0., 359.99, dra do $
    if (thisra GE rarange[0] AND thisra LE rarange[1]) then $
     oplot, [0.95,1]*zmax*cos((thisra-ra0)/!radeg), $
      [0.95,1]*zmax*sin((thisra-ra0)/!radeg), thick=thick
   dra = 5.
   for thisra = 0., 359.99, dra do $
    if (thisra GE rarange[0] AND thisra LE rarange[1]) then $
     oplot, [0.975,1]*zmax*cos((thisra-ra0)/!radeg), $
      [0.975,1]*zmax*sin((thisra-ra0)/!radeg), thick=thick

   theta = (rarange[0]-ra0)/!radeg
   xyouts, 0.5*deltaz*cos(theta) + 0.06*zmax*sin(theta), $
    0.5*deltaz*sin(theta) - 0.06*zmax*cos(theta), $
    'z=', $
    align=0.5, orient=rarange[0]-ramid+90., charthick=thick, charsize=csize
   for thisz=deltaz, deltaz*floor(zmax/deltaz), deltaz do $
    xyouts, thisz*cos(theta) + 0.06*zmax*sin(theta), $
     thisz*sin(theta) - 0.06*zmax*cos(theta), $
     strtrim(string(thisz,format='(f4.2)'),2), $
     align=0.5, orient=rarange[0]-ramid+90., charthick=thick, charsize=csize
   for thisz=deltaz, deltaz*floor(zmax/deltaz), deltaz do $
    oplot, thisz*cos(theta) + [0,-0.05*zmax*sin(theta)], $
     thisz*sin(theta) + [0,0.05*zmax*cos(theta)], $
     thick=thick
   dfpsclose

   return
end
;------------------------------------------------------------------------------
