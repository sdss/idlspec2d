;+
; NAME:
;   platecompare
;
; PURPOSE:
;   Interactive comparison of descrepant redshifts for the same objects.
;
; CALLING SEQUENCE:
;   platecompare, plate, [mjd=, psfile= ]
;
; INPUTS:
;   plate       - Plate number(s)
;
; OPTIONAL INPUTS:
;   mjd         - MJD for each plate number.  If specified, then this
;                 must have one MJD per plate number.  If not specified,
;                 then all MJD's associated with each plate are read.
;                 That list of MJD's comes from the PLATELIST command.
;   psfile     - If set, then send plot to a PostScript file instead of
;                to the SPLOT interactive widget.  The PostScript file name
;                can be set explicitly, e.g. with PSFILE='test.ps'.  Or if
;                you simply set this as a flag, e.g. with /PSFILE, then the
;                default file name is platecompare-pppp.ps,
;                where pppp is the first plate number.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Duplicate objects are found by matching positions to within 1 arcsec.
;
; EXAMPLES:
;   Compare data from different pluggings of the same plate:
;   IDL> platecompare, [406,406,406], mjd=[51817,51869,51876]
;
;   Compare date from two plates, 360 and 362, which are actually the
;   same tile on the sky:
;   IDL> platecompare, [360,360,362], mjd=[51780,51816,51999]
;
; BUGS:
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_angle_group
;   djs_icolor()
;   djs_oplot
;   djs_plot
;   platelist
;   readspec
;   sdss_flagname()
;
; REVISION HISTORY:
;   14-Aug-2001  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platecompare, plate, mjd=mjd, psfile=psfile

   if (n_params() LT 1) then begin
      print, 'Syntax - platecompare, plate, [mjd= ]'
      return
   endif

   if (keyword_set(mjd)) then begin
      if (n_elements(plate) EQ 1) then begin
         platevec = replicate(plate, n_elements(mjd))
      endif else if (n_elements(mjd) NE n_elements(plate)) then begin
         print, 'Number of elements in PLATE and MJD do not agree.'
         return
      endif else begin
         platevec = plate
      endelse
      mjdvec = mjd
   endif else begin
      for iplate=0, n_elements(plate)-1 do begin
         platelist, plist=plist
         ii = where(plist.plate EQ plate[iplate], ni)
         if (ni GT 0) then begin
            if (NOT keyword_set(platevec)) then $
             platevec = replicate(plate[iplate],ni) $
            else $
             platevec = [platevec, replicate(plate[iplate],ni)]
            if (NOT keyword_set(mjdvec)) then mjdvec = plist[ii].mjd $
             else platevec = [platevec, plist[ii].mjd]
         endif
      endfor
      if (NOT keyword_set(platevec)) then begin
         print, 'No matches found for plate ', plate
         return
      endif
   endelse

   cspeed = 2.99792458e5

   ; Read the redshift files
   readspec, platevec, mjd=mjdvec, zans=zansall, plug=plugall

   ; Group duplicate objects together
   ngroup = djs_angle_group(zansall.plug_ra, zansall.plug_dec, 1./3600., $
    gstart=gstart, gcount=gcount, gindx=gindx)
   print, 'Number of duplicate objects: ', ngroup
   if (ngroup EQ 0) then return

   ; Open the PostScript file
   if (keyword_set(psfile)) then begin
      if (size(psfile,/tname) EQ 'STRING') then psfilename = psfile $
       else psfilename = string(platevec[0], $
        format='("platecompare-",i4.4,".ps")')
      dfpsplot, psfilename, /color
   endif

   for igroup=0, ngroup-1 do begin
      thiszans = zansall[gindx[gstart[igroup]:gstart[igroup]+gcount[igroup]-1]]
      thisplug = plugall[gindx[gstart[igroup]:gstart[igroup]+gcount[igroup]-1]]

      vdiff = (max(thiszans.z) - min(thiszans.z)) * cspeed
      vbad = ((vdiff GT 250.) AND (thiszans[0].class EQ 'GALAXY')) $
          OR ((vdiff GT 500.) AND (thiszans[0].class EQ 'QSO')) $
          OR ((vdiff GT 50.) AND (thiszans[0].class EQ 'STAR'))
      vbad = vbad AND ((thiszans[0].zwarning AND 1) EQ 0) ; Not SKY fiber

      xrange = [3600, 9300]
      textcolor = 'green'

      if (vbad) then begin
         !p.multi = [0, 1, gcount[igroup]]
         if (thiszans[0].sn_median LT 1) then nsmooth = 7 $
          else if (thiszans[0].sn_median LT 3) then nsmooth = 5 $
          else if (thiszans[0].sn_median LT 5) then nsmooth = 3 $
          else nsmooth = 1

         for ii=0, gcount[igroup]-1 do begin
            primtarget = sdss_flagname('TARGET', thisplug[ii].primtarget, $
             /concat)
            sectarget = sdss_flagname('TTARGET', thisplug[ii].sectarget, $
             /concat)
            print, thiszans[ii].class, thiszans[ii].z
            title = string(thiszans[ii].plate, thiszans[ii].mjd, $
             thiszans[ii].fiberid, $
             format='("Plate ", i4, "-", i5, " Fiber ", i3)')
            cz = thiszans[ii].z * cspeed
            zstring = string(strtrim(thiszans[ii].class), $
             strtrim(' '+thiszans[ii].subclass), $
             format='(a, a)')
            if (abs(cz) LT 3000) then $
             zstring = zstring + string(cz, format='(" z=", f6.0)') + ' km/s' $
            else $
             zstring = zstring + string(thiszans[ii].z, format='(" z=", f8.5)')
            if (thiszans[ii].zwarning NE 0) then $
             zstring = zstring + $
              '  ZWARNING=' + strtrim(string(thiszans[ii].zwarning),2)
            readspec, thiszans[ii].plate, thiszans[ii].fiberid, $
             mjd=thiszans[ii].mjd, wave=wave, flux=objflux, synflux=synflux
            ytitle = 'F_\lambda'
            if (nsmooth GT 1) then begin
               objflux = smooth(objflux, nsmooth)
               synflux = smooth(synflux, nsmooth)
               ytitle = ytitle $
                + string(nsmooth, format='(" (nsmooth=", i2, ")")')
            endif
            yrange = minmax(synflux)
            if (yrange[0] EQ yrange[1]) then yrange = minmax(objflux)
            ymin = (1.2 * yrange[0] - 0.2 * yrange[1]) < 0
            ymax = -0.2 * yrange[0] + 1.2 * yrange[1]
            if (ymax EQ ymin) then ymax = ymin + 1
            yrange = [ymin, ymax]
            djs_plot, wave, objflux, $
             xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
             xtitle=zstring, ytitle=ytitle, title=title, charsize=2.5
            djs_oplot, wave, synflux, color='blue'
            xpos = 1.0 * !x.crange[0] + 0.0 * !x.crange[1]
            ypos = -0.03 * !y.crange[0] + 1.03 * !y.crange[1]
            xyouts, xpos, ypos, primtarget+' '+sectarget, $
             color=djs_icolor(textcolor)
         endfor
         if (NOT keyword_set(psfile)) then $
          cc = strupcase(get_kbrd(1))
      endif
   endfor

   if (keyword_set(psfile)) then dfpsclose

   return
end
;------------------------------------------------------------------------------
