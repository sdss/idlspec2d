;+
; NAME:
;   platecompare
;
; PURPOSE:
;   Interactive comparison of descrepant redshifts for the same objects.
;
; CALLING SEQUENCE:
;   platecompare, plate, [mjd=, run2d=, run1d=, psfile= ]
;
; INPUTS:
;   plate       - Plate number(s)
;
; OPTIONAL INPUTS:
;   mjd        - MJD for each plate number.  If specified, then this
;                must have one MJD per plate number.  If not specified,
;                then all MJD's associated with each plate are read.
;                That list of MJD's comes from the PLATELIST command.
;   run2d      - Names of 2D reductions if comparing different versions
;                of the same data.
;   run1d      - Names of 1D redutions if comparing different versions
;                of the same data; must have same number of elements as RUN2D;
;                if not set, then default to the same versions as RUN2D
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
;   Plot any objects that show a descrepency in redshifts between any
;   of the observations, or any objects that have ZWARNING set to zero
;   in one observation and non-zero in another.
;
; EXAMPLES:
;   Compare data from different pluggings of the same plate:
;   IDL> platecompare, [406,406,406], mjd=[51817,51869,51876]
;
;   Compare date from two plates, 360 and 362, which are actually the
;   same tile on the sky:
;   IDL> platecompare, [360,360,362], mjd=[51780,51816,51999]
;
;   Compare different reductions of the same data, in this case plate 401/51788:
;   IDL> platecompare, 401, mjd=51788, run2d=['v5_4_14','v5_4_30'], $
;    run1d=['v5_4_14','v5_4_30']
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
;   struct_append()
;
; REVISION HISTORY:
;   14-Aug-2001  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro platecompare, plate, mjd=mjd, run2d=run2d, run1d=run1d1, psfile=psfile

   if (n_params() LT 1) then begin
      print, 'Syntax - platecompare, plate, [mjd= ]'
      return
   endif

   charsize = 1.5

   if (keyword_set(run2d)) then begin
      if (keyword_set(run1d1)) then run1d = run1d1 $
       else run1d = run2d
      if (n_elements(run2d) NE n_elements(run1d)) then $
       message, 'Number of elements in RUN2D and RUN1D must agree!'
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
   if (NOT keyword_set(run2d)) then begin
      readspec, platevec, mjd=mjdvec, zans=zansall, plug=plugall
   endif else begin
      for itop=0, n_elements(run2d)-1 do begin
         readspec, platevec, mjd=mjdvec, zans=zans1, plug=plug1, $
          run2d=run2d[itop], run1d=run1d[itop]
         zansall = struct_append(zansall, zans1)
         plugall = struct_append(plugall, plug1)
         if (itop EQ 0) then begin
            all2d = replicate(run2d[itop], n_elements(zans1))
            all1d = replicate(run1d[itop], n_elements(zans1))
         endif else begin
            all2d = [all2d, replicate(run2d[itop], n_elements(zans1))]
            all1d = [all1d, replicate(run1d[itop], n_elements(zans1))]
         endelse
      endfor
   endelse

   if (NOT keyword_set(zansall)) then begin
      print, 'Data not found for plate ', platevec[0], ' MJD=', mjdvec[0]
      return
   endif

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
      if (keyword_set(run2d)) then begin
         thisr2d = all2d[gindx[gstart[igroup]:gstart[igroup]+gcount[igroup]-1]]
         thisr1d = all1d[gindx[gstart[igroup]:gstart[igroup]+gcount[igroup]-1]]
      endif

      vdiff = (max(thiszans.z) - min(thiszans.z)) * cspeed
      totalwarn = total(thiszans.zwarning NE 0)
      vbad = ((vdiff GT 250.) AND strmatch(thiszans[0].class,'GALAXY*')) $
          OR ((vdiff GT 500.) AND strmatch(thiszans[0].class,'QSO*')) $
          OR ((vdiff GT 50.) AND strmatch(thiszans[0].class,'STAR*')) $
          OR (totalwarn NE 0 AND totalwarn NE n_elements(thiszans))
      vbad = vbad AND ((thiszans[0].zwarning AND 1) EQ 0) ; Not SKY fiber

      xrange = [3600, 9300]
      textcolor = 'green'

      if (vbad) then begin
         !p.multi = [0, 1, gcount[igroup]]
         if (thiszans[0].sn_median[2] LT 1) then nsmooth = 7 $
          else if (thiszans[0].sn_median[2] LT 3) then nsmooth = 5 $
          else if (thiszans[0].sn_median[2] LT 5) then nsmooth = 3 $
          else nsmooth = 1

         for ii=0, gcount[igroup]-1 do begin
            primtarget = sdss_flagname('TARGET', thisplug[ii].primtarget, $
             /concat)
            sectarget = sdss_flagname('TTARGET', thisplug[ii].sectarget, $
             /concat)
            title = string(thiszans[ii].plate, thiszans[ii].mjd, $
             thiszans[ii].fiberid, $
             format='("Plate ", i4, "-", i5, " Fiber ", i3)')
            ; Replace undescores with double-underscores that will
            ; survive the conversion to TeX in DJS_PLOT
            if (keyword_set(run2d)) then $
             title = title + '  ' + repstr(thisr2d[ii],'_','__') $
              + ' ' + repstr(thisr1d[ii],'_','__')
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
            print, zstring
            if (keyword_set(run2d)) then read2d = thisr2d[ii] $
             else read2d = 0
            if (keyword_set(run2d)) then read1d = thisr1d[ii] $
             else read1d = 0
            readspec, thiszans[ii].plate, thiszans[ii].fiberid, $
             mjd=thiszans[ii].mjd, wave=wave, flux=objflux, synflux=synflux, $
             run2d=read2d, run1d=read1d
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
             xtitle=zstring, ytitle=ytitle, title=title, charsize=charsize
            djs_oplot, wave, synflux, color='blue'
            xpos = 0.95 * !x.crange[0] + 0.05 * !x.crange[1]
            ypos = 0.08 * !y.crange[0] + 0.92 * !y.crange[1]
            xyouts, xpos, ypos, primtarget+' '+sectarget, $
             color=djs_icolor(textcolor), charsize=charsize
         endfor
         if (NOT keyword_set(psfile)) then $
          cc = strupcase(get_kbrd(1))
      endif
   endfor

   if (keyword_set(psfile)) then dfpsclose

   return
end
;------------------------------------------------------------------------------
