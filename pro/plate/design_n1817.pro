; Design the special plate for the N1817 star cluster.

;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
function read_n1817

   p1 = design_read2mass('n1817.match')
   p2 = design_read2mass('n1817.j11')
   p1.priority = 3
   p2.priority = 2

   result = [p1, p2]
   return, result
end

;------------------------------------------------------------------------------
pro design_n1817

   epoch = 1998.
   racen = [78.0d, 77.8d, 78.2d, 78.4d]
   deccen = [16.7d,  16.7d, 16.7d, 16.7d]
   ; All magnitudes are in the range [8.2,18.1]
;   magmin = reverse([ 8.0, 10.8, 13.6])
;   magmax = reverse([11.0, 13.8, 16.6])
   magmin = reverse([ 8.0, 10.4, 12.8, 15.2])
   magmax = reverse([10.6, 13.0, 15.4, 18.2])
   guidemag = [10.5, 12.5]
   tilenums = [9231,9232,9233,9234]
   platenums = [801,1,2,3]
   matchdist = 2.0/3600. ; match distance in degrees

   ntile = n_elements(racen)

   ;----------
   ; Read stars from Eisenstein's lists which use 2MASS positions.
   ; Discard duplicates.

   stardata = read_n1817()
   junk = djs_angle_group(stardata.ra, stardata.dec, matchdist, $
    gstart=gstart, gindx=gindx)
   stardata = stardata[gindx[gstart]]

   ;----------
   ; Read the Tycho stars

;   tycdat = tyc_read(/small, epoch=epoch)
;   adiff = djs_diff_angle(tycdat.radeg, tycdat.dedeg, racen[0], deccen[0])
;   tycdat = tycdat[ where(adiff LT 1.5) ]

   ;----------
   ; For every Tycho star, find the match in the Eisenstein catalog.
   ; If there are objects not in Eisenstein, then add them.

;   junk = djs_angle_match(tycdat.radeg, tycdat.dedeg, $
;    stardata.ra, stardata.dec, dtheta=matchdist, mindx=mindx, mdist=mdist)
;
;   iadd = where(mindx EQ -1, nadd)
;   if (nadd GT 0) then begin
;      tycadd = design_starstruct(nadd)
;      tycadd.ra = tycdat[iadd].radeg
;      tycadd.dec = tycdat[iadd].dedeg
;      tycadd.mag = tyc_sdssmags(tycdat[iadd].bmv, tycdat[iadd].vmag)
;      tycadd.objtype = 'SERENDIPITY_MANUAL'
;      tycadd.priority = 3
;      stardata = [stardata, tycadd]
;   endif

   ;----------
   ; Select objects that may be good spectro-photo standards.
   ; Find stars with a match in the Tycho catalog with -0.1 < B-V < +0.1.

;   junk = djs_angle_match(tycdat.radeg, tycdat.dedeg, $
;    stardata.ra, stardata.dec, dtheta=matchdist, mindx=mindx, mdist=mdist)
;   indx1 = where(mindx NE -1)
;   indx2 = mindx[indx1]
;   iphoto = where(tycdat[indx1].bmv GT -0.1 AND tycdat[indx1].bmv LT 0.1)
;   if (iphoto[0] NE -1) then begin
;      stardata[indx2[iphoto]].objtype = 'SPECTROPHOTO_STD'
;      stardata[indx2[iphoto]].sectarget = 32L
;   endif

   ;----------
   ; Assign tile numbers based upon the magnitudes.
   ; Allow for overlapping magnitude slices, which assign some stars
   ; on several of the tiles.

   for itile=0, ntile-1 do begin
      iadd = where(stardata.mag[2] GE magmin[itile] $
               AND stardata.mag[2] LE magmax[itile], nadd)
      if (nadd EQ 0) then $
       message, 'No objects available for this tile'
      addstar = stardata[iadd]
      addstar.tilenum = tilenums[itile]
      addstar.holetype = 'OBJECT'
      if (itile EQ 0) then newdata = addstar $
       else newdata = [newdata, addstar]
   endfor

   ;----------
   ; Now include all stars within an appropriate magnitude range
   ; as possible guide stars.  These will be duplicate entries
   ; but with HOLETYPE = 'GUIDE'

   iadd = where(stardata.mag[2] GE guidemag[0] $
            AND stardata.mag[2] LE guidemag[1], nadd)
   if (nadd EQ 0) then $
    message, 'No guide stars available'
   addstar = stardata[iadd]
   addstar.tilenum = 0
   addstar.holetype = 'GUIDE'
   newdata = [newdata, addstar]

   stardata = newdata

;splot,stardata.ra,stardata.dec,ps=4,color='green'
;soplot,tycdat.radeg,tycdat.dedeg,ps=1,symsize=0.5

   design_multiplate, stardata, racen=racen, deccen=deccen, $
    tilenums=tilenumes, platenums=platenums, /addfund

end
;------------------------------------------------------------------------------
