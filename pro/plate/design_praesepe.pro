;
; Cluster center is around RA=130.0 deg, DEC=19.6 deg (J2000)

;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
function read_praesepe1, filename

   readcol, filename, $
    ra, dec, bmag, rmag, jmag, hmag, kmag, flag, $
    format='(D,D,F,F,F,F,F,L)'

   ; Fix some of the B magnitudes which are set to zero.
   ibad = where(bmag EQ 0, nbad)
   if (nbad GT 0) then begin
      bmag[ibad] = rmag[ibad] + 1.0 ; Approximate as B-R = 1.
   endif

   result = design_starstruct(n_elements(ra))
   result.ra = ra
   result.dec = dec
   result.mag = transpose( [[bmag], [bmag], [rmag], [rmag], [rmag]] )
   result.flag = flag

   return, result
end

;------------------------------------------------------------------------------
function read_praesepeall

   p3 = read_praesepe1('praesepe.j10')
   p2 = read_praesepe1('praesepe.match.j10')
   p1 = read_praesepe1('praesepe.j13')
   p3.priority = 3
   p2.priority = 2
   p1.priority = 1

   result = [p3, p2, p1]
   return, result
end

;------------------------------------------------------------------------------
pro design_praesepe

   epoch = 1998.
   racen = [130.0d, 129.8d, 130.2d]
   deccen = [19.6d,  19.6d, 19.6d]
   magmin = [ 6.0, 10.0, 13.0]
   magmax = [10.2, 13.2, 16.5]
   guidemag = [10.5, 12.5]

   ntile = n_elements(racen)

   ;----------
   ; Read stars from Eisenstein's lists which use 2MASS positions.
   ; Discard duplicates.

   stardata = read_praesepeall()
   junk = djs_angle_group(stardata.ra, stardata.dec, 5.0/3600., $
    gstart=gstart, gindx=gindx)
   stardata = stardata[gindx[gstart]]

   ;----------
   ; Read the Tycho stars

   tycdat = tyc_read(/small, epoch=epoch)
   adiff = djs_diff_angle(tycdat.radeg, tycdat.dedeg, racen[0], deccen[0])
   tycdat = tycdat[ where(adiff LT 1.5) ]

   ;----------
   ; For every Tycho star, find the match in the Eisenstein catalog.
   ; If there are objects not in Eisenstein, then add them.

   ntot = djs_angle_match(tycdat.radeg, tycdat.dedeg, $
    stardata.ra, stardata.dec, dtheta=5.0/3600., mindx=mindx, mdist=mdist)

   iadd = where(mindx EQ -1, nadd)
   if (nadd GT 0) then begin
      tycadd = design_starstruct(nadd)
      tycadd.ra = tycdat[iadd].radeg
      tycadd.dec = tycdat[iadd].dedeg
      ; Horrible approximations for magnitudes !!!???
      umag = 0.0 + tycdat[iadd].vmag + 1.0 * tycdat[iadd].bmv
      gmag = -0.149 + tycdat[iadd].vmag + 0.626 * tycdat[iadd].bmv
      rmag = gmag
      imag = gmag
      zmag = gmag
      tycadd.mag = transpose( [[umag],[gmag],[rmag],[imag],[zmag]] )
      tycadd.priority = 3
      stardata = [stardata, tycadd]
   endif

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
      addstar.tilenum = itile + 1
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

stop
;splot,stardata.ra,stardata.dec,ps=4,color='green'
;soplot,tycdat.radeg,tycdat.dedeg,ps=1,symsize=0.5

   design_multiplate, stardata, racen=racen, deccen=deccen

stop
end
;------------------------------------------------------------------------------
