;+
; NAME:
;   design_multiplate
;
; PURPOSE:
;   Routine to design a plate with several pointings.
;
; CALLING SEQUENCE:
;   design_multiplate, stardata, [ tilenums=, platenums=, racen=, deccen=, $
;    guidetiles= ]
;
; INPUTS:
;   stardata   - Structure with data for each star; must contain the
;                fields RA, DEC, MAG[5], PRIORITY, HOLETYPE, TILENUM.
;                HOLETYPE can be either 'OBJECT' or 'GUIDE'; objects with
;                other values are ignored.  Other elements in the structure
;                will be copied into the output plug-map files.
;                Stars with TILENUM=0 can only be used as guide stars.
;   racen      - RA center for each tile
;   deccen     - DEC center for each tile
;
; OPTIONAL INPUTS:
;   tilenums   - Array of tile numbers; default to the unique list of
;                tile numbers in STARDATA.TILENUM, but exclude the number zero.
;   platenums  - Array of plate numbers; default to the same as TILENUMS.
;   guidetiles - Tile number for each of the 11 guide fibers.  There exist
;                default values for the cases 1,2,3 or 4 tiles.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   There is no attempt to move plate centers such that a guide fiber
;   can be re-used in another pointing.
;
;   Objects are assigned according to their priority.  If none are specified
;   with STARDATA.PRIORITY, then random priorities between 1 and 100 are
;   assigned.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   cpbackup
;   djs_diff_angle()
;   djs_laxisgen()
;   plate_rotate
;   yanny_free
;   yanny_par()
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;   design_append()
;
; REVISION HISTORY:
;   21-Nov-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
; Search for conflicts between an existing list of drill holes and one
; more potential object.  Return the existing list with the new object
; appended if there was no conflict.

function design_append, newplug, oneplug

   platescale = 217.7358 ; mm/degree

   if (NOT keyword_set(newplug)) then return, oneplug

   ; Discard objects within 55 arcsec of existing objects.
   ; Do this based upon XFOCAL,YFOCAL positions.
   ; The closest two fibers can be is PLATESCALE * 55/3600(deg) = 3.32652 mm

   r2 = (newplug.xfocal - oneplug.xfocal)^2 $
        + (newplug.yfocal - oneplug.yfocal)^2
   mindist = min(sqrt(r2))
   if (mindist LT platescale*55./3600.) then return, newplug

   return, [newplug, oneplug]
end
;------------------------------------------------------------------------------
pro design_multiplate, stardata, tilenums=tilenums, platenums=platenums, $
 racen=racen, deccen=deccen, guidetiles=guidetiles

   if (NOT keyword_set(tilenums)) then begin
      tilenums = stardata.tilenum
      tilenums = tilenums[ uniq(tilenums, sort(tilenums)) ]
      ; Exclude any tile numbers of zero
      tilenums = tilenums[ where(tilenums NE 0) ]
   endif

   if (NOT keyword_set(platenums)) then platenums = tilenums

   ntile = n_elements(tilenums)
   if (n_elements(platenums) NE ntile) then $
    message, 'Number of elements in PLATENUMS must agree w/ number of tiles'
   if (n_elements(racen) NE ntile OR n_elements(deccen) NE ntile) then $
    message, 'Number of elements in RACEN,DECCEN must agree w/ number of tiles'

   if (NOT keyword_set(guidetiles)) then begin
      ; First assign tile numbers like they are 1-indexed
      if (ntile EQ 1) then begin
         guidetiles = [1,1,1,1,1,1,1,1,1,1,1]
      endif else if (ntile EQ 2) then begin
         ; For 1st pointing use guide fibers 2,4,6,8,9,10,11
         ; For 2nd pointing use guide fibers 1,3,5,7
         guidetiles = [2,1,2,1,2,1,2,1,1,1,1]
      endif else if (ntile EQ 3) then begin
         ; For 1st pointing use guide fibers 3,5,8,10,11
         ; For 2nd pointing use guide fibers 1,4,6
         ; For 3rd pointing use guide fibers 2,7,9
         guidetiles = [2,3,1,2,1,2,3,1,3,1,1]
      endif else if (ntile EQ 4) then begin
         ; For 1st pointing use guide fibers 3,5,8,10,11
         ; For 2nd pointing use guide fibers 1,6
         ; For 3rd pointing use guide fibers 4,7
         ; For 4th pointing use guide fibers 2,9
         guidetiles = [2,4,1,3,1,2,3,1,4,1,1]
      endif else begin
         message, 'GUIDETILES must be specified if using more than 4 tiles'
      endelse
      ; Now re-assigne the guide tile numbers to the actual tile numbers
      guidetiles = tilenums[guidetiles-1]
   endif
   if (N_elements(guidetiles) NE 11) then $
    message, 'GUIDETILES must be an 11-element array'

   if (NOT keyword_set(racen) OR NOT keyword_set(deccen)) then $
    message, 'RACEN,DECCEN must be specified'

   plugmaptfile = 'plPlugMapT-' + string(tilenums,format='(i4.4)') + '.par'
   plugmappfile = 'plPlugMapP-' + string(platenums,format='(i4.4)') + '.par'

   fakemag = 25.0 ; Magnitudes for all fake objects

   ;----------
   ; Read a template plugmap structure

   yanny_read, 'plPlugMapT-XXXX.par', pp, $
    hdr=plughdr, enums=plugenum, structs=plugstruct
   blankplug = *pp[0]
   yanny_free, pp
   struct_assign, {junk:0}, blankplug

   ;----------
   ; Set up info for guide fibers.
   ;
   ; The following info is from the "plate" product in the
   ; file "$PLATE_DIR/test/plParam.par".
   ;   XREACH,YREACH = Center of the fiber reach [mm]
   ;   RREACH = Radius of the fiber reach [mm]
   ;   XPREFER,YREACH = Preferred position for the fiber [mm]
   ; Note that the plate scale is approx 217.7358 mm/degree.
   ; Moving +RA is +XFOCAL, +DEC is +YFOCAL.

   gfiber = create_struct( $
    'xreach'   , 0.0, $
    'yreach'   , 0.0, $
    'rreach'   , 0.0, $
    'xprefer'  , 0.d, $
    'yprefer'  , 0.d )
   nguide = 11
   gfiber = replicate(gfiber, nguide)

   platescale = 217.7358 ; mm/degree
   guideparam = [[  1,  199.0,  -131.0,  165.0,  199.0,  -131.0 ], $
                 [  2,   93.0,  -263.0,  165.0,   93.0,  -263.0 ], $
                 [  3, -121.0,  -263.0,  165.0, -121.0,  -263.0 ], $
                 [  4, -227.0,  -131.0,  165.0, -227.0,  -131.0 ], $
                 [  5, -199.0,   131.0,  165.0, -199.0,   131.0 ], $
                 [  6,  -93.0,   263.0,  165.0,  -93.0,   263.0 ], $
                 [  7,  121.0,   263.0,  165.0,  121.0,   263.0 ], $
                 [  8,  227.0,   131.0,  165.0,  227.0,   131.0 ], $
                 [  9,   14.0,   131.0,  139.5,   14.0,    65.0 ], $
                 [ 10,  -14.0,  -131.0,  165.0,  -14.0,   -65.0 ], $
                 [ 11,   93.0,  -131.0,  139.5,   93.0,  -131.0 ] ]
   gfiber.xreach = transpose(guideparam[1,*])
   gfiber.yreach = transpose(guideparam[2,*])
   gfiber.rreach = transpose(guideparam[3,*])
   gfiber.xprefer = transpose(guideparam[4,*])
   gfiber.yprefer = transpose(guideparam[5,*])

   ;---------------------------------------------------------------------------
   ; MAIN LOOP OVER POINTING NUMBER
   ;---------------------------------------------------------------------------

   for itile=0, ntile-1 do begin

      allplug = 0

      thistilenum = tilenums[itile]
      thisracen = racen[itile]
      thisdeccen = deccen[itile]

      ;----------
      ; For this pointing, find the approximate RA,DEC coordinates
      ; for the optimal guide fiber positioning.

      dphi = gfiber.xprefer / platescale
      dtheta = gfiber.yprefer / platescale
      plate_rotate, 0.0, 0.0, thisracen, thisdeccen, $
       dphi, dtheta, guidera, guidedec

      ;----------
      ; Loop through each guide fiber

      for iguide=0, nguide-1 do begin
         if (guidetiles[iguide] EQ thistilenum) then begin

            ;----------
            ; Case where the guide fiber is on this pointing.
            ; Assign the nearest available guide fiber

print, 'Assigning real guide fiber number ', iguide+1
            indx = where(strtrim(stardata.holetype,2) EQ 'GUIDE')
            adiff = djs_diff_angle(guidera[iguide], guidedec[iguide], $
             stardata[indx].ra, stardata[indx].dec)
            junk = min(adiff, imin)
            addplug = blankplug
            struct_assign, stardata[indx[imin]], addplug
            addplug.throughput = 2L^31-1 ; Maximum value

         endif else begin

            ;----------
            ; Case where the guide fiber will not be used on this pointing.
            ; Assign a bunch of possible guide fibers near the optimal position.
            ; This will be a grid of 5x5 positions separated by 3 arcmin.

            ngrid = [5,5]
            dphi = (djs_laxisgen(ngrid, iaxis=0) - (ngrid[0]-1)/2) * 3./60.
            dtheta = (djs_laxisgen(ngrid, iaxis=1) - (ngrid[1]-1)/2) * 3./60.
            plate_rotate, 0.0, 0.0, guidera[iguide], guidedec[iguide], $
             dphi[*], dtheta[*], gridra, griddec

            addplug = replicate(blankplug, ngrid[0]*ngrid[1])
            addplug.ra = gridra
            addplug.dec = griddec
            addplug.throughput = 0L
            addplug.mag[*] = fakemag

         endelse

         addplug.holetype = 'GUIDE'
         addplug.objtype = 'NA'
;         addplug.fiberid = iguide + 1 ; ???
         addplug.sectarget = 64L

         if (NOT keyword_set(allplug)) then allplug = addplug $
          else allplug = [allplug, addplug]

      endfor ; End loop over guide fibers

      ;----------
      ; For this pointing, add all objects.

      indx = where(stardata.tilenum EQ thistilenum $
       AND strtrim(stardata.holetype,2) EQ 'OBJECT', nadd)
      if (nadd EQ 0) then $
       message, 'No objects found for this pointing'
      addplug = replicate(blankplug, nadd)
      struct_assign, stardata[indx], addplug

      addplug.holetype = 'OBJECT'
      addplug.objtype = 'SERENDIPITY_MANUAL'
      if ((where(tag_names(stardata) EQ 'PRIORITY'))[0] NE -1) then $
       addplug.throughput = (stardata[indx].priority > 1L) < (2L^31-2) $
      else $
       addplug.throughput = long(randomu(24680, nadd) * 100) + 1
      addplug.sectarget = 2L^24

      allplug = [allplug, addplug]

      ;----------
      ; For this pointing, add 800 fake skies randomly distributed.

      nadd = 800
      dphi = 3.0 * randomu(2345, nadd) - 1.5
      dtheta = 3.0 * randomu(3456, nadd) - 1.5
      plate_rotate, 0.0, 0.0, thisracen, thisdeccen, $
       dphi, dtheta, gridra, griddec

      addplug = replicate(blankplug, nadd)
      addplug.ra = gridra
      addplug.dec = griddec

      addplug.holetype = 'COHERENT_SKY'
      addplug.objtype = 'NA'
      addplug.sectarget = 16L
      addplug.throughput = 0L
      addplug.mag[*] = fakemag

      allplug = [allplug, addplug]

      ;----------
      ; For this pointing, add a grid of fake spectro-photo and reddening stars.
      ; This will be a grid of 13x13 positions separated by 13 arcmin.
      ; Alternate between calling them SPECTROPHOTO_STD or REDDEN_STD.

      ngrid = [13,13]
      dphi = (djs_laxisgen(ngrid, iaxis=0) - (ngrid[0]-1)/2) * 13./60.
      dtheta = (djs_laxisgen(ngrid, iaxis=1) - (ngrid[1]-1)/2) * 13./60.
      plate_rotate, 0.0, 0.0, thisracen, thisdeccen, $
       dphi[*], dtheta[*], gridra, griddec

      addplug = replicate(blankplug, ngrid[0]*ngrid[1])
      addplug.ra = gridra
      addplug.dec = griddec

      addplug.holetype = 'OBJECT'
      ieven = where(lindgen(ngrid[0]*ngrid[1]) MOD 2 EQ 0)
      iodd = where(lindgen(ngrid[0]*ngrid[1]) MOD 2 EQ 1)
      addplug[ieven].objtype = 'SPECTROPHOTO_STD'
      addplug[ieven].sectarget = 32L
      addplug[iodd].objtype = 'REDDEN_STD'
      addplug[iodd].sectarget = 2L
      addplug.throughput = 0L
      addplug.mag[*] = fakemag

      allplug = [allplug, addplug]

      ;----------
      ; Remove objects more than 1.49 deg from the center.
      ; Also remove objects within 68 arcsec of the center hole.
      ; (This is a QUALITY hole added somewhere in the PLATE product.)

      adiff = djs_diff_angle(thisracen, thisdeccen, allplug.ra, allplug.dec)
      allplug = allplug[ where(adiff LT 1.49 AND adiff GT 68./3600.) ]

      ;----------
      ; Resolve conflicts.
      ; First re-sort everything from highest priority to lowest, so
      ; that we can always keep the 1st object in the list in the event
      ; of a conflict.

      isort = reverse(sort(allplug.throughput))
      allplug = allplug[isort]

      ; Loop through all objects, discarding close pairs (within 55 arcsec)
      ; Always keep the 1st of any close group of objects.
      ; This is a stupidly slow N^2 implementation...

      iobj = 0L
      while (iobj LT n_elements(allplug)-1) do begin
         adiff = djs_diff_angle(allplug[iobj].ra, allplug[iobj].dec, $
          allplug.ra, allplug.dec)
         igood = where(adiff GT 55.0/3600.0 $
          OR lindgen(n_elements(allplug)) EQ iobj, ngood)
         allplug = allplug[igood]
         iobj = iobj + 1L
      endwhile

      ;----------
      ; Write the plPlugMapT file for this tile...

      outhdr = ['completeTileVersion   v1_0', $
                'tileId ' + string(thistilenum), $
                'raCen ' + string(thisracen), $
                'decCen ' + string(thisdeccen) ]
      for jtile=0, ntile-1 do begin
         outhdr = [outhdr, $
          'raCen'+strtrim(string(jtile),2)+ '  '+string(racen[jtile]) ]
         outhdr = [outhdr, $
          'decCen'+strtrim(string(jtile),2)+ '  '+string(deccen[jtile]) ]
      endfor

      yanny_write, plugmaptfile[itile], ptr_new(allplug), hdr=outhdr, $
       enums=plugenum, structs=plugstruct

   endfor ; End loop over pointing number

   ;---------------------------------------------------------------------------
   ; RUN "makePlates" IN THE SDSS "PLATE" PRODUCT.
   ; The required inputs are the plPlugMapT-$TILE.par files,
   ; plus plPlan.par, plObs.par, plParam.par.
   ;---------------------------------------------------------------------------

; ???
print, 'Now run "makePlates" in the "plate" product'
stop

   ;---------------------------------------------------------------------------
   ; RE-COMBINE THE PLUGMAP FILES
   ;---------------------------------------------------------------------------

   ;----------
   ; Read the plugMapP files and track which plate number each object is from
   ; by storing the plate number in the FIBERID field.

   hdrarr = replicate(ptr_new(), ntile)
   for itile=0, ntile-1 do begin
      yanny_read, plugmappfile[itile], pp, $
       hdr=plughdr, enums=plugenum, structs=plugstruct
      hdrarr[itile] = ptr_new(plughdr)
      thisplate = yanny_par(plughdr, 'plateId')
      (*pp[0]).fiberid = thisplate ; Store the plate number in FIBERID
      if (itile EQ 0) then allplug = *pp[0] $
       else allplug = [allplug, *pp[0]]
      if (itile EQ 0) then platenums = thisplate $
       else platenums = [platenums, thisplate]
      yanny_free, pp
   endfor

   ;----------
   ; Select the 11 non-fake guide stars

   iguide = where(allplug.holetype EQ 'GUIDE' AND allplug.throughput GT 0)
   if (n_elements(iguide) NE 11) then $
    message, 'The number of guide fibers is wrong.'
   newplug = allplug[iguide]

   ;----------
   ; Select all real objects and then fake skies if we run out of real objects.
   ; Sort by priority (in THROUGHPUT field), but then add tiny random numbers
   ; in order to randomize between objects with the same priority.

   iobj = where(allplug.holetype EQ 'OBJECT', nobj)
   isort = reverse(sort(allplug[iobj].throughput + randomu(2468,nobj)))
   iobj = iobj[isort]
   nadded = 0
   for ii=0, nobj-1 do begin
      if (nadded LT 640) then begin
         nbefore = n_elements(newplug)
         newplug = design_append(newplug, allplug[iobj[ii]])
         if (n_elements(newplug) GT nbefore) then nadded = nadded + 1
       endif
   endfor
   if (nadded LT 640) then $
    message, 'Fewer than 640 objects'

   ;---------------------------------------------------------------------------
   ; WRITE THE MODIFIED PLUGMAPP FILES.
   ; This combines objects from the different tiles.
   ; This overwrites files that already exist.
   ;---------------------------------------------------------------------------

   for itile=0, n_elements(platenums)-1 do begin

      thisplate = platenums[itile]
      modplug = newplug

      ;----------
      ; Keep only the guide fibers on this tile. ???


      ;----------
      ; Objects that are actually on other tiles are renamed to sky fibers
      ; on this plate.

      indx = where(modplug.holetype EQ 'OBJECT' $
       AND modplug.fiberid NE thisplate)
      if (indx[0] NE -1) then begin
         modplug[indx].holetype = 'COHERENT_SKY'
         modplug[indx].objtype = 'NA'
         modplug[indx].primtarget = 0L
         modplug[indx].sectarget = 16L
      endif else begin
         message, 'No objects to use as guide fibers for plate '+string(thisplate)
      endelse

      ;----------
      ; Add the one quality hole which we already have, which is the
      ; center hole

      iqual = where(allplug.holetype EQ 'QUALITY' $
       AND allplug.fiberid EQ thisplate, nqual)
      if (nqual NE 1) then $
      message, 'We expect 1 quality hole already (the center hole)'
      modplug = [modplug, allplug[iqual]]

      ; Must set the following to avoid PLATE from crashing!!
      modplug.spectrographid = -9999
      modplug.fiberid = -9999
      modplug.throughput = -9999

      cpbackup, plugmappfile[itile]

      yanny_write, plugmappfile[itile], ptr_new(modplug), $
       hdr=*(hdrarr[itile]), enums=plugenum, structs=plugstruct
   endfor

   ;----------
; Also need to keep QUALITY and ALIGNMENT holes !!!???

   ;---------------------------------------------------------------------------
   ; CREATE THE DRILL FILES FROM THE 1ST PLATE.
   ; (Drill files from any of the plates would be identical.)
   ; Run the code fiberPlates, makeFanuc, makeDrillPos, use_cs3.
   ; The fiberPlates code selects the guide stars and sky fibers from
   ; those available, and renames COHERENT_SKY/NA objects to OBJECT/SKY.
   ;---------------------------------------------------------------------------

; ???
print, 'Now run "fiberPlates" in the "plate" product'
stop

   return
end

