;+
; NAME:
;   design_multiplate
;
; PURPOSE:
;   Routine to design a plate with several pointings.
;
; CALLING SEQUENCE:
;   design_multiplate, stardata, [ tilenums=, platenums=, racen=, deccen=, $
;    guidetiles=, apotemperature=apotemperature, /addfund, /norename ]
;
; INPUTS:
;   stardata   - Structure with data for each star; must contain the
;                fields RA, DEC, MAG[5], PRIORITY, HOLETYPE, TILENUM.
;                HOLETYPE can be either 'OBJECT' or 'GUIDE'; objects with
;                other values are ignored.  Other elements in the structure
;                will be copied into the output plug-map files.
;                If OBJTYPE is not passed in this structure, then it is set to
;                'SERENDIPITY_MANUAL' for all HOLETYPE='OBJECT'.
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
;   apotemperature - Design temperature for APO; default to 5 deg C.
;   addfund    - If set, then add one more tile+plate with 9 or 10 holes
;                in a line of declination from guide fiber #11 for the
;                nearest of the 5 SDSS fundamental standard stars.
;                We try to place 10 of these holes, spaced by 1 arcmin, but
;                one can be knocked out by the guide-fiber alignment hole.
;                This tile is given the TILEID=max(TILEID)+1, PLATEID=0.
;   norename   - The default is to rename the output plugMapP files to
;                plate names like 800,800B,800C,... if the first plate
;                number is 800.  Set this keyword to keep the names as
;                they are passed in PLATENUMS.
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
;   assigned.  We reserve priorities of 0 and [2^31-2,2^31] for internal
;   purposes.
;
; EXAMPLES:
;
; BUGS:
;   These SDSS primary standards coordinates were not quite correct
;   in the design of the initial set of special plates.  See the comments
;   in the code.
;
; PROCEDURES CALLED:
;   concat_dir()
;   cpbackup
;   current_mjd()
;   djs_diff_angle()
;   djs_laxisgen()
;   plate_rotate
;   yanny_free
;   yanny_par()
;   yanny_read
;   yanny_write
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

function design_append, newplug, newplatearr, oneplug, oneplate

   platescale = 217.7358 ; mm/degree

   ; If this is the 1st object in the list, then we can always keep it
   if (NOT keyword_set(newplug)) then begin
      newplatearr = oneplate
      return, oneplug
   endif

   ; Discard objects within 55 arcsec of existing objects.
   ; Do this based upon XFOCAL,YFOCAL positions.
   ; The closest two fibers can be is PLATESCALE * 55/3600(deg) = 3.32652 mm

   r2 = (newplug.xfocal - oneplug.xfocal)^2 $
        + (newplug.yfocal - oneplug.yfocal)^2
   mindist = min(sqrt(r2))
   if (mindist LT platescale*55./3600.) then return, newplug

   ; Discard objects within 7.0 mm of existing guide fibers

   iguide = where(strtrim(newplug.holetype) EQ 'GUIDE')
   if (iguide[0] NE -1) then begin
      r2 = (newplug[iguide].xfocal - oneplug.xfocal)^2 $
           + (newplug[iguide].yfocal - oneplug.yfocal)^2
      mindist = min(sqrt(r2))
      if (mindist LT 7.0) then return, newplug
   endif

   newplatearr = [newplatearr, oneplate]
   return, [newplug, oneplug]
end
;------------------------------------------------------------------------------
pro design_multiplate, stardata, tilenums=tilenums, platenums=platenums, $
 racen=racen, deccen=deccen, guidetiles=guidetiles, $
 apotemperature=apotemperature, addfund=addfund, norename=norename

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
      endif else if (ntile EQ 5) then begin
         ; For 1st pointing use guide fibers 3,5,8,10,11
         ; For 2nd pointing use guide fibers 1,6
         ; For 3rd pointing use guide fibers 4,7
         ; For 4th pointing use guide fibers 2
         ; For 5th pointing use guide fibers 9
         guidetiles = [2,4,1,3,1,2,3,1,5,1,1]
      endif else begin
         message, 'GUIDETILES must be specified if using more than 5 tiles'
      endelse
      ; Now re-assigne the guide tile numbers to the actual tile numbers
      guidetiles = tilenums[guidetiles-1]
   endif
   if (N_elements(guidetiles) NE 11) then $
    message, 'GUIDETILES must be an 11-element array'

   if (NOT keyword_set(racen) OR NOT keyword_set(deccen)) then $
    message, 'RACEN,DECCEN must be specified'

   if (N_elements(apotemperature) EQ 0) then $
    apotemperature = 5.0

   plugmaptfile = 'plPlugMapT-' + string(tilenums,format='(i4.4)') + '.par'
   plugmappfile = 'plPlugMapP-' + string(platenums,format='(i4.4)') + '.par'

   fakemag = 25.0 ; Magnitudes for all fake objects
   maxpriority = 2L^31 - 1 ; Maximum value; this is the value for GUIDE stars
   fundpriority = maxpriority - 2 ; Priority for fundamental standards
   paramdir = concat_dir(getenv('IDLSPEC2D_DIR'), 'examples')

   if (keyword_set(addfund)) then begin
      if ((where(tilenums EQ 0))[0] NE -1) then $
       message, 'PLATEID of 0 is reserved for the fundamental star tile'
      fundplatenum = 0
      fundtilenum = max(tilenums) + 1
   endif

   ;----------
   ; Set up the data for the fundamental standards, and apply the
   ; proper motion corrections from epoch 2000 to the current epoch.
   ; The BD magnitudes are those in the file 'metaFC_bd17isonlyfund.fit'
   ; checked into the "mtstds" product on 16 June 2000.
   ; The HD magnitudes are from Fukugita et al (1996) AJ 111, 1748.

   pfund = $
   { name:         '',  $
     ra:         0.0d,  $
     dec:        0.0d,  $
     pm_ra:      0.0 ,  $
     pm_dec:     0.0 ,  $
     vmag:       0.0 ,  $
     mag:        fltarr(5) }
   funddat = replicate(pfund, 5)
   funddat.name  = ['HD_19445','BD+21o607','HD_84937','BD+26o2606','BD+17o4708']
; These first set of coordinates, which are not quite correct, were used
; in the design of the initial set of special plates.
;   funddat.ra    =[  47.10663,   63.64800, 147.23375,  222.25961,   332.87167 ]
;   funddat.dec   =[  26.33094,   22.35119,  13.74426,   25.70750,    18.09139 ]
;   funddat.pm_ra =[    -0.210,      0.425,     0.373,     -0.009,       0.512 ]
;   funddat.pm_dec=[    -0.830,       9.22,    -0.774,     -0.346,       0.060 ]
   funddat.ra    =[  47.10579,   63.64573, 147.23510,  222.25980,   332.88109 ]
   funddat.dec   =[  26.33096,   22.35116,  13.74427,   25.70250,    18.09277 ]
   funddat.pm_ra =[    -0.210,      0.425,     0.373,     -0.009,       0.512 ]
   funddat.pm_dec=[    -0.830,       9.22,    -0.774,     -0.346,       0.060 ]
   funddat.vmag  = [      8.05,       10.0,      8.28,       9.72,        9.47]
   funddat[0].mag = [ 9.08 , 8.23 , 7.92 , 7.82,  7.79 ] ; HD_19445
   funddat[1].mag = [10.289, 9.395, 9.114, 9.025, 9.017] ; BD+21 (V=10.0)
   funddat[2].mag = [ 9.32 , 8.46 , 8.23 , 8.16,  8.16 ] ; HD_84937
   funddat[3].mag = [10.761, 9.891, 9.604, 9.503, 9.486] ; BD+26 (V=9.72)
   funddat[4].mag = [10.560, 9.640, 9.350, 9.250, 9.230] ; BD+17 (V=9.47)

   thismjd = current_mjd()
   funddat.ra = funddat.ra $
    + (funddat.pm_ra / 3600.) * (thismjd - 51544)/365.24
   funddat.dec = funddat.dec $
    + (funddat.pm_dec / 3600.) * (thismjd - 51544)/365.24

   ;----------
   ; Read a template plugmap structure

   yanny_read, filepath('plPlugMapT-XXXX.par', root_dir=paramdir), pp, $
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
      print
      print, 'Working on TILE=', thistilenum, ' PLATE=', platenums[itile]

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

            print, 'Assigning real guide fiber number ', iguide+1, $
             ' on tile ', thistilenum
            indx = where(strtrim(stardata.holetype,2) EQ 'GUIDE')
            adiff = djs_diff_angle(guidera[iguide], guidedec[iguide], $
             stardata[indx].ra, stardata[indx].dec)
            junk = min(adiff, imin)
            addplug = blankplug
            struct_assign, stardata[indx[imin]], addplug
            addplug.throughput = maxpriority ; Maximum value

            ; Save the coordinates of the real guide star #11
            ; and the corresponding plate center.
            if (iguide+1 EQ 11) then begin
               gra1 = thisracen
               gdec1 = thisdeccen
               gra2 = stardata[indx[imin]].ra
               gdec2 = stardata[indx[imin]].dec
            endif

         endif else begin

            ;----------
            ; Case where the guide fiber will not be used on this pointing.
            ; Assign a bunch of possible guide fibers near the optimal position.
            ; This will be a grid of 7x7 positions separated by 3 arcmin.

            ngrid = [7,7]
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
         addplug.sectarget = 64L

         ; The "fiberPlates" code will overwrite the FIBERID.  I don't
         ; want it to change them, like decide that guide fiber #11 is
         ; really better as #1.  So store the FIBERID that I want in
         ; PRIMTARGET, then re-assign it later.
         addplug.fiberid = iguide + 1
         addplug.primtarget = iguide + 1

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
      if ((where(tag_names(stardata) EQ 'OBJTYPE'))[0] NE -1) then $
       addplug.objtype = stardata[indx].objtype
      if ((where(tag_names(stardata) EQ 'PRIORITY'))[0] NE -1) then $
       addplug.throughput = (stardata[indx].priority > 1L) < (maxpriority-3) $
      else $
       addplug.throughput = long(randomu(24680, nadd) * 100) + 1
;      addplug.sectarget = 2L^24 ; This would be the serendipity flag

      ; Call any objects serendipity if not called anything else...
      for i=0, nadd-1 do begin
         if (strtrim(addplug[i].objtype) EQ '') then begin
            addplug[i].objtype = 'SERENDIPITY_MANUAL'
            addplug[i].primtarget = 2L^24
         endif
      endfor

      allplug = [allplug, addplug]

      ;----------
      ; For this pointing, add 1500 fake skies randomly distributed.

      nadd = 1500
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
      ; --> Actually, don't call anything a reddening standard, since
      ;     we don't need them for the code to run.

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
; Don't call anything a reddening standard, since we don't need them.
addplug.objtype = 'SPECTROPHOTO_STD'
addplug.sectarget = 32L
;      addplug[ieven].objtype = 'SPECTROPHOTO_STD'
;      addplug[ieven].sectarget = 32L
;      addplug[iodd].objtype = 'REDDEN_STD'
;      addplug[iodd].sectarget = 2L
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

      print, 'Before resolving conflicts: Number(HOLETYPE=OBJECT) = ', $
       n_elements(where(strtrim(allplug.holetype) EQ 'OBJECT'))
      iobj = 0L
      while (iobj LT n_elements(allplug)-1) do begin
         adiff = djs_diff_angle(allplug[iobj].ra, allplug[iobj].dec, $
          allplug.ra, allplug.dec)
         igood = where(adiff GT 55.0/3600.0 $
          OR lindgen(n_elements(allplug)) EQ iobj, ngood)
         allplug = allplug[igood]
         iobj = iobj + 1L
      endwhile
      print, 'After resolving conflicts:  Number(HOLETYPE=OBJECT) = ', $
       n_elements(where(strtrim(allplug.holetype) EQ 'OBJECT'))

      ;----------
      ; Trim to only 600 real objects at most, or "fiberPlates" will crash
      ; (excluding spectro-photo standards).  That routine expects 600
      ; objects at most + 32 skies + 8 spectro-photos.

      iobj = where(strtrim(allplug.holetype) EQ 'OBJECT' $
       AND strtrim(allplug.objtype) NE 'SPECTROPHOTO_STD', nobj)
      if (nobj GT 600) then begin
         qgood = lindgen(n_elements(allplug)) + 1B
         qgood[iobj[600:nobj-1]] = 0 ; Mark these last ones as bad
         allplug = allplug[ where(qgood) ]
      endif
      print, 'After trimming:             Number(HOLETYPE=OBJECT) = ', $
       n_elements(where(strtrim(allplug.holetype) EQ 'OBJECT'))

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
   ; ADD HOLES FOR FUNDAMENTAL STANDARD
   ;
   ; Add a line of holes north/south of guide fiber #11 with very high priority.
   ; Append another tile+plate for this.
   ;---------------------------------------------------------------------------

   if (keyword_set(addfund)) then begin
      ntile = ntile + 1
      tilenums = [tilenums, fundtilenum]
      platenums = [platenums, fundplatenum]
      plugmaptfile = [plugmaptfile, $
       'plPlugMapT-' + string(fundtilenum,format='(i4.4)') + '.par']
      plugmappfile = [plugmappfile, $
       'plPlugMapP-' + string(fundplatenum,format='(i4.4)') + '.par']

      ; Find the nearest fundamental standard
      adiff = djs_diff_angle(funddat.ra, funddat.dec, gra1, gdec1)
      mindiff = min(adiff, ifund)
      print, 'Selecting the fundamental star ', funddat[ifund].name, $
       ' dAngle=', mindiff, ' degrees'

      ; Pretend that we had a pointing with guide#11 (GRA2,GDEC2) at the center
      ; position and an object at GRA1,GDEC1.  Move to a pointing with the
      ; fundamental star (FUNDDAT[IFUND].RA,DEC) in the center, and solve for
      ; the "object position" -- which will actually be the plate center.

      print, 'Put fundamental star fibers near guide #11 at ', gra2, gdec2
;      plate_rotate, gra2, gdec2, funddat[ifund].ra, funddat[ifund].dec, $
;       gra1, gdec2, fundracen, funddeccen
      plate_newcenter, gra1, gdec1, fundracen, funddeccen, $
       gra2, gdec2, funddat[ifund].ra, funddat[ifund].dec

      racen = [racen, fundracen]
      deccen = [deccen, funddeccen]

      ;----------
      ; Copy the last plug-map data, calling everything priority zero.
      ; Rotate all those coordinates from the last plate center
      ; to the fundamental-star plate center.

      allplug.throughput = 0
      allplug.mag[*] = fakemag
      plate_rotate, thisracen, thisdeccen, fundracen[0], funddeccen[0], $
       allplug.ra, allplug.dec, tmpra, tmpdec
      allplug.ra = tmpra
      allplug.dec = tmpdec

      ;----------
      ; Add a line of holes in declination offset from guide fiber #11.

      decoffset =  [-5,-4,-3,-2,-1,1,2,3,4,5]/60.0
      nadd = n_elements(decoffset)
      addplug = replicate(blankplug, nadd)
      addplug.ra = funddat[ifund].ra
      addplug.dec = funddat[ifund].dec + decoffset
      addplug.mag = funddat[ifund].mag
      addplug.throughput = fundpriority ; Very high priority for these holes
      addplug.holetype = 'OBJECT'
      addplug.objtype = 'SERENDIPITY_MANUAL' ; This will be changed later
                                             ; to SPECTROPHOTO_STD
      allplug = [allplug, addplug]

      ;----------
      ; Write the plPlugMapT file for this tile...

      outhdr = ['completeTileVersion   v1_0', $
                'tileId ' + string(fundtilenum), $
                'raCen ' + string(fundracen), $
                'decCen ' + string(funddeccen) ]
      for jtile=0, ntile-1 do begin
         outhdr = [outhdr, $
          'raCen'+strtrim(string(jtile),2)+ '  '+string(racen[jtile]) ]
         outhdr = [outhdr, $
          'decCen'+strtrim(string(jtile),2)+ '  '+string(deccen[jtile]) ]
      endfor

      yanny_write, plugmaptfile[itile], ptr_new(allplug), hdr=outhdr, $
       enums=plugenum, structs=plugstruct
   endif

   ;---------------------------------------------------------------------------
   ; RUN "makePlates" IN THE SDSS "PLATE" PRODUCT.
   ; The required inputs are the plPlugMapT-$TILE.par files,
   ; plus plPlan.par, plObs.par, plParam.par.
   ; The fiberPlates code selects the guide stars and sky fibers from
   ; those available, and renames COHERENT_SKY/NA objects to OBJECT/SKY.
   ; It also generates the ALIGNMENT holes for each GUIDE fiber.
   ;---------------------------------------------------------------------------

   ;----------
   ; Create the file "plPlan.par" in the current directory.

   cd, current=thisdir
   cd, thisdir
   plhdr = '# Created on ' + systime()
   plhdr = [plhdr, "parametersDir " + paramdir]
   plhdr = [plhdr, "parameters    " + "plParam.par"]
   plhdr = [plhdr, "plObsFile     " + "plObs.par"]
   plhdr = [plhdr, "outFileDir    " + thisdir]
   plhdr = [plhdr, "tileDir       " + thisdir]
   yanny_write, 'plPlan.par', hdr=plhdr

   ;----------
   ; Create the file "plObs.par" in the current directory.

   plhdr = '# Created on ' + systime()
   plhdr = [plhdr, "plateRun special"]
   plstructs = ["typedef struct {", $
                "   int plateId;", $
                "   int tileId;", $
                "   float temp;", $
                "   float haMin;", $
                "   float haMax;", $
                "   int mjdDesign", $
                "} PLOBS;"]
   plobs = create_struct(name='PLOBS', $
    'PLATEID'  ,  0L, $
    'TILEID'   ,  0L, $
    'TEMP'     , 0.0, $
    'HAMIN'    , 0.0, $
    'HAMAX'    , 0.0, $
    'MJDDESIGN',  0L)
   plobs = replicate(plobs, ntile)
   plobs.plateid = platenums
   plobs.tileid = tilenums
   plobs.temp = apotemperature
   plobs.mjddesign = thismjd
   yanny_write, 'plObs.par', ptr_new(plobs), hdr=plhdr, structs=plstructs

   print
   print, 'In the "plate" product run the following commands:"
   print, '   makePlates'
   print, '   fiberPlates -skipBrightCheck'
   print, 'Then type ".cont" in IDL to continue.'
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
;      (*pp[0]).fiberid = thisplate ; Store the plate number in FIBERID
      if (itile EQ 0) then allplug = *pp[0] $
       else allplug = [allplug, *pp[0]]
      if (itile EQ 0) then platenums = thisplate $
       else platenums = [platenums, thisplate]
      if (itile EQ 0) then platearr = replicate(thisplate, n_elements(*pp[0])) $
       else platearr = [platearr, replicate(thisplate, n_elements(*pp[0]))]
      yanny_free, pp
   endfor

   ;----------
   ; Select the 11 non-fake guide stars

   iguide = where(allplug.holetype EQ 'GUIDE' AND allplug.throughput GT 0)
   if (n_elements(iguide) NE 11) then $
    message, 'The number of guide fibers is wrong.'
   newplug = allplug[iguide]
   newplatearr = platearr[iguide]

   ; I've stored the correct FIBERID's in the PRIMTARGET field...
   correctid = newplug.primtarget
   newplug.fiberid = correctid
   newplug.primtarget = 0

   ;----------
   ; Find the alignment hole corresponding to each of these guide stars

   for ii=0, 10 do begin
      ; The line below matches based upon PLATE and FIBERID, though that
      ; FIBERID isn't neccessarily the correct one; we force it to be correct.
      jj = where(allplug.holetype EQ 'ALIGNMENT' $
       AND platearr EQ platearr[iguide[ii]] $
       AND allplug.fiberid EQ allplug[iguide[ii]].fiberid, nj)
      if (nj NE 1) then $
       message, 'Wrong number of alignment holes'
      addplug = allplug[jj]
      addplug.fiberid = correctid[ii] ; Force the correct FIBERID
      newplug = [newplug, addplug]
      newplatearr = [newplatearr, platearr[jj]]
   endfor

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
         newplug = design_append(newplug, newplatearr, $
          allplug[iobj[ii]], platearr[iobj[ii]])
         if (n_elements(newplug) GT nbefore) then nadded = nadded + 1
       endif
   endfor
   if (nadded LT 640) then $
    message, 'Fewer than 640 objects'

   ;----------
   ; If the /ADDFUND flag is set, then recast those fundamental standard holes
   ; as such.  Call them SPECTROPHOTO_STD, and give them all the ra,dec
   ; of the star.  This last thing may be a bit confusing, because then
   ; the xfocal,yfocal positions are inconsistent with ra,dec.  But whatever
   ; data we get through these holes (using a smear) will be of that star.

;   indx = where(newplug.throughput EQ fundpriority, nindx)
   indx = where(newplatearr EQ fundplatenum, nindx) ; Should be same as above
   if (nindx GT 0) then begin
      print, 'Renaming ', nindx, ' fundamental standard holes'
      newplug[indx].objtype = 'SPECTROPHOTO_STD'
;      newplug[indx].ra = funddat[ifund].ra
;      newplug[indx].dec = funddat[ifund].dec
   endif

   ;---------------------------------------------------------------------------
   ; DELETE INTERMEDIATE FILES
   ;---------------------------------------------------------------------------

   ;----------
   ; Delete the existing plPlugMapP,plOverlay,plFanuc,plMeas,plDrillPos files

   print, 'Removing old plPlugMapP,plOverlay,plFanuc,plMeas,plDrillPos files'
   for itile=0, n_elements(platenums)-1 do begin
      spawn, '\rm -f ' + plugmappfile[itile]
      platestr = string(platenums[itile],format='(i4.4)')
      spawn, '\rm -f ' + 'plOverlay-'+platestr+'.ps'
      spawn, '\rm -f ' + 'plFanuc-'+platestr+'.ps'
      spawn, '\rm -f ' + 'plMeas-'+platestr+'.ps'
      spawn, '\rm -f ' + 'plDrillPos-'+platestr+'.ps'
   endfor

   ;---------------------------------------------------------------------------
   ; CHANGE THE PLATE NAMES (unless /NORENAME is set)
   ;---------------------------------------------------------------------------

   if (NOT keyword_set(norename)) then begin
      for itile=0, n_elements(platenums)-1 do begin
         pointingname = string(byte(65+itile))
         thisplatename = string(platenums[0],format='(i4.4)') + pointingname

         ; The first tile retains its name, the others get letters appended.
         if (itile GT 0) then $
          plugmappfile[itile] = 'plPlugMapP-' + thisplatename + '.par'

         ; Change the "plateId" in the Yanny header
; Actually don't change this, since it might break SOP!!
;         junk = yanny_par(*(hdrarr[itile]), 'plateId', indx=indx)
;         (*(hdrarr[itile]))[indx] = 'plateId ' + thisplatename

         ; Add a Yanny keyword "pointing" that is A,B,C,...
         *(hdrarr[itile]) = [*(hdrarr[itile]), 'pointing ' + pointingname]
      endfor
   endif

   ;---------------------------------------------------------------------------
   ; SET THE FIBERID'S FOR PLOTTING OVERLAY FILES
   ; These are done with 7 bundles across the lowest part, then 9, 9, 7.
   ;---------------------------------------------------------------------------

   print, 'Setting negative FIBERIDs for "makePlots"'
   iobj = where(newplug.holetype EQ 'OBJECT')
   if (n_elements(iobj) NE 640) then $
    message, 'Fewer than 640 OBJECT positions on this plate!'

;      isort = sort(newplug[iobj].yfocal)
;      newplug[iobj[isort]].fiberid = -lindgen(640) - 1

   newid = lindgen(20) + 1
   isort = sort(newplug[iobj].yfocal)

   ; Loop over 4 swaths in declination
   ysort = sort(newplug[iobj].yfocal)
   for iyswath=0, 3 do begin
      ; Select which objects are in this horizontal band
      if (iyswath EQ 0) then begin
         indx = iobj[ysort[0:139]]
         nxswath = 7
      endif else if (iyswath EQ 1) then begin
         indx = iobj[ysort[140:319]]
         nxswath = 9
      endif else if (iyswath EQ 2) then begin
         indx = iobj[ysort[320:499]]
         nxswath = 9
      endif else if (iyswath EQ 3) then begin
         indx = iobj[ysort[500:639]]
         nxswath = 7
      endif

      ; For those objects in this horizontal band, sort by X
      indx = indx[ sort(newplug[indx].xfocal) ]

      ; Loop through short stretch of X in this horizontal band,
      ; then sort in those stretches by Y
      for ixswath=0, nxswath-1 do begin
         ii = indx[ixswath*20:ixswath*20+19]
         ii = ii[ sort(newplug[ii].yfocal) ]
         newplug[ii].fiberid = -newid ; (set these ID's as negative numbers)
         newid = newid + 20
      endfor
   endfor

   ;---------------------------------------------------------------------------
   ; WRITE THE MODIFIED PLUGMAPP FILES.
   ; This combines objects from the different tiles.
   ; This overwrites files that already exist.
   ;---------------------------------------------------------------------------

   for itile=0, n_elements(platenums)-1 do begin

      thisplate = platenums[itile]
      modplug = newplug
      junk = where(modplug.holetype EQ 'OBJECT' $
       AND newplatearr EQ thisplate, ct)
      print, 'Number of objects on plate ', thisplate, ' = ', ct

      ;----------
      ; Keep only the guide fibers on this tile. ???

      ;----------
      ; Objects that are actually on other tiles are renamed to sky fibers
      ; on this plate.

      iobj = where(modplug.holetype EQ 'OBJECT' $
       AND newplatearr NE thisplate)
      if (indx[0] NE -1) then begin
         modplug[iobj].holetype = 'OBJECT'
         modplug[iobj].objtype = 'SKY'
         modplug[iobj].mag = fakemag
         modplug[iobj].primtarget = 0L
         modplug[iobj].sectarget = 16L
      endif else begin
         message, 'No objects to use as guide fibers for plate '+string(thisplate)
      endelse

      ;----------
      ; Rotate the positions of fibers from other tiles to be where they
      ; will actually be pointing on this tile.

      for itile2=0, n_elements(platenums)-1 do begin
         if (itile2 NE itile) then begin
            indx = where(newplatearr EQ platenums[itile2])
            if (indx[0] NE -1) then begin
               print, 'Rotating ', n_elements(indx), ' objects from plate ', $
                platenums[itile2], ' to plate ', thisplate
               plate_rotate, racen[itile2], deccen[itile2], $
                racen[itile], deccen[itile], $
                modplug[indx].ra, modplug[indx].dec, tmpra, tmpdec
               modplug[indx].ra = tmpra
               modplug[indx].dec = tmpdec
            endif
         endif
      endfor

      ;----------
      ; Alignment holes seem to always be given RA=0, DEC=0

      indx = where(strtrim(modplug.holetype) EQ 'ALIGNMENT')
      if (indx[0] EQ -1) then $
       message, 'Alignment holes not found'
      modplug[indx].ra = 0
      modplug[indx].dec = 0

      ;----------
      ; Add the one quality hole which we already have, which is the
      ; center hole

      iqual = where(allplug.holetype EQ 'QUALITY' $
       AND platearr EQ thisplate, nqual)
      if (nqual NE 1) then $
      message, 'We expect 1 quality hole already (the center hole)'
      modplug = [modplug, allplug[iqual]]

      ; Must set the following to avoid "fiberPlates" from crashing!!
      modplug.spectrographid = -9999
      modplug.throughput = -9999

      yanny_write, plugmappfile[itile], ptr_new(modplug), $
       hdr=*(hdrarr[itile]), enums=plugenum, structs=plugstruct
   endfor

   ;---------------------------------------------------------------------------
   ; CREATE THE DRILL FILES FROM THE 1ST PLATE.
   ; (Drill files from any of the plates would be identical.)
   ; Run the code makeFanuc, makeDrillPos, use_cs3.
   ;---------------------------------------------------------------------------

   ;----------
   ; Modify the "plObs.par" file, so that we generate drill files only
   ; from the first tile/plate.

   yanny_read, 'plObs.par', plobs, hdr=plhdr, structs=plstructs
   plobs = (*plobs)[0]
   yanny_write, 'plObs.par', ptr_new(plobs), hdr=plhdr, structs=plstructs

   print
   print, 'In the "plate" product run the following commands:"
   print, '   makeFanuc'
   print, '   makeDrillPos'
   print, '   use_cs3'
   print, '   makePlots -skipBrightCheck'
   print, 'Then you are done!'

   return
end
;------------------------------------------------------------------------------
