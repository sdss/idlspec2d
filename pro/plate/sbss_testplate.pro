; Generate a test plate for Bright Stars as a plPlugMapT files
; that can go directly through the "plate" product code.
; Call it tile #9001, plate #9001.

; Exclusion zone for guide fibers = 6.995 mm = 115.654 arcsec
; Exclusion zone for alignment holes = 3.327 mm = 55.0079 arcsec

pro sbss_testplate

   plateid = 9001L
   tileid = 9001L
   setenv, 'HIPPARCOS_DIR=/u/schlegel/hipparcos'
;   setenv, 'BOSS_SPECTRO_REDUX=/u/dss/spectro'

   ;----------
   ; Read plug-map info for plate 360, which we will use for selecting
   ; spectro-photo stars and for a template of the plPlugMap file format.

;   readspec, 360, plug=plug360
   yanny_read, '/u/dss/astrolog/51780/plPlugMapM-0360-51774-01.par', $
    pp, hdr=plughdr, enums=plugenum, structs=plugstruct
   plug360 = *pp[0]

   ;----------
   ; Create a blank plug-map

   blankplug = plug360[0]
   struct_assign, {junk:0}, blankplug

   ;---------------------------------------------------------------------------
   ; GUIDE FIBERS
   ;---------------------------------------------------------------------------

   ;----------
   ; Select initial guesses for pointing centers:
   ;   Four at DEC=+55
   ;   Four at DEC=+53
   ;   One at the position of plates 360 and 362

   racen =  [265.0, 267.0, 269.0, 271.0, 265.0, 267.0, 269.0, 271.0, 264.3484]
   deccen = [ 55.0,  55.0,  55.0,  55.0,  53.0,  53.0,  53.0,  53.0,  54.53896]
   npoint = n_elements(racen)

   ;----------
   ; Make list of which pointing each guide fiber is assigned to

   gfiber = create_struct( $
    'pointing' , 0L, $
    'fibernum' , 0L, $
    'ra'       , 0.d, $
    'dec'      , 0.d, $
    'hip'      , 0L , $
    'tyc'      , 0L , $
    'vmag'     , 0.0, $
    'xreach'   , 0.0, $
    'yreach'   , 0.0, $
    'rreach'   , 0.0, $
    'xprefer'  , 0.d, $
    'yprefer'  , 0.d )
   gfiber = replicate(gfiber, 11)
;   gfiber.pointing = [0,0,0,0,1,1,2,4,4,5,8]
   ; Put the last two guide fibers (which are bigger) in the 0th pointing
   gfiber.pointing = [5,8,0,0,1,1,2,4,4,0,0]
   igolden = 9 ; This says that gfiber[igolden] is the golden guide fiber
   gfiber.fibernum = lindgen(11) + 1

   ;----------
   ; Read the list of possible guide stars from the Tycho catalog.
   ; Trim to some subset now for simplicity in a boneheaded manner.

   tycdat = tyc_read(/small)
   indx = where(tycdat.radeg GT min(racen)-2.0 $
            AND tycdat.radeg LT max(racen)+2.0 $
            AND tycdat.dedeg GT min(deccen)-2.0 $
            AND tycdat.dedeg LT max(deccen)+2.0 $
            AND tycdat.vmag GE 10.0 $
            AND tycdat.vmag LE 12.5 )
   tycdat = tycdat[indx]

   ;----------
   ; Assign a golden fiber to each pointing, and
   ; tweak the pointing centers.

   for ipoint=0, npoint-1 do begin
      adiff = djs_diff_angle(racen[ipoint], deccen[ipoint], $
       tycdat.radeg, tycdat.dedeg)
      mindiff = min(adiff, ibest)
print, 'Shifting pointing number ', ipoint, ' by ', mindiff, ' deg'

      ; The plate center will be 10 arc min **above** this star
      ; Check this!!!
      racen[ipoint] = tycdat[ibest].radeg
      deccen[ipoint] = tycdat[ibest].dedeg + 10.0/60

      ; The golden fiber is assigned to pointing #0, but will
      ; also be used for all other pointings.
      if (ipoint EQ 0) then begin
         gfiber[igolden].pointing = 0
         gfiber[igolden].ra = tycdat[ibest].radeg
         gfiber[igolden].dec = tycdat[ibest].dedeg
         gfiber[igolden].vmag = tycdat[ibest].vmag
         gfiber[igolden].hip = tycdat[ibest].hip
         gfiber[igolden].tyc = tycdat[ibest].tyc
      endif
   endfor

   ;----------
   ; The following info is from the "plate" product in the
   ; file "$PLATE_DIR/test/plParam.par".
   ;   XREACH,YREACH = Center of the fiber reach [mm]
   ;   RREACH = Radius of the fiber reach [mm]
   ;   XPREFER,YREACH = Preferred position for the fiber [mm]
   ; Note that the plate scale is approx 217.7358 mm/degree.
   ; Moving +RA is +XFOCAL, +DEC is +YFOCAL.

   platescale = 217.7358
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

   ;----------
   ; Position each guide fiber (except we've already done the golden one)

   for iguide=0, 10 do begin
      if (iguide NE igolden) then begin

         ; This pointing is centered here...
         thisracen = racen[gfiber[iguide].pointing]
         thisdeccen = deccen[gfiber[iguide].pointing]

         ; This guide fiber will now want to be at
         ; the plate position (gfiber[iguide].xprefer, gfiber[iguide].yprefer)
         ; which we convert to coordinates (thisra1,thisdec1)
         dphi = gfiber[iguide].xprefer / platescale
         dtheta = gfiber[iguide].yprefer / platescale
         thisra1 = thisracen + dphi / cos(thisdeccen/!radeg) ; approximate
         thisdec1 = thisdeccen + dtheta

         ; Now find the Tycho star closest to (thisra,thisdec)
         ; Give a huge penalty if more than 1.49 deg from the plate center
         adiff = djs_diff_angle(thisra1, thisdec1, $
                                tycdat.radeg, tycdat.dedeg) $
          + 999 * (djs_diff_angle(thisracen, thisdeccen, $
                   tycdat.radeg, tycdat.dedeg) GT 1.49)

         mindiff = min(adiff, ibest)
print,'Guide star #', iguide+1, ' offset from optimal = ', mindiff, ' deg'
; Need to check that this position is within the paramters!!!???

         ; Fill in the information for this guide fiber
         gfiber[iguide].ra = tycdat[ibest].radeg
         gfiber[iguide].dec = tycdat[ibest].dedeg
         gfiber[iguide].vmag = tycdat[ibest].vmag
         gfiber[iguide].hip = tycdat[ibest].hip
         gfiber[iguide].tyc = tycdat[ibest].tyc
      endif
   endfor

   ;----------
   ; Create the output structure and put the guide fibers into it

   addplug = replicate(blankplug, 11)
   addplug.objid = transpose([ $
    [replicate(plateid,11)], $ ; Plate number
    [gfiber.pointing]      , $ ; Pointing number
    [gfiber.hip]           , $ ; Hipparcos ID
    [gfiber.tyc]           , $ ; Tycho ID
    [replicate(0L,11)]])       ; (Unused)
   addplug.holetype = 'GUIDE'
   addplug.ra = gfiber.ra
   addplug.dec = gfiber.dec
   ; Use the V magnitude for all 5 magnitudes in this file
   addplug.mag = gfiber.vmag ## replicate(1L, 5)
   addplug.objtype = 'NA'
   addplug.fiberid = gfiber.fibernum ; These are numbered 1 through 11
   addplug.primtarget = gfiber.pointing ; Put pointing num here too
   addplug.sectarget = 64L
   allplug = addplug

   ;---------------------------------------------------------------------------
   ; SCIENCE FIBERS
   ;---------------------------------------------------------------------------

   ;----------
   ; Select the Hipparcos (science) targets
   ; Trim to some subset now for simplicity in a boneheaded manner.
   ; For now, trim to V <= 11.0.

   hipdat = hip_read(/small)
   indx = where(hipdat.radeg GT min(racen)-2.0 $
            AND hipdat.radeg LT max(racen)+2.0 $
            AND hipdat.dedeg GT min(deccen)-2.0 $
            AND hipdat.dedeg LT max(deccen)+2.0 $
            AND hipdat.vmag NE  0.0 $
            AND hipdat.vmag LE 11.0 )
   hipdat = hipdat[indx]

   ;----------
   ; Loop through each pointing and assign science targets
   ; Don't do the last pointing, which is reserved for plate 360 objects

   for ipoint=0, npoint-2 do begin
      adiff = djs_diff_angle(racen[ipoint], deccen[ipoint], $
       hipdat.radeg, hipdat.dedeg)
      iadd = where(adiff LT 1.49, nadd)
      if (nadd EQ 0) then $
       message, 'Crap!  No stars in this pointing!'

      addplug = replicate(blankplug, nadd)
      addplug.objid = transpose([ $
       [replicate(plateid,nadd)], $  ; Plate number
       [replicate(ipoint,nadd)], $     ; Pointing number
       [hipdat[iadd].hip], $           ; Hipparcos ID
       [hipdat[iadd].tyc], $           ; Tycho ID
       [replicate(0L,nadd)]])          ; (Unused)
      addplug.holetype = 'OBJECT'
      addplug.ra = hipdat[iadd].radeg
      addplug.dec = hipdat[iadd].dedeg
      ; Use the V magnitude for all 5 magnitudes in this file
      addplug.mag = hipdat[iadd].vmag ## replicate(1L, 5)
      addplug.objtype = 'STAR_BHB'
      addplug.primtarget = replicate(ipoint,nadd) ; Put pointing num here too

      allplug = [allplug, addplug]
   endfor

   ;----------
   ; Assign science targets to the last pointing from plate 360
   ; This will include the 8 SPECTROPHOTO_STD and 8 REDDEN_STD

   iadd = where( strtrim(plug360.objtype,2) EQ 'SPECTROPHOTO_STD' $
    OR strtrim(plug360.objtype,2) EQ 'REDDEN_STD', nadd)
;    OR (strtrim(plug360.holetype,2) EQ 'OBJECT' AND plug360.mag[2] LT 16.5), nadd)
   addplug = plug360[iadd]
   addplug.primtarget = replicate(npoint-1,nadd) ; Put pointing num here

   ; Check the positions and modify if necessary -- what a kludge!!!???
   ; This happens because the plate center isn't exactly where we want it.
   adiff = djs_diff_angle(racen[npoint-1], deccen[npoint-1], $
    addplug.ra, addplug.dec)
   for i=0, n_elements(addplug)-1 do begin
      if (adiff[i] GT 1.49) then begin
print, 'Kludging an object position!!'
         if (addplug[i].ra GT racen[npoint-1]) then $
          addplug[i].ra = addplug[i].ra - 0.25 $
         else $
          addplug[i].ra = addplug[i].ra + 0.25
         if (addplug[i].dec GT deccen[npoint-1]) then $
          addplug[i].dec = addplug[i].dec - 0.25 $
         else $
          addplug[i].dec = addplug[i].dec + 0.25
      endif
   endfor

   allplug = [allplug, addplug]

   ;---------------------------------------------------------------------------
   ; ADD A LINE OF HOLES IN DECLINATION NEXT TO THE GOLDEN GUIDE FIBER
   ; This will be used for something like "bd17smear".
   ;---------------------------------------------------------------------------

   ;----------
   ; Add exactly 10 such holes.
   ; Say this is also pointing #0.

   nadd = 10
   ipoint = 0
   thisra = replicate(gfiber[igolden].ra, nadd)
   thisdec = gfiber[igolden].dec + [-5,-4,-3,-2,-1,1,2,3,4,5]/60.0
   addplug = replicate(blankplug, nadd)

   addplug.objid = transpose([ $
    [replicate(plateid,nadd)], $ ; Plate number
    [replicate(ipoint,nadd)] , $ ; Pointing number
    [replicate(0L,nadd)]     , $ ; (Unused)
    [replicate(0L,nadd)]     , $ ; (Unused)
    [replicate(0L,nadd)]])       ; (Unused)
   addplug.holetype = 'OBJECT'
   addplug.ra = thisra
   addplug.dec = thisdec
   ; Use the V magnitude for this guide fiber
   addplug.mag = gfiber[igolden].vmag ## replicate(1L, 5)
   addplug.objtype = 'STAR_BHB'
   addplug.primtarget = replicate(ipoint,nadd) ; Put pointing num here too

   allplug = [allplug, addplug]

   ;---------------------------------------------------------------------------
   ; ADD A BUNCH OF SKY FIBERS (in random positions)
   ;---------------------------------------------------------------------------

   nadd = 1000
   ipoint = 0
   addplug = replicate(blankplug, nadd)
   addplug.ra = 3.0 * randomu(2345, nadd) - 1.5 + racen[0]
   addplug.dec = 3.0 * randomu(3456, nadd) - 1.5 + deccen[0]
   addplug.mag = 25.0 ## replicate(1L, 5) ; Call these 25th-mag objects
   ; Set to COHERENT_SKY,NA for a plPlugMapT file
   addplug.holetype = 'COHERENT_SKY'
   addplug.objtype = 'NA'
   ; Set to OBJECT,SKY for a plPlugMapP file
;   addplug.holetype = 'OBJECT'
;   addplug.objtype = 'SKY'
   addplug.primtarget = replicate(ipoint,nadd) ; Put pointing num here too
   addplug.sectarget = 16L

   ; Trim away all skies that are not within 1.49 deg of the plate center
   adiff = djs_diff_angle(racen[ipoint], deccen[ipoint], $
    addplug.ra, addplug.dec)
   ikeep = where(adiff LT 1.49)
   addplug = addplug[ikeep]

   allplug = [allplug, addplug]

   ;---------------------------------------------------------------------------
   ; ROTATE COORDINATES -> All to the 0th pointing
   ;---------------------------------------------------------------------------

   ;----------
   ; First, save the real coordinates in the fields STARL,EXPL !!
   ; This is just a bizarre solution to avoid losing this information.

   allplug.starl = allplug.ra
   allplug.expl = allplug.dec

   ;----------
   ; Loop through each pointing number (saved in the PRIMTARGET field!)

   for ipoint=0, max(allplug.primtarget) do begin
      ; Select objects in this pointing
      indx = where(allplug.primtarget EQ ipoint)

      ; Now these objects must be rotated from their actual plate
      ; center at RACEN[IPOINT],DECCEN[IPOINT] to RACEN[0],DECCEN[0]
      plate_rotate, racen[ipoint], deccen[ipoint], racen[0], deccen[0], $
       allplug[indx].ra, allplug[indx].dec, ra_new, dec_new
      allplug[indx].ra = ra_new
      allplug[indx].dec = dec_new
   endfor

   ;---------------------------------------------------------------------------
   ; RESOLVE CONFLICTS
   ;---------------------------------------------------------------------------

; Resolve conflicts with center hole ??? Separate by 68 arc sec

   ;----------
   ; Loop through all objects, discarding close pairs (within 55 arcsec)
   ; Always keep the 1st of any close group of objects.
   ; This is a stupidly slow N^2 implementation...

   iobj = 0L
   while (iobj LT n_elements(allplug)-1) do begin
      adiff = djs_diff_angle(allplug[iobj].ra, allplug[iobj].dec, $
       allplug.ra, allplug.dec)
      igood = where(adiff GT 55.0/3600.0 $
       OR lindgen(n_elements(allplug)) EQ iobj, ngood)
;      if (ngood LT n_elements(allplug)) then $
;       print, 'Discarding ', n_elements(allplug)-ngood, ' objects'
      allplug = allplug[igood]
      iobj = iobj + 1L
   endwhile

   ;----------
   ; Here's the dumbest thing to do -- pick the first 651 objects,
   ; which will be 11 guide + 640 objects.

;   allplug = allplug[0:650]
   allplug = allplug[0:699] ; Add extra skies in case of collisions !!!???

   ;---------------------------------------------------------------------------
   ; WRITE THE FILE
   ;---------------------------------------------------------------------------

   ; Put a running number in the last OBJID field -- for Yanny's benefit
   allplug.objid[4] = lindgen(n_elements(allplug)) + 1

   outhdr = ['completeTileVersion   v1_0', $
             'tileId ' + string(tileid), $
             'raCen ' + string(racen[0]), $
             'decCen ' + string(deccen[0]) ]
   for ipoint=0, npoint-1 do begin
      outhdr = [outhdr, $
       'raCen'+strtrim(string(ipoint),2)+ '  '+string(racen[ipoint]) ]
      outhdr = [outhdr, $
       'decCen'+strtrim(string(ipoint),2)+ '  '+string(deccen[ipoint]) ]
   endfor
   outfile = 'plPlugMapT-' + string(plateid,format='(i4.4)') + '.par'
   yanny_write, outfile, ptr_new(allplug), hdr=outhdr, $
    enums=plugenum, structs=plugstruct

   adiff = djs_diff_angle(racen[0], deccen[0], allplug.ra, allplug.dec)
print, 'Min distance from plate center = ', min(adiff)*3600.0, ' arcsec'
print, 'Max distance from plate center = ', max(adiff), ' deg'

stop

splot,allplug.ra,allplug.dec,ps=4
i=where(strtrim(allplug.objtype,2) eq 'REDDEN_STD')
soplot,allplug[i].ra,allplug[i].dec,ps=4,color='red'
j=where(strtrim(allplug.objtype,2) eq 'SPECTROPHOTO_STD')
soplot,allplug[j].ra,allplug[j].dec,ps=4,color='green'
k=where(strtrim(allplug.holetype,2) eq 'GUIDE')
for kk=0,10 do sxyouts, allplug[k[kk]].ra, allplug[k[kk]].dec, $
 strtrim(string(allplug[k[kk]].fiberid),2),color='blue', charsize=2

end
