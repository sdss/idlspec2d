;+
; NAME:
;   design_multiplate
;
; PURPOSE:
;   Routine to design a plate with several pointings.
;
; CALLING SEQUENCE:
;   design_multiplate, stardata, [ tilenums=, racen=, deccen=, $
;    guidetiles= ]
;
; INPUTS:
;   stardata   - Structure with data for each star; must contain the
;                fields RA, DEC, MAG[5], PRIORITY, HOLETYPE, TILENUM.
;                HOLETYPE can be either 'OBJECT' or 'GUIDE'; objects with
;                other values are ignored.
;                Stars with TILENUM=0 can only be used as guide stars.
;   racen      - RA center for each tile
;   deccen     - DEC center for each tile
;
; OPTIONAL INPUTS:
;   tilenums   - Array of tile numbers; default to the unique list of
;                tile numbers in STARDATA.TILENUM, but exclude the number zero.
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
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_avsigclip()
;
; INTERNAL SUPPORT ROUTINES:
;   spbiasgen1
;
; REVISION HISTORY:
;   21-Nov-2001  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro design_multiplate, stardata, tilenums=tilenums, $
 racen=racen, deccen=deccen, guidetiles=guidetiles

   if (NOT keyword_set(tilenums)) then begin
      tilenums = stardata.tilenum
      tilenums = tilenums[ uniq(tilenums, sort(tilenums)) ]
      ; Exclude any tile numbers of zero
      tilenums = tilenums[ where(tilenums NE 0) ]
   endif

   ntile = n_elements(tilenums)
   if (n_elements(racen) NE ntile OR n_elements(deccen) NE ntile) then $
    message, 'Number of elements in RACEN,DECCEN must agree w/ number of tiles'

   if (NOT keyword_set(guidetiles)) then begin
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
   endif
   if (N_elements(guidetiles) NE 11) then $
    message, 'GUIDETILES must be an 11-element array'

   if (NOT keyword_set(racen) OR NOT keyword_set(deccen)) then $
    message, 'RACEN,DECCEN must be specified'

   fakemag = 25.0 ; Magnitudes for all fake objects

   ;----------
   ; Read a template plugmap structure

   yanny_read, 'plPlugMapT-XXXX.par', pp, $
    hdr=plughdr, enums=plugenum, structs=plugstruct
   blankplug = *pp[0]
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

splog, 'Assigning real guide fiber number ', iguide+1
            indx = where(strtrim(stardata.holetype,2) EQ 'GUIDE')
            adiff = djs_diff_angle(guidera[iguide], guidedec[iguide], $
             stardata[indx].ra, stardata[indx].dec)
            junk = min(adiff, imin)
            addplug = blankplug
            struct_assign, stardata[indx[imin]], addplug
            addplug.throughput = 2L

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
      addplug.throughput = 1L
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
      addplug.objtype = 'SKY'
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
      ; Remove objects more than 1.49 deg from the center

      adiff = djs_diff_angle(thisracen, thisdeccen, allplug.ra, allplug.dec)
      allplug = allplug[ where(adiff LT 1.49) ]

      ;----------
      ; Resolve conflicts

      ; First re-sort everything from highest priority to lowest, so
      ; that we can always keep the 1st object in the list in the event
      ; of a conflict.

      isort = reverse(sort(allplug.throughput))
      allplug = allplug[isort]
      allplug.throughput = 0 ; Zero this out now that we've used it.

      ; Resolve conflicts with center hole ??? Separate by 68 arc sec ???

      ; Loop through all objects, discarding close pairs (within 55 arcsec)
      ; Always keep the 1st of any close group of objects.
      ; This is a stupidly slow N^2 implementation...

      iobj = 0L
      while (iobj LT n_elements(allplug)-1) do begin
         adiff = djs_diff_angle(allplug[iobj].ra, allplug[iobj].dec, $
          allplug.ra, allplug.dec)
         igood = where(adiff GT 55.0/3600.0 $
          OR lindgen(n_elements(allplug)) EQ iobj, ngood)
;         if (ngood LT n_elements(allplug)) then $
;          print, 'Discarding ', n_elements(allplug)-ngood, ' objects'
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

      outfile = 'plPlugMapT-' + string(thistilenum,format='(i4.4)') + '.par'
      yanny_write, outfile, ptr_new(allplug), hdr=outhdr, $
       enums=plugenum, structs=plugstruct
stop

   endfor ; End loop over pointing number

stop
   return
end

