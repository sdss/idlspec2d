;+
; NAME:
;   design_plate
;
; PURPOSE:
;   Routine to design a single plate.
;
; CALLING SEQUENCE:
;   design_plate, stardata, [ racen=, deccen=, tilenum=, platenum=, $
;    airtemp=, nstd=, nminsky=, nextra= ]
;
; INPUTS:
;   stardata   - Structure with data for each star; must contain the
;                fields RA, DEC, MAG[5], HOLETYPE, OBJTYPE.
;                HOLETYPE should be either 'OBJECT' or 'GUIDE'.
;                Special values of OBJTYPE are 'SKY' and 'SPECTROPHOTO_STD'.
;                Other elements in the structure (such as PRIMTARGET,...)
;                will be copied into the output plug-map files.
;   racen      - RA center for tile
;   deccen     - DEC center for tile
;
; OPTIONAL INPUTS:
;   tilenum    - Tile number; default to 1.
;   platenum   - Plate number; default to the same as TILENUMS.
;   airtemp    - Design temperature for APO; default to 5 deg C.
;   nstd       - Number of spectro-photo standards; default to 8.
;                This will be split with NSTD/2 standards on the North
;                half of the plate, and the same number on the South.
;   nminsky    - Minimum number of sky fibers; default to 32.
;   nextra     - Number of extra fibers to assign, beyond the 640 science
;                fibers.  In principle, this should be 0.  In practice,
;                the "fiberPlates" code in the PLATE product will either
;                crash or multiply-assign fibers to the same object unless
;                there are extra objects.  Default to 10.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Priorities can be specified with PRIORITY in the input structure;
;   otherwise, they are assigned randomly.
;
;   We always assign guide stars first, then spectro-photo standards,
;   then objects, then skies.  So the priorities are relevant only
;   within each of those categories.
;
;   This script generates the following files:
;     plObs.par
;     plPlan.par
;     plPlugMapT-$TILE.par
;   The commands from the SDSS "plate" product generate the following files:
;     makePlates - Generate plPlugMapP-$PLATE.par
;     fiberPlates - Generate plOverlay-$PLATE.par, and **overwrite**
;                   the file plPlugMapP-$PLATE.par
;     makeFanuc - Generate plFanuc-$PLATE.par
;     makeDrillPos - Generate plMeas-$PLATE.par, plDrillPos-$PLATE.par
;     use_cs3 - Generates no files.  Simply reports fiber collisions.
;     makePlots - Generate plOverlay-$PLATE.ps
;
; EXAMPLES:
;
; BUGS:
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
;   radec_to_xyfocal
;   design_append()
;
; REVISION HISTORY:
;   14-Jan-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
; Approximately transform RA,DEC -> XFOCAL,YFOCAL.

pro radec_to_xyfocal, ra, dec, xfocal, yfocal, $
 racen=racen, deccen=deccen, airtemp=airtemp

   platescale = 217.7358 ; mm/degree

   radec_to_munu, ra, dec, mu, nu, node=racen-90, incl=deccen
   xfocal = platescale * (mu - racen)
   yfocal = platescale * nu

   return
end
;------------------------------------------------------------------------------
; Search for conflicts between an existing list of drill holes and one
; more potential object.  Return the existing list with the new object
; appended if there was no conflict.

function design_append, allplug, oneplug, nadd=nadd

   if (n_elements(oneplug) NE 1) then $
    message, 'ONEPLUG must contain only one element'
   if (NOT keyword_set(oneplug)) then oneplug = 0
   nadd = 0L

   platescale = 217.7358 ; mm/degree
   if (strtrim(oneplug.holetype) EQ 'GUIDE') then morespace = 7.5 $ ; in mm
    else morespace = 0.

   ;----------
   ; If this is the 1st object in the list, then we can always keep it

   if (NOT keyword_set(allplug)) then begin
      nadd = 1L
      return, oneplug
   endif

   ;----------
   ; Discard objects within 68 arcsec of the center hole or more
   ; then 1.49 deg from the plate center (pad to 75 arcsec and 1.485 deg).

   thisrad = sqrt(oneplug.xfocal^2 + oneplug.yfocal^2)
   if (thisrad LE platescale*75./3600.+morespace $
    OR thisrad GE platescale*1.485-morespace) then return, allplug

   ;----------
   ; Discard objects within 55 arcsec of existing objects. (Pad to 70 arcsec.)
   ; Do this based upon XFOCAL,YFOCAL positions.
   ; The closest two fibers can be is PLATESCALE * 55/3600(deg) = 3.32652 mm

   if (keyword_set(allplug)) then begin
      r2 = (allplug.xfocal - oneplug.xfocal)^2 $
           + (allplug.yfocal - oneplug.yfocal)^2
      mindist = min(sqrt(r2))
      if (mindist LT platescale*70./3600.+morespace) then return, allplug
   endif

   ;----------
   ; Discard objects within 7.0 mm of existing guide fibers.
   ; (Pad this to 10 mm just to be extra careful, since the downstream PLATE
   ; code is so stupid.)

   if (keyword_set(allplug)) then begin
      iguide = where(strtrim(allplug.holetype) EQ 'GUIDE', ct)
      if (ct GT 0) then begin
         r2 = (allplug[iguide].xfocal - oneplug.xfocal)^2 $
              + (allplug[iguide].yfocal - oneplug.yfocal)^2
         mindist = min(sqrt(r2))
         if (mindist LT 10.0+morespace) then return, allplug
      endif
   endif

   nadd = 1L
   return, [allplug, oneplug]
end

;------------------------------------------------------------------------------
pro design_plate, stardata1, tilenums=tilenum, platenums=platenum, $
 racen=racen, deccen=deccen, airtemp=airtemp, nstd=nstd, nminsky=nminsky, $
 nextra=nextra

   if (NOT keyword_set(tilenum)) then tilenum = 1
   if (NOT keyword_set(platenum)) then platenum = tilenum
   if (NOT keyword_set(racen) OR NOT keyword_set(deccen)) then $
    message, 'RACEN,DECCEN must be specified'
   if (n_elements(airtemp) EQ 0) then airtemp = 5.0
   if (NOT keyword_set(nstd)) then nstd = 8L
   if (NOT keyword_set(nminsky)) then nminsky = 32L
   if (n_elements(nextra) EQ 0) then nextra = 10L
   ntot = 651L + nextra ; Number of guide + science fibers

   plugmaptfile = 'plPlugMapT-' + string(tilenum,format='(i4.4)') + '.par'
   plugmappfile = 'plPlugMapP-' + string(platenum,format='(i4.4)') + '.par'
   maxpriority = 2L^31 - 1 ; Maximum value; this is the value for GUIDE stars
   paramdir = concat_dir(getenv('IDLSPEC2D_DIR'), 'examples')

   ;----------
   ; If the priorities of targets are not specified in the input structure,
   ; then assign random priorities between 1 and 100.

   if ((where(tag_names(stardata1) EQ 'PRIORITY'))[0] NE -1) then $
    priority = (stardata1.priority > 1L) < (maxpriority-2) $
   else $
    priority = long(randomu(24680, n_elements(stardata1)) * 100) + 1

   ;----------
   ; Add the tags XFOCAL,YFOCAL to the structure of object data.

   radec_to_xyfocal, stardata1.ra, stardata1.dec, xfocal, yfocal, $
    racen=racen, deccen=deccen, airtemp=airtemp
   xydata = replicate(create_struct('XFOCAL', 0L, 'YFOCAL', 0L), $
    n_elements(stardata1))
   xydata.xfocal = xfocal
   xydata.yfocal = yfocal
   stardata = struct_addtags(stardata1, xydata)

   ;----------
   ; Read a template plugmap structure

   blankplug = (yanny_readone( $
    filepath('plPlugMapT-XXXX.par', root_dir=paramdir), pp, $
    hdr=plughdr, enums=plugenum, structs=plugstruct))[0]
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

   ;----------
   ; Add guide fibers

   for iguide=0, nguide-1 do begin
      ;----------
      ; Assign the nearest available guide fiber

      print, 'Assigning guide fiber number ', iguide+1

      nadd1 = 0
      while (nadd1 EQ 0) do begin
         indx = where(strtrim(stardata.holetype,2) EQ 'GUIDE' $
          AND priority GT 0, ct)
         if (ct EQ 0) then $
          message, 'Ran out of guide stars!'
         adiff = djs_diff_angle( $
          gfiber[iguide].xprefer, gfiber[iguide].yprefer, $
          stardata[indx].xfocal, stardata[indx].yfocal)
         junk = min(adiff, ibest)

         addplug = blankplug
         struct_assign, stardata[indx[ibest]], addplug
         addplug.holetype = 'GUIDE'
         addplug.objtype = 'NA'
         addplug.sectarget = 64L

         allplug = design_append(allplug, addplug, nadd=nadd1)

         priority[indx[ibest]] = 0 ; Don't try to target again
      endwhile
   endfor

   ;----------
   ; Add spectro-photo standards

   for nsouth=-1, 1, 2 do begin ; First do South, then North half of plate
      nadd = 0L
      while (nadd LT nstd/2.) do begin
         indx = where(strtrim(stardata.holetype,2) EQ 'OBJECT' $
          AND strtrim(stardata.objtype,2) EQ 'SPECTROPHOTO_STD' $
          AND nsouth*stardata.yfocal GE 0 $
          AND priority GT 0, ct)
         if (ct EQ 0) then $
          message, 'Ran out of spectro-photo stars!'

         junk = max(priority[indx], ibest)
         addplug = blankplug
         struct_assign, stardata[indx[ibest]], addplug
         addplug.holetype = 'OBJECT'
         addplug.objtype = 'SPECTROPHOTO_STD'
         addplug.sectarget = 32L

         allplug = design_append(allplug, addplug, nadd=nadd1)
         nadd = nadd + nadd1

         priority[indx[ibest]] = 0 ; Don't try to target again
      endwhile
   endfor

   ;----------
   ; Add objects that are neither sky nor objects

   while (n_elements(allplug) LT ntot - nminsky) do begin
      indx = where(strtrim(stardata.holetype,2) EQ 'OBJECT' $
       AND strtrim(stardata.objtype,2) NE 'SKY' $
       AND priority GT 0, ct)
      if (ct EQ 0) then $
       message, 'Ran out of non-sky object targets!'

      junk = max(priority[indx], ibest)
      addplug = blankplug
      struct_assign, stardata[indx[ibest]], addplug
      addplug.holetype = 'OBJECT'

      allplug = design_append(allplug, addplug)

      priority[indx[ibest]] = 0 ; Don't try to target again
   endwhile

   ;----------
   ; Add sky fibers

   while (n_elements(allplug) LT ntot) do begin
      indx = where(strtrim(stardata.holetype,2) EQ 'OBJECT' $
       AND strtrim(stardata.objtype,2) EQ 'SKY' $
       AND priority GT 0, ct)
      if (ct EQ 0) then $
       message, 'Ran out of sky targets!'

      junk = max(priority[indx], ibest)
      addplug = blankplug
      struct_assign, stardata[indx[ibest]], addplug
      addplug.holetype = 'COHERENT_SKY'
      addplug.objtype = 'NA'
      addplug.sectarget = 16L

      allplug = design_append(allplug, addplug)

      priority[indx[ibest]] = 0 ; Don't try to target again
   endwhile

   ;----------
   ; Write the plPlugMapT file

   outhdr = ['completeTileVersion   v1_0', $
             'tileId ' + string(tilenum), $
             'raCen ' + string(racen), $
             'decCen ' + string(deccen) ]

   allplug.xfocal = -999 ; These values will be computed by PLATE
   allplug.yfocal = -999 ; These values will be computed by PLATE
   allplug.spectrographid = -999L ; These values will be computed by PLATE
   allplug.fiberid = -999L ; These values will be computed by PLATE
   allplug.throughput = 1L ; These values will be computed by PLATE
   yanny_write, plugmaptfile, ptr_new(allplug), hdr=outhdr, $
    enums=plugenum, structs=plugstruct

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
    'PLATEID'  , platenum, $
    'TILEID'   , tilenum, $
    'TEMP'     , airtemp, $
    'HAMIN'    , 0.0, $
    'HAMAX'    , 0.0, $
    'MJDDESIGN', current_mjd())
   yanny_write, 'plObs.par', ptr_new(plobs), hdr=plhdr, structs=plstructs

   print
   print, 'In the "plate" product run the following commands:"
   print, '   makePlates'
   print, '   fiberPlates -skipBrightCheck'
   print, '   makeFanuc'
   print, '   makeDrillPos'
   print, '   use_cs3'
   print, '   makePlots -skipBrightCheck'
   print, 'Then you are done!'

   return
end
;------------------------------------------------------------------------------
