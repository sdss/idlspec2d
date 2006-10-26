;+
; NAME:
;   design_plate
;
; PURPOSE:
;   Routine to design a single plate.
;
; CALLING SEQUENCE:
;   design_plate, stardata, [ racen=, deccen=, tilenum=, platenum=, $
;    airtemp=, nstd=, nminsky=, ngtarg=, /southern ]
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
;   nstd       - Number of spectro-photo standards; default to 16.
;                This will be split with NSTD/2 standards on the North
;                half of the plate, and the same number on the South.
;   nminsky    - Minimum number of sky fibers; default to 32.
; COMMENTS:
;   All non-SKY and non-SPECTROPHOTO_STD objects will be put on the
;   plate. This routine will choose which of the given SKY and
;   SPECTROPHOTO_STD targets to put on the plate.
;
;   If non-SKY and non-SPECTROPHOTO_STD objects collide with one
;   another at all (are < 55'' apart) this routine will halt and
;   report an error.
;
;   If non-SKY and non-SPECTROPHOTO_STD objects are within 68 arcsec
;   of the center, or are outside 1.49 deg, this routine will halt and
;   report an error.
;
;   We always assign then objects, then guide stars, then
;   spectro-photo standards, then skies.  
;
;   This script generates the following files:
;     plObs.par
;     plPlan.par
;     plPlugMapT-$TILE.par
;     plPlugMapP-$PLATE.par
;   The commands from the SDSS "plate" product generate the following files:
;     makePlates - Generate plPlugMapP-$PLATE.par
;     fiberPlates - Generate plOverlay-$PLATE.par, and **overwrite**
;                   the file plPlugMapP-$PLATE.par
;     makeFanuc - Generate plFanuc-$PLATE.par
;     makeDrillPos - Generate plMeas-$PLATE.par, plDrillPos-$PLATE.par
;     use_cs3 - Generates no files.  Simply reports fiber collisions.
;     makePlots - Generate plOverlay-$PLATE.ps
;
;   When there is a problem, the "fiberPlates" script outputs error
;   messages like:
;     collision at 554 = {307 125 6 160 47} OBJECT with...
;     collision at 565 = {307 125 2 168 4015} COHERENT_SKY with...
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   concat_dir()
;   current_mjd()
;   djs_diff_angle()
;   splog
;   yanny_readone()
;   yanny_write
;
; INTERNAL SUPPORT ROUTINES:
;   approx_radec_to_xyfocal
;   design_append()
;
; REVISION HISTORY:
;   14-Jan-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
; Approximately transform RA,DEC -> XFOCAL,YFOCAL.
pro approx_radec_to_xyfocal, ra, dec, xfocal, yfocal, $
                             racen=racen, deccen=deccen, airtemp=airtemp

platescale = 217.7358           ; mm/degree

radec_to_munu, ra, dec, mu, nu, node=racen-90, incl=deccen
xfocal = platescale * (mu - racen)
yfocal = platescale * nu

return
end
;-------------------------------------------------------
; Check geometrical conditions
function geometry_check, data, racen, deccen, guide=guide

platesize=1.49
platecen=68./3600.

holeok=bytarr(n_elements(data))

; check for things in the plate and not too near the center
spherematch, racen, deccen, data.ra, data.dec, $
  platesize, m1, m2, d12
if(m2[0] eq -1) then return
holeok[m2]=1 
iclose=where(d12 lt platecen, nclose)
if(nclose gt 0) then $
  holeok[m2[iclose]]=0

return, holeok

end
;---------------------------------------------------------
; return whether each new data collides with an old hole
function collision_check, holedata, newdata, guide=guide

holesize=holesize(guide=guide)
guidesize=holesize(/guide)

; find which new objects don't collide with old
holeok=bytarr(n_elements(newdata))+1

; first check old guide stars
iguide=where(strtrim(holedata.holetype,2) eq 'GUIDE', nguide)
if(nguide gt 0) then begin
    spherematch, holedata[iguide].ra, holedata[iguide].dec, $
      newdata.ra, newdata.dec, 0.5*(holesize+guidesize), $
      m1, m2, d12
    if(m2[-1] gt -1) then $
      holeok[m2]=0
endif

; first check new non-guide stars
inotguide=where(strtrim(holedata.holetype,2) ne 'GUIDE', nnotguide)
if(nnotguide gt 0) then begin
    spherematch, holedata.ra, holedata.dec, $
      newdata[inotguide].ra, newdata[inotguide].dec, $ 
      0.5*(holesize+guidesize), m1, m2, d12
    if(m2[-1] gt -1) then $
      holeok[inotguide[m2]]=0
endif

return, holeok

end
;------------------------------------------------------------------------------
; Append to current list. (This used to check for collisions, etc; I
; am keeping it as a separate function in case we need to do something
; special --- MRB)
function design_append, allplug, oneplug, nadd=nadd

if (n_elements(oneplug) NE 1) then $
  message, 'ONEPLUG must contain only one element'
if (NOT keyword_set(oneplug)) then oneplug = 0
nadd = 0L

;----------
; If this is the 1st object in the list, then we can always keep it
if (NOT keyword_set(allplug)) then begin
    nadd = 1L
    return, oneplug
endif

;----------
; Or just append it
nadd = 1L
return, [allplug, oneplug]

end

;------------------------------------------------------------------------------
; Main program
pro design_plate, stardata1, tilenums=tilenum, platenums=platenum, $
                  racen=racen, deccen=deccen, airtemp=airtemp, nstd=nstd, $
                  nminsky=nminsky

if (NOT keyword_set(tilenum)) then tilenum = 1
if (NOT keyword_set(platenum)) then platenum = tilenum
if (NOT keyword_set(racen) OR NOT keyword_set(deccen)) then $
  message, 'RACEN,DECCEN must be specified'
if (n_elements(airtemp) EQ 0) then airtemp = 5.0
if (NOT keyword_set(nstd)) then nstd = 16L
if (NOT keyword_set(nminsky)) then nminsky = 32L
ntot = 640L                     ; Number of fibers
nsci = ntot-nstd-nminsky
used=bytarr(n_elements(stardata))

;----------
; Set up outputs
plugmaptfile = 'plPlugMapT-' + string(tilenum,format='(i4.4)') + '.par'
plugmappfile = 'plPlugMapP-' + string(platenum,format='(i4.4)') + '.par'
paramdir = concat_dir(getenv('IDLSPEC2D_DIR'), 'examples')
blankplug = (yanny_readone(filepath('plPlugMapT-XXXX.par', $
                                    root_dir=paramdir), pp, $
                           hdr=plughdr, enums=plugenum, $
                           structs=plugstruct))[0]
struct_assign, {junk:0}, blankplug

;----------
; Add the tags XFOCAL,YFOCAL to the structure of object data.
approx_radec_to_xyfocal, stardata1.ra, stardata1.dec, xfocal, yfocal, $
  racen=racen, deccen=deccen, airtemp=airtemp
xydata = replicate(create_struct('XFOCAL', 0L, 'YFOCAL', 0L), $
                   n_elements(stardata1))
xydata.xfocal = xfocal
xydata.yfocal = yfocal
stardata = struct_addtags(stardata1, xydata)

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
   
platescale = 217.7358           ; mm/degree
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
; Add objects that are neither sky nor stds
; First check that geometric conditions are met:
;   no 55'' collisions
;   nothing closer than 68'' to center
;   nothing further than 1.49 deg from center
indx = where(strtrim(stardata.holetype,2) EQ 'OBJECT' $
             AND strtrim(stardata.objtype,2) NE 'SKY', ct)
if(ct ne nsci) then $
  message, 'Not enough non-sky, non-std objects? That must be wrong.'
holeok=geometry_check(stardata[indx], racen, deccen)
iok=where(holeok, nok)
if(nok ne ct) then $
  message, 'Some non-sky, non-std objects not on plate!'
idec=random_decollide(stardata[indx].ra, stardata[indx].dec, seed=seed, $
                      ndec=ndec)
if(ndec ne ct) then $
  message, 'Non-sky, non-std objects collide with one another!'
for i=0L, ct-1L do begin
    addplug = blankplug
    struct_assign, stardata[indx[i]], addplug
    addplug.holetype = 'OBJECT'
    allplug = design_append(allplug, addplug, nadd=nadd1)
    used[indx[i]]=1
endfor
iused=where(used)

;----------
; Add guide fibers
indx=where(strtrim(stardata.holetype,2) EQ 'GUIDE', ct)
if(ct eq 0) then $
  message, 'No guide fibers! Abort!'
holeok=geometry_check(stardata[indx], racen, deccen, /guide)
iok=where(holeok, nok)
if(nok eq 0) then $
  message, 'No guide stars on plate!'
indx=indx[iok]
collok=collision_check(stardata[iused], stardata[indx], /guide)
iok=where(collok, nok)
if(nok eq 0) then $
  message, 'No guide stars do not collide!'
indx=indx[iok]
idec=random_decollide(stardata[indx].ra, stardata[indx].dec, seed=seed, $
                      ndec=ndec, /guide)
guide=bytarr(n_elements(stardata))
guide[indx[idec]]=1

for iguide=0, nguide-1 do begin
;----------
; Assign the nearest available guide fiber(s) to this guide position.
    
    print, 'Assigning guide fiber number ', iguide+1
    
    indx = where(used eq 0 and guide gt 0, ct)
    if (ct EQ 0) then $
      message, 'No guide stars for guide #', iguide
    if (ct GT 0) then begin
        adiff = djs_diff_angle(gfiber[iguide].xprefer, $
                               gfiber[iguide].yprefer, $
                               stardata[indx].xfocal, $
                               stardata[indx].yfocal)
        
        junk = min(adiff, ibest)
        
        addplug = blankplug
        struct_assign, stardata[indx[ibest]], addplug
        addplug.holetype = 'GUIDE'
        addplug.objtype = 'NA'
        addplug.sectarget = 64L
        
        allplug = design_append(allplug, addplug, nadd=nadd1)
        
        used[indx[ibest]] = 1   ; Don't try to target again
    endif
endfor
iused=where(used)

;----------
; Add spectro-photo standards
indx=where(strtrim(stardata.holetype,2) EQ 'OBJECT' AND $
           strtrim(stardata.objtype,2) EQ 'SPECTROPHOTO_STD' AND $
           used EQ 0, ct)
if(ct eq 0) then $
  message, 'No spectro-photo standards! Abort!'
holeok=geometry_check(stardata[indx], racen, deccen)
iok=where(holeok, nok)
if(nok eq 0) then $
  message, 'No spectro-photo standards stars on plate!'
indx=indx[iok]
collok=collision_check(stardata[iused], stardata[indx])
iok=where(collok, nok)
if(nok eq 0) then $
  message, 'No spectro-photo standards do not collide!'
indx=indx[iok]
idec=random_decollide(stardata[indx].ra, stardata[indx].dec, seed=seed, $
                      ndec=ndec)
std=bytarr(n_elements(stardata))
std[indx[idec]]=1

for nsouth=-1, 1, 2 do begin ; First do South, then North half of plate
    nadd = 0L
    while (nadd LT nstd/2.) do begin
        indx = where(std gt 0 AND used eq 0 AND $
                     AND nsouth*stardata.yfocal GE 0, ct)
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
        
        used[indx[ibest]] = 1 ; Don't try to target again
    endwhile
endfor
iused=where(used)

;----------
; Add sky fibers
indx=where(strtrim(stardata.holetype,2) EQ 'OBJECT' AND $
           strtrim(stardata.objtype,2) EQ 'SKY' AND $
           used EQ 0, ct)
if(ct eq 0) then $
  message, 'No sky! Abort!'
holeok=geometry_check(stardata[indx], racen, deccen)
iok=where(holeok, nok)
if(nok eq 0) then $
  message, 'No sky on plate!'
indx=indx[iok]
collok=collision_check(stardata[iused], stardata[indx])
iok=where(collok, nok)
if(nok eq 0) then $
  message, 'No sky do not collide!'
indx=indx[iok]
idec=random_decollide(stardata[indx].ra, stardata[indx].dec, seed=seed, $
                      ndec=ndec)
sky=bytarr(n_elements(stardata))
sky[indx[idec]]=1

while (n_elements(allplug) LT ntot) do begin
    indx = where(sky gt 0 and used eq 0, ct)
    if (ct EQ 0) then $
      message, 'Ran out of sky targets!'
    
    junk = max(priority[indx], ibest)
    addplug = blankplug
    struct_assign, stardata[indx[ibest]], addplug
    addplug.holetype = 'COHERENT_SKY'
    addplug.objtype = 'NA'
    addplug.sectarget = 16L
    
    allplug = design_append(allplug, addplug)
    
    used[indx[ibest]] = 1       ; Don't try to target again
endwhile

;----------
; Write the plPlugMapT file
outhdr = ['completeTileVersion   v1_0', $
          'tileId ' + string(tilenum), $
          'raCen ' + string(racen), $
          'decCen ' + string(deccen) ]

allplug.xfocal = -999         ; These values will be computed by PLATE
allplug.yfocal = -999         ; These values will be computed by PLATE
allplug.spectrographid = -999L ; These values will be computed by PLATE
allplug.fiberid = -999L       ; These values will be computed by PLATE
allplug.throughput = 1L       ; These values will be computed by PLATE
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

;----------
; Run the SDSS "PLATE" code

print
print, 'In the "plate" product run the following commands:"'
print, '   makePlates'
print, '   fiberPlates -skipBrightCheck'
print, '   makeFanuc'
print, '   makeDrillPos'
print, '   use_cs3'
print, '   makePlots -skipBrightCheck'
print
;   setupplate = 'setup plate'
setupplate = 'setup -r /u/schlegel/plate plate' ; ???
spawn, setupplate +'; echo "makePlates" | plate'
spawn, setupplate +'; echo "fiberPlates -skipBrightCheck" | plate'
spawn, setupplate +'; echo "makeFanuc" | plate'
spawn, setupplate +'; echo "makeDrillPos" | plate'
spawn, setupplate +'; echo "use_cs3" | plate'
spawn, setupplate +'; echo "makePlots -skipBrightCheck" | plate'

;----------
; Read the final plPlugMapP file

plugmap = yanny_readone(plugmappfile)
junk = where(strmatch(plugmap.holetype,'GUIDE*'), nguide)
junk = where(strmatch(plugmap.holetype,'COHERENT_SKY*'), ncoherent)
junk = where(strmatch(plugmap.holetype,'OBJECT*'), nobject)
junk = where(strmatch(plugmap.objtype,'SKY*'), nsky)

splog, 'Final NGUIDE = ', nguide
splog, 'Final NCOHERENT_SKY = ', ncoherent
splog, 'Final NSKY = ', nsky
splog, 'Final non-SKY objects = ', nobject - nsky

return
end
;------------------------------------------------------------------------------
