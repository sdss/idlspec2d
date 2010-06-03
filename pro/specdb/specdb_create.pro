;+
; NAME:
;   specdb_create
;
; PURPOSE:
;   Create spectro database.
;
; CALLING SEQUENCE:
;   specdb_create, [ topdir, astrolog= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   topdir     - Search for reduced data files in TOPDIR/MJD/*.
;                This should be an absolute file path, and we default to
;                '/home/data/2d_test'.
;   astrolog   - Search for plug-map files in PLUGDIR/MJD/*.
;                If not specified, then look for these files in the same
;                place as the reduced data.
;
; OUTPUT:
;
; COMMENTS:
;   Look for the input files in:
;     $TOPDIR/MJD/2d/spSpec2d-b1-*.fits    - pppp=plate
;     $SPECLOG/MJD/plPlugMapM-pppp-mmmmm-aa.par - Plug map files, aa=mapper
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   headfits()
;   mrdfits()
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;   spsummary_create()
;
; REVISION HISTORY:
;   27-Mar-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function spsummary_create, ct

   ftemp = create_struct( name='SPSUMMARY', $
    'PLATE', 0L, $
    'FIBERID', 0L, $
    'MJD', 0L, $
    'OBJID', intarr(5), $
    'RA', 0D, $
    'DEC', 0D, $
    'MAG', fltarr(5), $
    'OBJTYPE', '', $
    'PRIMTARGET', 0L, $
    'SECTARGET', 0L, $
    'Z', 0.0 )

   spsummary = replicate(ftemp, ct)

   return, spsummary
end

;------------------------------------------------------------------------------

pro specdb_create, topdir, astrolog=astrolog

   if (NOT keyword_set(topdir)) then begin
      if ((findfile('/usr/sdss/data05/spectro/2d_test'))[0] NE '') then $
       topdir = '/usr/sdss/data05/spectro/rawdata' $
      else if ((findfile('/home/data/2d_test'))[0] NE '') then $
       topdir = '/home/data/2d_test' $
      else if ((findfile('rawdata'))[0] NE '') then $
       topdir = '' $
      else begin
        print, 'Must specify TOPDIR'
        return
      endelse
   endif

   ;----------
   ; Find the list reduced plate numbers

   cd, topdir, current=olddir
   allplates = findfile()

   ;----------
   ; Loop over each reduced plate

   for iplate=0, N_elements(allplates)-1 do begin

      plateid = allplates[iplate]
      platenum = fix(plateid)
;      platestr = string(plateid, format='(i4.4)')

      spfiles = findfile( filepath('spSpec2d-b1-*.fits', $
       root_dir=plateid, subdirectory='2d'), count=nfile )

      ;----------
      ; Loop over each reduced file for this plate

      for ifile=0, nfile-1 do begin
print,spfiles[ifile]

         sp_one = spsummary_create(640)

         ; From HDU=2, get the plug map information
         plugmap  = mrdfits(spfiles[ifile], 2)
         struct_assign, plugmap, sp_one
         sp_one.objtype = strtrim(sp_one.objtype,2) ; Trim trailing spaces

         ; From the header, get the MJD
         hdr = headfits(spfiles[ifile])
         sp_one.mjd = long( sxpar(hdr, 'MJD') )

         ; Random numbers for z for now... ???
         sp_one.z = randomu(1234,640)

         if (keyword_set(sp_all)) then $
          sp_all = [sp_all, sp_one] $
         else $
          sp_all = sp_one

      endfor
   endfor

   cd, olddir
stop

   return
end
;------------------------------------------------------------------------------
