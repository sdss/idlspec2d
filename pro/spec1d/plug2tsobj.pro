;+
; NAME:
;   plug2tsobj
;
; PURPOSE:
;   Construct a tsObj structure that matches all entries in a plugmap structure
;
; CALLING SEQUENCE:
;   tsobj = plug2tsobj(plateid, [ra, dec, plugmap=, dmin= ])
;
; INPUTS:
;   plateid    - Plate number; this can be either a scalar, in which case
;                the same plate is used for all objects, or a vector.
;
; OPTIONAL INPUTS:
;   ra         - Array of right ascension (degrees)
;   dec        - Array of declination (degrees)
;   plugmap    - Plug map structure, which must contain RA, DEC.
;                This must be set if RA and DEC are not set.
;   dmin       - Minimum separation between input position and position
;                of the closest match; default to 60.0 arcsec.
;
; OUTPUTS:
;   tsobj      - tsObj structure, sorted such that each entry corresponds
;                to each entry in the PLUGMAP structure; return 0 if the
;                tsObj file was not found on disk.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The calibObj files are assumed to be in the directory $SPECTRO_DATA/calibobj.
;   These files were to have only the 640 objects for each plate.
;   But since plates can be re-plugged, we must re-sort these
;   files to match the object ordering in the plug-map structure.
;
; EXAMPLES:
;   Read the plug-map for plate 306, fibers 1 to 10, then construct the
;   calibObj structure:
;   > readspec, 306, indgen(10)+1, plug=plug
;   > tsobj = plug2tsobj(306,plugmap=plug)
;
; BUGS:
;   I should match objects based upon the MATCHRA,MATCHDEC in calibObj ???
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   fits_read
;   mrdfits
;   splog
;
; REVISION HISTORY:
;   25-Jun-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function plug2tsobj, plateid, ra, dec, plugmap=plugmap, dmin=dmin

   root_dir = getenv('SPECTRO_DATA')
   if (NOT keyword_set(root_dir)) then $
    message, 'Environment variable SPECTRO_DATA must be set!'

   if (keyword_set(plugmap) $
    AND (NOT keyword_set(ra) OR NOT keyword_set(dec))) then begin
      ra = plugmap.ra
      dec = plugmap.dec
   endif

   if (NOT keyword_set(dmin)) then dmin = 60.0

   ;----------
   ; If PLATEID is a vector, then sort by plate number and call this routine
   ; iteratively, once for each plate number.

   if (n_elements(plateid) GT 1) then begin
      platenums = plateid[ uniq(plateid, sort(plateid)) ]
      nplate = n_elements(platenums)
      for iplate=0, nplate-1 do begin
         indx = where(plateid EQ platenums[iplate])
         tsobj1 = plug2tsobj( platenums[iplate], ra[indx], dec[indx] )

         if (iplate EQ 0) then begin
            tsobj = replicate(tsobj1[0], n_elements(ra))
            tmpobj = tsobj[0]
         endif

         for i=0, n_elements(tsobj1)-1 do begin
            struct_assign, tsobj1[i], tmpobj
            tsobj[indx[i]] = tmpobj
         endfor
      endfor
      return, tsobj
   endif

   platestr = strtrim(string(fix(plateid[0])),2)
;   filename = 'tsObj*-*' + platestr + '.fit*'
   platestr = string(plateid[0],format='(i4.4)')
   filename = 'calibPlateP-' + platestr + '.fits'

   ; Select the first matching file if there are several
;   filename = (findfile(filepath(filename, root_dir=root_dir, $
;    subdirectory='plates')))[0]
   filename = (findfile(filepath(filename, root_dir=root_dir, $
    subdirectory='calibobj')))[0]
   if (NOT keyword_set(filename)) then begin
      print, 'WARNING: calibObj file not found for plate ' + platestr
      return, 0
   endif

   ; Make certain that the file exists and is valid
   message = 0
   fits_read, filename, junk, /no_abort, message=message
   if (keyword_set(message)) then tstemp = 0 $ ; File is invalid FITS file
    else tstemp = mrdfits(filename, 1)
   if (NOT keyword_set(tstemp)) then begin
      print, 'WARNING: calibObj file is empty: ' + filename
      return, 0
   endif

   tsobj1 = tstemp[0]
   struct_assign, {junk:0}, tsobj1 ; Zero-out this new structure
   tsobj = replicate(tsobj1, n_elements(ra))

   ;----------
   ; Find the tsObj-file entry for each plug-map entry by matching
   ; the RA,DEC positions on the sky.  Insist that they agree to 1 arcsec.

   for iplug=0, n_elements(ra)-1 do begin
      ; Assume that this object is non-existent if RA=DEC=0
      if (ra[iplug] NE 0 AND dec[iplug] NE 0) then begin
         adist = djs_diff_angle(tstemp.ra, tstemp.dec, ra[iplug], dec[iplug])
         thismin = min(adist, imin)
         if (thismin GT dmin/3600.) then $
          splog, 'Warning: No matches to within ', dmin, ' arcsec at RA=', $
           ra[iplug], ' DEC=', dec[iplug] $
         else $
          tsobj[iplug] = tstemp[imin]
      endif
   endfor

   return, tsobj
end
;------------------------------------------------------------------------------
