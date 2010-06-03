;+
; NAME:
;   plug2tsobj
;
; PURPOSE:
;   Read photoPosPlate structure of imaging data for a plate
;
; CALLING SEQUENCE:
;   tsobj = plug2tsobj(plateid, [ra, dec, mjd=, indir=, $
;    dmin=, /silent ])
;
; INPUTS:
;   plateid    - Plate number (scalar)
;
; OPTIONAL INPUTS:
;   ra         - Array of right ascension (degrees)
;   dec        - Array of declination (degrees)
;   mjd        - MJD of observation; if set and there is a photoPosPlate file
;                for that plate+MJD, then no coordinates need be specified
;                with RA, DEC
;   indir      - Input directory name for file; default to
;                $BOSS_SPECTRO_REDUX/$RUN2D/$PLATE
;   dmin       - Minimum separation between input position and position
;                of the closest match using MATCHRA,MATCHDEC in the calibObj
;                files; default to 2.0 arcsec.
;
; OUTPUTS:
;   tsobj      - tsObj structure; return 0 if no photoPosPlate files found
;
; OPTIONAL OUTPUTS:
;   silent     - Make the call to MRDFITS silent, but still log warning
;                messages if the calibObj file is not found or is empty.
;
; COMMENTS:
;   Prefer file
;    $BOSS_SPECTRO_REDUX/$RUN2D/$PLATE/$RESOLVE_DIR/photoPosPlate-$PLATE-$MJD.fits
;   If not found, then search for file
;    $BOSS_SPECTRO_REDUX/calibobj/calibPlateP-$PLATE.fits
;   The calibPlateP files must be resorted.
;
;   Print a warning message for non-matched objects only if they are not skies,
;   according to the PLUGMAP structure.
;
; EXAMPLES:
;   Read the plug-map for plate 3690 MJD 55182:
;   > tsobj = plug2tsobj(3690, mjd=55182)
;
; BUGS:
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
function plug2tsobj, plateid, ra, dec, mjd=mjd, indir=indir1, $
 dmin=dmin, silent=silent

   if (n_elements(plateid) NE 1) then $
    message, 'PLATEID must be a scalar!'
   if (keyword_set(mjd)) then begin
      if (n_elements(mjd) NE n_elements(plateid)) then $
       message, 'Number of elements in PLATEID and MJD must agree!'
   endif

   platestr = string(plateid, format='(i4.4)')
   if (n_elements(indir1) GT 0) then indir = indir1 $
    else indir = getenv('BOSS_SPECTRO_REDUX')+'/'+platestr+'/'+strtrim(run2d1,2)

   if (keyword_set(ra)) then begin
      if (n_elements(ra) NE n_elements(dec)) then $
       message, 'Number of elements in RA and DEC must agree!'
   endif

   if (NOT keyword_set(dmin)) then dmin = 2.0

   ;----------
   ; First look for photoPosPlate files (if MJD set)

   if (keyword_set(mjd)) then begin
      mjdstr = string(mjd, format='(i5.5)')
      shortname = 'photoPosPlate-'+platestr+'-'+mjdstr+'.fits'

      ; Look in the output RUN2D directory first, then any subdirectories
      ; if not found
      filename = filepath(shortname, root_dir=indir)
      filename = (findfile(filename, count=ct))[0]
      qsorted = 1B
   endif

   ;----------
   ; Next look for calibPlateP file

   if (keyword_set(filename) EQ 0) then begin
      filename = 'calibPlateP-' + platestr + '.fits'
      filename = (findfile(filepath(filename, root_dir=indir)))[0]
      if (keyword_set(filename)) then begin
         qsorted = 0B
      endif
   end

   ;----------
   ; Read the file

   if (NOT keyword_set(filename)) then begin
      splog, 'WARNING: photoPosPlate file not found for plate ' + platestr
      return, 0
   endif
   splog, 'Reading object calibration file: ' + filename

   ; Make certain that the file exists and is valid
   message = 0
   fits_read, filename, junk, /no_abort, message=message
   if (keyword_set(message)) then tstemp = 0 $ ; File is invalid FITS file
    else tstemp = mrdfits(filename, 1, silent=silent)
   if (NOT keyword_set(tstemp)) then begin
      splog, 'WARNING: calibObj file is empty: ' + filename
      return, 0
   endif

   ;----------
   ; Sort the file if necessary

;   if (keyword_set(ra)) then qsorted = 0B
   if (qsorted) then return, tstemp

   tsobj1 = tstemp[0]
   struct_assign, {junk:0}, tsobj1 ; Zero-out this new structure
   tsobj = replicate(tsobj1, n_elements(ra))

   ;----------
   ; Match objects by positions on the sky

   for iplug=0, n_elements(ra)-1 do begin
      ; Assume that this object is non-existent if RA=DEC=0
      if (ra[iplug] NE 0 AND dec[iplug] NE 0) then begin
         adist = djs_diff_angle(tstemp.matchra, tstemp.matchdec, $
          ra[iplug], dec[iplug])
         thismin = min(adist, imin)
         if (thismin GT dmin/3600.) then begin
            if (keyword_set(plugmap)) then begin
               if (strmatch(plugmap[iplug].objtype,'SKY*') EQ 0) then $
                splog, 'Warning: Unmatched OBJTYPE=' $
                 + strtrim(plugmap[iplug].objtype,2) $
                 + ' at RA=', ra[iplug], ',DEC=', dec[iplug]
            endif else begin
                splog, 'Warning: Unmatched ' $
                 + ' at RA=', ra[iplug], ',DEC=', dec[iplug]
            endelse
         endif else begin
            tsobj[iplug] = tstemp[imin]
         endelse
      endif
   endfor

   return, tsobj
end
;------------------------------------------------------------------------------
