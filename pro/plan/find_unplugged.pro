;+
; NAME:
;   find_unplugged
;
; PURPOSE:
;   Find all unplugged fibers that still seemed to get a spectrum
;   in the reduced data.
;
; CALLING SEQUENCE:
;   find_unplugged, [ mjd_start ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd_start   - If specified, then only search plates with MJDs beyond
;                 this date.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Identify all spectra in the reduced data that are labelled as NOPLUG
;   in the plug-map file, but still have a wavelength coverage of
;   log(wave) > 0.15 in the reduced spectrum.  These may be fibers that
;   were plugged, but missed by the fiber-mapper.
;
;   A companion procedure to this is PLUGFILE_UNMAPPED.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   platelist
;   readspec
;   splog
;
; REVISION HISTORY:
;   10-Sep-2003  Written by David Schlegel, Princeton (not checked in then)
;-
;------------------------------------------------------------------------------
pro find_unplugged, mjstart

   if (NOT keyword_set(mjstart)) then mjstart = 0
   platelist, plist=plist
   itrim = where(plist.mjd GT mjstart AND strmatch(plist.status1d,'Done*'))
   plist = plist[itrim]

   for i=0L, n_elements(plist)-1 do begin
      readspec, plist[i].plate, mjd=plist[i].mjd, plug=plug, zans=zans, /silent
      ibad = where((zans.anyandmask AND 1) NE 0 $
       AND zans.wcoverage GT 0.15, nbad)
      if (nbad GT 0) then begin
         for j=0L, nbad-1 do $
          splog, plist[i].plate, plist[i].mjd, ibad[j]+1, $
           zans[ibad[j]].wcoverage
      endif
   endfor

   return
end
;------------------------------------------------------------------------------
