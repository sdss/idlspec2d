;+
; NAME:
;   plug2tsobj
;
; PURPOSE:
;   Construct a tsObj structure that matches all entries in a plugmap structure
;
; CALLING SEQUENCE:
;   tsobj = plug2tsobj(plate, plugmap)
;
; INPUTS:
;   plate      - Plate number.
;   plugmap    - Plug map structure, which must contain at least RA, DEC.
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   tsobj      - tsObj structure, sorted such that each entry corresponds
;                to each entry in the PLUGMAP structure.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   At present, this looks for the tsObj files in '/data/spectro/plates'
;   that were constructed by Fermi to have only the entries for each plate.
;   But since plates can be re-plugged, we must re-sort these files to
;   match the object ordering in the plug-map structure.
;
; EXAMPLES:
;   Read the plug-map for plate 306, fibers 1 to 10, then construct the
;   tsObj structure:
;   > readspec, 306, indgen(10)+1, plug=plug
;   > tsobj = plug2tsobj(306,plug)
;
; BUGS:
;
; PROCEDURES CALLED:
;   mrdfits
;
; REVISION HISTORY:
;   25-Jun-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function plug2tsobj, plate, plugmap

   platestr = strtrim(string(plate),2)
   dirname = '/data/spectro/plates'
   filename = 'tsObj-*-' + platestr + '.fit'

   ; Select the first matching file if there are several
   filename = (findfile(filepath(filename, root_dir=dirname)))[0]
   if (NOT keyword_set(filename)) then $
    message, 'tsObj file not found'

   tstemp = mrdfits(filename, 1, structyp='TSOBJ')

   ;----------
   ; Find the tsObj-file entry for each plug-map entry by matching
   ; the RA,DEC positions on the sky.  Insist that they agree to 1 arcsec.

   cosdec = plugmap.dec * !DPI / 180
   for iplug=0, n_elements(plugmap)-1 do begin
      adist = djs_diff_angle(tstemp.ra, tstemp.dec, $
       plugmap[iplug].ra, plugmap[iplug].dec)
      dmin = min(adist, imin)
      if (dmin GT 1./3600) then $
       message, 'No matches to within 1 arcsec'
      if (iplug EQ 0) then $
       tsobj = tstemp[imin] $
      else $
       tsobj = [tsobj, tstemp[imin]]
   endfor

   return, tsobj
end
;------------------------------------------------------------------------------
