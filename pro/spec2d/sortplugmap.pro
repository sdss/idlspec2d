;+
; NAME:
;   sortplugmap
;
; PURPOSE:
;   Cull and order the plugmap entries to those that match the object fibers
;
; CALLING SEQUENCE:
;  plugsort = sortplugmap(plugmap, [ spectrographid, fibermask=, nfibers= ])
;
; INPUTS:
;   plugmap         - The full plug map structure read in from plPlugMapM
;
; OPTIONAL KEYWORDS:
;   spectrographid  - The spectrograph number, either 1 or 2;
;                     if not set (or 0), then return all object fibers
;   fibermask       - Byte array with bits set for unknown fibers
;   nfibers         - Number of fibers to locate; default 320 if SPECTROGRAPHID
;                     is set, or to 640 otherwise.
;
; OUTPUTS:
;   plugsort        - Culled plugmap structure matching plugged fibers
;
; OPTIONAL OUTPUTS:
;   fibermask       - Modified.
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   fibermask_bits()
;
; REVISION HISTORY:
;   19-Oct-1999  Written by S. Burles, Chicago
;-
;------------------------------------------------------------------------------
function sortplugmap, plugmap, spectrographid, fibermask=fibermask, $
 nfibers=nfibers

   if (n_params() LT 1) then begin
      print, 'Syntax: plugsort = sortplugmap(plugmap, [ spectrographid, $
      print, ' fibermask=, nfibers= ])'
      return, -1
   endif

   if (NOT keyword_set(nfibers)) then begin
      if (keyword_set(spectrographid)) then nfibers = 320 $
       else nfibers = 640
   endif
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfibers) $
    else if (n_elements(fibermask) NE nfibers) then $
     message, 'Number of elements in FIBERMASK do not match NFIBER' 
   
   plugsort = replicate(plugmap[0], nfibers)
   plugsort.holetype = 'OBJECT'
   plugsort.objtype = 'NA'

   if (keyword_set(spectrographid)) then begin
      possible = where( plugmap.holetype EQ 'OBJECT' $
       AND plugmap.fiberid GT nfibers*(spectrographid-1) $
       AND plugmap.fiberid LE nfibers*spectrographid )
   endif else begin
      possible = where( plugmap.holetype EQ 'OBJECT' $
       AND plugmap.fiberid GT 0)
   endelse

   if (possible[0] EQ -1) then $
    message, 'No fibers found in plugmap!'

   nmissing = nfibers - n_elements(possible)
   splog, 'Number of missing fibers: ', nmissing
   
   if (keyword_set(spectrographid)) then begin
      place = plugmap[possible].fiberid - (spectrographid-1)*nfibers - 1
   endif else begin
      place = plugmap[possible].fiberid - 1
   endelse
   plugsort[place] = plugmap[possible]

   ;-------------------------------------------------
   ; Set the appropriate fibermask bit if a fiber not found in plugmap file.
   ; Do this by first setting all bits to 1, then unsetting the good ones.

   fibermask = fibermask OR fibermask_bits('NOPLUG')
   fibermask[place] = fibermask[place] - fibermask_bits('NOPLUG')

   return, plugsort
end
