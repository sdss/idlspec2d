;+
; NAME:
;   sortplugmap
;
; PURPOSE:
;   Cull the nTrace plugmap entries which match the fibers in the
;     current image
;
; CALLING SEQUENCE:
;  plugsort = sortplugmap(plugmap, spectrographid, fibermask, nFibers=nFibers)
;
; INPUTS:
;   plugmap         - The full Plug map structure read in from plPlugMapM
;   spectrographid  - The spectrograph number
;
; OPTIONAL KEYWORDS:
;   nFiebrs         - number of fibers to locate (default 320)
;
; OUTPUTS:
;   plugsort        - culled plugmap structure matching current traces
;
; OPTIONAL OUTPUTS:
;   fibermask       - returns an array of byte with bits set for unknown fibers
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
function sortplugmap, plugmap, spectrographid, fibermask, nFibers=nFibers

   if (NOT keyword_set(nFibers)) then nFibers=320
   
   plugsort = replicate({plugmapobj}, nFibers)
   plugsort.holetype = 'OBJECT'
   plugsort.objtype = 'NA'

   possible = where( plugmap.holetype EQ 'OBJECT' AND $
         plugmap.fiberId GT nFibers*(spectrographid-1) AND $
         plugmap.fiberId LE nFibers*spectrographid )

   if (possible[0] EQ -1) then $
    message, 'No fibers found in plugmap!'

   missing = nFibers - n_elements(possible)
   splog, missing, ' fibers seem to be broken'
   
   place = plugmap[possible].fiberId - (spectrographid-1)*nFibers - 1
   plugsort[place] = plugmap[possible]

   ;-------------------------------------------------
   ;  Fibermask is set to 1 if not found in plugmap file

   fibermask = bytarr(nFibers) OR fibermask_bits('NOPLUG')
 
   ;-------------------------------------------------
   ;  Reset the good plugmap entries to zero, 
   ;  Maybe not the best way to do this

   fibermask[place] = 0    

   return, plugsort
end
