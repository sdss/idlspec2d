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
