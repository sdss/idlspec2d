function sortplugmap, plugmap, spectrographid, nFibers=nFibers

   if (NOT keyword_set(nFibers)) then nFibers=320
   
   plugsort = replicate({plugmapobj}, nFibers)
   plugsort.holetype = 'OBJECT'
   plugsort.objtype = 'NA'

   possible = where( plugmap.holetype EQ 'OBJECT' AND $
         plugmap.fiberId GT nFibers*(spectrographid-1) AND $
         plugmap.fiberId LE nFibers*spectrographid )

   if (possible[0] EQ -1) then $
    message, 'No fibers found in plugmap!'

   if n_elements(possible) NE nFibers then begin
      missing = nFibers - n_elements(possible)
      print, 'SORTPLUGMAP: ', missing, ' fibers seem to be broken'
   endif
   
   place = plugmap[possible].fiberId - (spectrographid-1)*nFibers - 1
   plugsort[place] = plugmap[possible]

   return, plugsort
end
