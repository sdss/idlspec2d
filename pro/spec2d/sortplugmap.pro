function sortplugmap, plugmap, camcol, nFibers=nFibers

	if (NOT keyword_set(nFibers)) then nFibers=320
	
	if (camcol EQ 1 OR camcol EQ 4) then $	
	   specId = 0 $
	else if (camcol EQ 2 OR camcol EQ 3) then $	
	   specId = 1 $
	else message, 'camcol is not valid'

	plugsort = replicate({plugmapobj}, nFibers)
	plugsort.holetype = 'OBJECT'
	plugsort.objtype = 'NA'

	possible = where(plugmap.holetype EQ 'OBJECT' AND $
		   plugmap.fiberId GT nFibers*specId AND $
		   plugmap.fiberId LE nFibers*(specId+1))

	
	if (possible[0] EQ -1) then $
	  message, 'no fibers found in plugmap!'

	if n_elements(possible) NE nFibers then begin
	  missing = nFibers - n_elements(possible)
	  print, 'SORTPLUGMAP: ', missing, ' fibers seem to be broken'
	endif
	
	place = plugmap[possible].fiberId - specId*nFibers - 1
	plugsort[place] = plugmap[possible]

	return, plugsort
end
	
	
	
