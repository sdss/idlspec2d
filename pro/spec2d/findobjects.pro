function findobjects, plugmap, extras

	base = ['GALAXY', 'QSO', 'STD', 'STAR', 'SKY']

	answer = ptrarr(5)

	for i=0,4 do $
            answer[i] = ptr_new(where(strpos(plugmap.objtype, base[i]) NE -1))

	return,answer
end	
