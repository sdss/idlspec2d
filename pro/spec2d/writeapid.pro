pro writeapid, plugmap, filename

	get_lun, ilun
	openw, ilun, filename

	nfibers = (size(plugmap))[1]

	for i=0,nfibers-1 do begin
	   temp = plugmap[i]
	   guess = 1
	   descrip = ' gal'
	   if (temp.objtype EQ 'SKY') then begin
	      descrip = ' sky'
              guess = 0
           endif
	   if (temp.objtype EQ 'SPECTROPHOTO_STD') then begin
              descrip = ' F'
              guess = 2
           endif

;	   printf, ilun, string(format='(i3,i2,a,a7,5(f9.4))', $
;                 i+1,guess,' ',descrip,temp.mag)
	   printf, ilun, string(format='(i3,i2,a,a16,5(f9.4))', $
                 temp.fiberId, guess,' ',temp.objtype, temp.mag)
	endfor

	free_lun, ilun
	return
end

