pro writeapid, tt, filename

	get_lun, ilun
	openw, ilun, filename

	nfibers = (size(tt))[1]

	for i=0,nfibers-1 do begin
	   temp = tt[i].plugmap
	   guess = 1
	   descrip = ' gal'
	   if (temp.objtype EQ 'SKY') then begin
              guess = 0
	      descrip = ' sky'
           endif
	   if (temp.objtype EQ 'SPECTROPHOTO_STD') then descrip = ' F'

	   printf, ilun, string(format='(i3,i2,a,a7,5(f9.4))', $
                 i+1,guess,' ',descrip,temp.mag)
;	   printf, ilun, string(format='(i3,i2,a,a7,5(f9.4))', $
;                 i+1,guess,' ',temp.objtype,temp.mag)
	endfor

	free_lun, ilun
	return
end

