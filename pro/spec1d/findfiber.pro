function findfiber, in, out

	nfiber = n_elements(in)

      
	for i=0, nfiber - 1 do begin

	  tt = where(in[i].objid[0] EQ out.run AND in[i].objid[1] EQ out.rerun AND $
                     in[i].objid[2] EQ out.camcol AND in[i].objid[3] EQ out.field AND $
                     in[i].objid[4] EQ out.id, n)

          if n NE 1 then begin
            tt = 0
            print, i, 'object not found'
          endif

	  if i EQ 0 then flist = tt[0] $
          else flist = [flist, tt[0]]

	endfor

    return, flist
end

	
