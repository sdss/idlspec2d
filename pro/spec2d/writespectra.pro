
pro writespectra, tt, hdr, filebase

	ntrace = (size(tt))[1]
	npix = (size(tt.counts))[1]
	
	for i=0,ntrace -1 do begin
	  thishdr = hdr
	  name = filebase + string(format = '(i3.3)', i+1) + '.fit'

	  flux = tt[i].counts
	  err = flux*0.0 - 1.0
	  good = where(tt[i].invvar GT 0.0)
	  if (good[0] NE -1) then $
            err(good) = 1.0/sqrt(tt[i].invvar[good])

	  sxaddpar, thishdr, 'OBJID', string(format='(5(i))', $
                tt[i].plugmap.objid)
	  sxaddpar, thishdr, 'MAG', string(format='(5(f8.3))', $
                tt[i].plugmap.mag)
	  sxaddpar, thishdr, 'RA', tt[i].plugmap.ra
	  sxaddpar, thishdr, 'DEC', tt[i].plugmap.dec
	  sxaddpar, thishdr, 'OBJTYPE', tt[i].plugmap.objtype
	  sxaddpar, thishdr, 'PRIMTARG', tt[i].plugmap.primtarget
	  sxaddpar, thishdr, 'SECTARGE', tt[i].plugmap.sectarget
	  sxaddpar, thishdr, 'FIBERID', tt[i].plugmap.fiberId
	  
;
;	Wavelengths
;


	  nparams = n_elements(tt[i].coeff) - 2 
	  sxaddpar, thishdr, 'NWORDER', nparams
	  sxaddpar, thishdr, 'WFITTYPE', tt[i].func
	  sxaddpar, thishdr, 'PIXMIN', tt[i].coeff[0]
	  sxaddpar, thishdr, 'PIXMAX', tt[i].coeff[1]
	  link = string(tt[i].coeff[0]) + string(tt[i].coeff[1])
	 
	  for j=0,nparams-1 do begin
	    keyw = 'COEFF'+string(format='(i1)',j) 
	    sxaddpar, thishdr, keyw, tt[i].coeff[j+2]
	    link = link + string(tt[i].coeff[j+2])
	  endfor

	
	   
	  sxaddpar, thishdr, 'WAT1_001', $
               'wtype=multispec label=Wavelength units=Angstroms'

	  link2 = $
   'wtype=multispec spec1 = "'+string(format = '(i)', i+1)+' 1 1'+ $
    string(tt[i].coeff[2])+string(tt[i].coeff[3]) + $
     string(npix) + ' 0.0 0.0 0.0 0.5 0.0 2 ' + string(nparams) + link

	  length = strlen(link2)
	  place=0
	  for j=0,length-1,60 do begin
	    place = place + 1
	    keyw = 'WAT2_'+string(format='(i3.3)',place) 
	    if (j + 59 LT length) then $
	      sxaddpar, thishdr, keyw, strmid(link2,j,60) $
	    else  sxaddpar, thishdr, keyw, strmid(link2,j)
	  endfor
	    
          sxaddpar, thishdr, 'WAXMAP01', '1 0 0 0 2 0 '

	  writefits, name, [[flux],[err]], thishdr
	endfor

	return
end
	   

	  
	  
