function readidlout, expres, disp, sig

      if (NOT keyword_set(expres)) then expres = 'idl*fit'

      files = findfile(expres)
      if (files[0] EQ '') then begin
        print, 'no files found'
        return, -1
      endif
 
;
;	Check headers first
;      

      nfiles = n_elements(files)
      npix = lonarr(nfiles)
      disp = dblarr(2,nfiles)

      for i=0,nfiles - 1 do begin
        hdr = headfits(files[i])
        npix[i] = sxpar(hdr, 'NAXIS1')
        disp[0,i] = sxpar(hdr, 'COEFF0')
        disp[1,i] = sxpar(hdr, 'COEFF1')
      endfor

      flux = fltarr(max(npix),nfiles)
      if (ARG_PRESENT(sig)) then sig = fltarr(max(npix),nfiles)

      for i=0,nfiles - 1 do begin
        data = readfits(files[i], /silent)
        flux[0:npix[i]-1,i] = data[*,0]
        if (ARG_PRESENT(sig)) then sig[0:npix[i]-1,i] = data[*,1]
      endfor
      
      return, flux
end
       	
	

