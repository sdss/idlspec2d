
pro fillplugmap, plugmap, hdr

      tempobjid = intarr(5)
      reads, sxpar(hdr, 'OBJID'), tempobjid
      plugmap.objid = tempobjid

      tempmag = fltarr(5)
      reads, sxpar(hdr, 'MAG'), tempmag
      plugmap.mag = tempmag

      plugmap.objtype = sxpar(hdr, 'OBJTYPE')
      plugmap.ra = sxpar(hdr, 'RA')
      plugmap.dec = sxpar(hdr, 'DEC')
      plugmap.xfocal = sxpar(hdr, 'XFOCAL')
      plugmap.yfocal = sxpar(hdr, 'YFOCAL')
      plugmap.fiberid = sxpar(hdr, 'FIBERID') 
      plugmap.z1d     = sxpar(hdr, 'Z') 
      plugmap.z1d_error = sxpar(hdr, 'Z_ERR') 
      plugmap.z1d_conf  = sxpar(hdr, 'Z_CONF') 
      plugmap.z1d_status = sxpar(hdr, 'Z_STATUS') 
      plugmap.fiberid = sxpar(hdr, 'FIBERID') 
      plugmap.primtarget = sxpar(hdr, 'PRIMTARG')
      plugmap.sectarget = sxpar(hdr, 'SECTARGE')
return
end

pro readidlout, flux, sig=sig, wave=wave, expres=expres, plugmap=plugmap

      if (size(expres,/tname) EQ 'UNDEFINED') then expres = 'spMerge2d*fits'

      files = findfile(expres)
      if (files[0] EQ '') then begin
        print, 'no files found'
        return
      endif
 
      if (ARG_PRESENT(plugmap)) then begin
         shortplugmap =  { $
           OBJID         :  intarr(5), $
           RA            :  0.0d, $  
           DEC           :  0.0d, $ 
           MAG           :  fltarr(5), $
           OBJTYPE       :  ' ', $
           XFOCAL        :  0.0d, $
           YFOCAL        :  0.0d, $
           FIBERID       :  -1, $
           z1d           : 0.0 , $
           z1d_error     : 0.0 , $
           z1d_conf      : 0.0 , $
           z1d_status    :  -9 , $
           PRIMTARGET    :  -1, $
           SECTARGET     :  -1 }
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
      if (ARG_PRESENT(wave)) then wave = fltarr(max(npix),nfiles)
      if (ARG_PRESENT(plugmap)) then $
         plugmap = replicate(shortplugmap,nfiles)


      for i=0,nfiles - 1 do begin
        data = readfits(files[i], hdr, /silent)
       
        flux[0:npix[i]-1,i] = data[*,0]
        if (ARG_PRESENT(sig)) then begin
           if ((size(data))[2] EQ 4) then sig[0:npix[i]-1,i] = data[*,2] $
           else sig[0:npix[i]-1,i] = data[*,1]
        endif

        if (ARG_PRESENT(wave)) then wave[*,i] = $
             findgen(max(npix))*disp[1,i] + disp[0,i]

        if (ARG_PRESENT(plugmap)) then begin
            fillplugmap, shortplugmap, hdr
            plugmap[i] = shortplugmap
        endif
      endfor
      
      return
end
       	
	

