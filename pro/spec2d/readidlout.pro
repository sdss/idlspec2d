
pro fillplugmap, plugmap, hdr

      tempobjid = intarr(5)
      reads, sxpar(hdr, 'OBJID'), tempobjid
      plugmap.objid = tempobjid

      tempmag = fltarr(5)
      reads, sxpar(hdr, 'MAG'), tempmag
      plugmap.mag = tempmag

      plugmap.objtype = sxpar(hdr, 'OBJTYPE')
      plugmap.ra = sxpar(hdr, 'RAOBJ')
      plugmap.dec = sxpar(hdr, 'DECOBJ')
      plugmap.xfocal = sxpar(hdr, 'XFOCAL')
      plugmap.yfocal = sxpar(hdr, 'YFOCAL')
      plugmap.plateid  = sxpar(hdr, 'PLATEID') 
      plugmap.fiberid = sxpar(hdr, 'FIBERID') 
      plugmap.z1d     = sxpar(hdr, 'Z') 
      plugmap.z1d_error = sxpar(hdr, 'Z_ERR') 
      plugmap.z1d_conf  = sxpar(hdr, 'Z_CONF') 
      plugmap.z1d_status = sxpar(hdr, 'Z_STATUS') 
      plugmap.primtarget = sxpar(hdr, 'PRIMTARG')
      plugmap.sectarget = sxpar(hdr, 'SECTARGE')
return
end

pro readidlout, flux, sig=sig, wave=wave, expres=expres, plugmap=plugmap, $
       onlyplug=onlyplug, files=files

      if (size(expres,/tname) EQ 'UNDEFINED') then expres = 'spMerge2d*fits'

      if (NOT keyword_set(files)) then begin
        files = findfile(expres)
        if (files[0] EQ '') then begin
;        print, 'no files found, trying another way'
          check = str_sep(expres, '*')
          ncheck = n_elements(check)
          allfiles = findfile('')
          nfiles = n_elements(allfiles)
          pos = lonarr(ncheck, nfiles)
          for i=0, ncheck -1 do pos[i,*] = strpos(allfiles, check[i])
          if (ncheck GE 2) then $
            diff = (pos[1:ncheck-1,*] - pos[0:ncheck-2,*]) $
                  *(pos[0:ncheck-2,*] NE -1) $
          else diff = pos

          goodfiles = where(total(diff LE 0,1) EQ 0)
          if (goodfiles[0] EQ -1) then begin
            print, 'This way did not work either, returning'
            return
          endif else files = allfiles[goodfiles]
        endif
      endif
 
      if (ARG_PRESENT(plugmap)) then begin
         shortplugmap =  { HDRPLUG,  $
           OBJID         :  intarr(5), $
           RA            :  0.0d, $  
           DEC           :  0.0d, $ 
           MAG           :  fltarr(5), $
           OBJTYPE       :  ' ', $
           XFOCAL        :  0.0d, $
           YFOCAL        :  0.0d, $
           PLATEID       :  -1, $
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

      if (ARG_PRESENT(plugmap)) then $
         plugmap = replicate(shortplugmap,nfiles)

      for i=0,nfiles - 1 do begin
        hdr = sdsshead(files[i])
        npix[i] = sxpar(hdr, 'NAXIS1')
        disp[0,i] = sxpar(hdr, 'COEFF0')
        disp[1,i] = sxpar(hdr, 'COEFF1')
        
        if (keyword_set(onlyplug)) then begin
            fillplugmap, shortplugmap, hdr
            plugmap[i] = shortplugmap
        endif
      endfor

      if (keyword_set(onlyplug)) then return

      flux = fltarr(max(npix),nfiles)
      if (ARG_PRESENT(sig)) then sig = fltarr(max(npix),nfiles)
      if (ARG_PRESENT(wave)) then wave = fltarr(max(npix),nfiles)


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
       	
	

