pro quickwave, rawfile, arcname, flatname

   nin = n_elements(rawfile)
   nout = n_elements(arcname)
   nflat = n_elements(flatname)

  
   if (nin NE nout) then begin
       print, 'Number of files IN is different then number of files OUT"
       return
   endif

   if (nin NE nflat) then begin
       print, 'Number of files IN is different then number of FLAT files"
       return
   endif

   for i=0,nin - 1 do begin

     sdssproc, rawfile[i], arcimg, hdr=archdr, color=color

     camname = strtrim(sxpar(archdr, 'CAMERAS'),2)

     tset = mrdfits(flatname[i],1)
     traceset2xy, tset, ycen, xcen 
 
     ; boxcar extract
     flux = extract_boxcar(arcimg, xcen, radius = 3.0)

     ; estimate fluxivar
     fluxivar = 1.0/(abs(flux) + 10.0)

     fitarcimage, flux, fluxivar, xpeak, ypeak, wset, aset=aset, $
        lampfile=lampfile, fibermask=tmp_fibmask, bestcorr=bestcorr, $
        color=color

     mwrfits, wset[i], arcname[i], /create

   endfor
end
