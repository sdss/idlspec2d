pro quickwave, rawfile, arcname, flatname

   sdssproc, rawfile, arcimg, hdr=archdr, color=color


   camname = strtrim(sxpar(archdr, 'CAMERAS'),2)

   tset = mrdfits(flatname,1)
   traceset2xy, tset, ycen, xcen 
 
   ; boxcar extract
   flux = extract_boxcar(rawfile, xcen, radius = 3.0)

   ; estimate fluxivar
   fluxivar = 1.0/(abs(flux) + 1.0)

   fitarcimage, flux, fluxivar, xpeak, ypeak, wset, aset=aset, $
      lampfile=lampfile, fibermask=tmp_fibmask, bestcorr=bestcorr

