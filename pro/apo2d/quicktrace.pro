function quicktrace, filename, tsetfile, plugmapfile, nbin=nbin

   if (NOT keyword_set(nbin)) then nbin=8

   ; Read in image
   sdssproc, filename, flatimg, hdr=hdr, spectrographid=spectrographid
   camname = strtrim(sxpar(hdr, 'CAMERAS'),2)

   yanny_read, plugmapfile, pdata
   plugmap = *pdata[0]
   yanny_free, pdata
   plugsort = sortplugmap(plugmap, spectrographid, fibermask=fibermask)

   dims = size(flatimg, /dimens)
   ncol = dims[0]
   nrow = dims[1]

   if (nrow MOD nbin NE 0) then begin
      print, 'Binning does not match'
      rstruct = 0
   endif else begin
      nsmallrow = nrow / nbin
      smallimg = djs_median(reform(flatimg,ncol,nbin,nsmallrow),2)

      xsol = trace320crude(smallimg, yset=ycen, maxshifte=1.40, $
                           fibermask=fibermask)
      ngfiber = total(fibermask EQ 0)

      xy2traceset, ycen*nbin + (nbin-1.0)/2.0, $
       xsol, tset, ncoeff=5, maxdev=0.1, xmin=0, xmax=nrow-1

      traceset2xy, tset, pixnorm, xfit

      ; Write traceset to fits file
      mwrfits, tset, tsetfile, /create
      mwrfits, plugsort, tsetfile
      mwrfits, fibermask, tsetfile

      traceset2xy, tset, xx, yy
      xmin = min(yy)
      xmax = max(yy)
      rstruct = create_struct('FLAVOR', 'flat', $
                              'CAMERA', camname, $
                              'NGOODFIBER', ngfiber, $
                              'XMIN', xmin, $
                              'XMAX', xmax )
   endelse

   return, rstruct
end
