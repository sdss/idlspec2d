
pro quicktrace, rawflatfield, tsetfile, plugmapfile, nbin=nbin

;	Read image in

     nin = n_elements(rawflatfield)
     nout = n_elements(tsetfile)

     if (nin NE nout) then begin
       print, 'Number of files IN is different then number of files OUT"
       return
     endif

     if (NOT keyword_set(nbin)) then nbin=8

     for i=0,nin - 1 do begin
       sdssproc, rawflatfield[i], flatimg, spectrographid=spectrographid

       yanny_read, plugmapfile, pdata
       plugmap = *pdata[0]
       yanny_free, pdata
       plugsort = sortplugmap(plugmap, spectrographid, fibermask)

       ncol = (size(flatimg,/dimen))[0]
       nrow = (size(flatimg,/dimen))[1]

       if (nrow mod nbin NE 0) then begin
         print, 'binning does not match'
       endif else begin
         nsmallrow = nrow / nbin
         smallimg = djs_median(reform(flatimg,ncol,nbin,nsmallrow),2)

         xsol = trace320crude(smallimg, yset=ycen, maxshifte=1.40, $
                              fibermask=fibermask)

         xy2traceset, ycen*nbin + (nbin-1.0)/2.0, $
           xsol, tset, ncoeff=5, maxdev=0.1, xmin =0, xmax = nrow-1

         traceset2xy, tset, pixnorm, xfit

;
;	Write traceset to fits file
;

         mwrfits, tset, tsetfile[i], /create
         mwrfits, plugsort, tsetfile[i]
         mwrfits, fibermask, tsetfile[i]

       endelse
     endfor

     return
end
	
	

	
