function quicktrace, filename, tsetfile, plugmapfile, nbin=nbin

   if (NOT keyword_set(nbin)) then nbin=8

   ;----------
   ; Read in image

   sdssproc, filename, flatimg, flativar, hdr=flathdr, $
    nsatrow=nsatrow, fbadpix=fbadpix, $
    spectrographid=spectrographid, camname=camname

   ;-----
   ; Decide if this flat is bad

   qbadflat = reject_flat(flatimg, flathdr, nsatrow=nsatrow, fbadpix=fbadpix)
   if (qbadflat) then begin
      splog, 'ABORT: Unable to reduce flat'
      return, 0
   endif

   ;----------
   ; Read in the plug map file, and sort it

   yanny_read, plugmapfile, pdata
   plugmap = *pdata[0]
   yanny_free, pdata
   plugsort = sortplugmap(plugmap, spectrographid, fibermask=fibermask)

   ;----------
   ; Compute the trace set, but binning every NBIN rows for speed

   dims = size(flatimg, /dimens)
   ncol = dims[0]
   nrow = dims[1]

   if (nrow MOD nbin NE 0) then begin
      splog, 'ABORT: Unable to bin at ', nbin
      return, 0
   endif

   nsmallrow = nrow / nbin
   smallimg = djs_median(reform(flatimg,ncol,nbin,nsmallrow),2)

   xsol = trace320crude(smallimg, yset=ycen, maxshifte=1.40, $
                        fibermask=fibermask)
   ngfiber = total(fibermask EQ 0)

   xy2traceset, ycen*nbin + (nbin-1.0)/2.0, $
    xsol, tset, ncoeff=5, maxdev=0.1, xmin=0, xmax=nrow-1

   traceset2xy, tset, pixnorm, xfit

   ;----------
   ; Boxcar extract

   flux = quickboxcar(flatimg, flativar, tset=tset, fluxivar=fluxivar)

   ;----------
   ; Argon anyone?

   nocrs = median(flux,3)   ; 3x3 median filter
   noargon = nocrs
   ntrace = (size(nocrs))[2]
   for i=0,ntrace -1 do noargon[*,i] = median(noargon[*,i],21)
   argonlevel = total(nocrs - noargon,1)
   djs_iterstat, argonlevel, median=medargon, sigma=sigargon
   argonsn = medargon / (sigargon /sqrt(ntrace))

   ; argonsn = 5 looks to be about normal, argonsn > 10 should throw warning
   if argonsn GT 10.0 then $
      splog,'WARNING: Argon is present, significance is :', argonsn

   ;----------
   ; Write traceset to FITS file

   mwrfits, flux, tsetfile, /create
   mwrfits, fluxivar, tsetfile
   mwrfits, tset, tsetfile
   mwrfits, plugsort, tsetfile
   mwrfits, fibermask, tsetfile

   ;----------
   ; Construct the returned structure

   traceset2xy, tset, xx, yy
   xmin = min(yy)
   xmax = max(yy)
   rstruct = create_struct('NGOODFIBER', ngfiber, $
                           'XMIN', xmin, $
                           'XMAX', xmax )

   return, rstruct
end
