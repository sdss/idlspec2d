;+
; NAME:
;   quicktrace
;
; PURPOSE:
;   Trace and boxcar extract an SDSS flat field image
;
; CALLING SEQUENCE:
;   rstruct = quicktrace (filename, tsetfile, plugmapfile, [ nbin=, $
;    /do_lock ] )
;
; INPUTS:
;   filename   - Flat-field filename
;   tsetfile   - Output FITS file
;   plugmapfile- Yanny parameter file with plugmap data, copied to an HDU
;                in the output TSETFILE for use by other routines.
;
; OPTIONAL INPUTS:
;   nbin       - Sub-sampling of row numbers for measuring the spatial
;                profile widths; default to 16 for every 16-th row.
;   do_lock    - Optional keyword for SDSSPROC
;
; OUTPUT:
;   rstruct    - Results to be added html file upon completion
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   apo_checklimits()
;   extract_image
;   fileandpath()
;   findfile()
;   fitflatwidth()
;   mwrfits
;   quickboxcar()
;   readplugmap()
;   reject_flat()
;   rmfile
;   sdssproc
;   sortplugmap()
;   splog
;   trace320crude
;   traceset2xy
;   xy2traceset
;
; REVISION HISTORY:
;   3-Apr-2000  Written by S. Burles & D. Schlegel, APO
;  28-Feb-2002  Modified to do full tracing, speed difference is critical
;  26-Jul-2009  Added keyword for nfibers, KD
;-
;------------------------------------------------------------------------------
function quicktrace, filename, tsetfile, plugmapfile, nbin=nbin, $
 do_lock=do_lock

   if (NOT keyword_set(nbin)) then nbin = 16

   ;----------
   ; Dispose of any pre-existing tsetfile for this exposure.

   if (keyword_set(findfile(tsetfile))) then rmfile, tsetfile

   ;----------
   ; Read in image

   sdssproc, filename, flatimg, flativar, hdr=flathdr, $
    nsatrow=nsatrow, fbadpix=fbadpix, $
    spectrographid=spectrographid, camname=camname, do_lock=do_lock

   ; load in configuration parameters from object
   configuration=obj_new('configuration', sxpar(flathdr, 'MJD'))

   ;-----
   ; Decide if this flat is bad

   qbadflat = reject_flat(flatimg, flathdr, nsatrow=nsatrow, fbadpix=fbadpix, $
    percent80thresh=configuration->spcalib_reject_calib_percent80thresh())

   if (qbadflat) then begin
      splog, 'ABORT: Unable to reduce flat'
      return, 0
   endif

   ;----------
   ; Read in the plug map file, and sort it

   plugmap = readplugmap(plugmapfile, spectrographid, /deredden, /apotags, $
    hdr=hdrplug, fibermask=fibermask)
   cartid = long(yanny_par(hdrplug, 'cartridgeId'))

   ;----------
   ; Compute the trace set, but binning every NBIN rows for speed
   ; This is not necessary any more, and it doesn't account for bad columns

   dims = size(flatimg, /dimens)
   ncol = dims[0]
   nrow = dims[1]

;   if (nrow MOD nbin NE 0) then begin
;      splog, 'ABORT: Unable to bin at ', nbin
;      return, 0
;   endif
;
;   nsmallrow = nrow / nbin
;   smallimg = djs_median(reform(flatimg,ncol,nbin,nsmallrow),2)

; Needed to increase MAXDEV from 0.15 to not drop end fibers on bundles ???
   if (strmid(camname,0,1) EQ 'b') then color = 'blue' $
    else color = 'red'
   ; Set the maxdev to twice what it would be for optimal extraction...
   xsol = trace320crude(flatimg, flativar, yset=ycen, maxdev=0.30, $
    fibermask=fibermask, cartid=cartid, xerr=xerr, flathdr=flathdr, $ 
    padding=configuration->spcalib_trace320crude_padding() ) 
   ; Consider a fiber bad only if any of the following mask bits are set,
   ; but specifically not if BADTRACE is set.
   badbits = sdss_flagval('SPPIXMASK','NOPLUG') $
    OR sdss_flagval('SPPIXMASK','BADFLAT')
   ngfiber = total((fibermask AND badbits) EQ 0)

;   xy2traceset, ycen, xsol, tset, ncoeff=7, maxdev=0.1
; The following with XERR takes 10X longer than the above, but won't crash ???
     outmask = 0
     ; Ignore values whose central point falls on a bad pixel
     inmask = flativar[xsol,ycen] GT 0
     xy2traceset, ycen, xsol, tset, $
      ncoeff=configuration->spcalib_xy2traceset_ncoeff(), $
      maxdev=0.1, outmask=outmask, /double, xerr=xerr, inmask=inmask

   ;----------
   ; Boxcar extract

   flux = quickboxcar(flatimg, flativar, tset=tset, fluxivar=fluxivar)

   ;----------
   ; Optimal-extraction of a sparse number of rows simply to measure the
   ; profile widths, and trigger warnings if the spectrographs look
   ; out-of-focus in the spatial dimension.

   traceset2xy, tset, rownums, xcen
   yrow = lindgen(long(nrow/nbin)) * nbin

   ;   mjd=sxpar(flathdr,'MJD')
   ;   if mjd gt 55055 then yrow=yrow[nrow/nbin/4.:3*nrow/nbin/4.] ;for
   ;   BOSS, helps with xsig by using middle 1/2 of image ; for example
   ;   use only[64:193] ; moving this into fitflatwidth instead
   sigma = 1.0

   extract_image, flatimg, flativar, xcen, sigma, $
    tempflux, tempfluxivar, proftype=1, wfixed=[1,1], yrow=yrow, $
    highrej=5, lowrej=5, npoly=10, ansimage=ansimage, relative=1

   widthset = fitflatwidth(tempflux, tempfluxivar, ansimage, fibermask, $
    ncoeff=5, sigma=sigma, medwidth=medwidth,/quick)

   if (apo_checklimits('flat', 'XSIGMA', camname, max(medwidth)) $ 
    EQ 'red') then $
    splog, 'WARNING: Median spatial widths = ' $
    + string(medwidth,format='(4f5.2)') + ' pix (Left Bottom Top Right)'

   ;----------
   ; Look for Argon lines (or any other emission lines) in the flat-fields,
   ; none of which should be there.

   nocrs = median(flux,3)   ; 3x3 median filter
   noargon = nocrs
   ntrace = (size(nocrs))[2]
   for i=0,ntrace -1 do noargon[*,i] = median(noargon[*,i],15)
   argonlevel = total(nocrs - noargon,1)
   djs_iterstat, argonlevel, median=medargon, sigma=sigargon
   argonsn = medargon / (sigargon /sqrt(ntrace))

   ; ARGONSN = 5 looks to be about normal, argonsn > 15 should throw warning

;   if (argonsn GT 15.0) then $
;    splog,'WARNING: Emission lines (Argon?) in flats at significance=', argonsn

   ;----------
   ; Write traceset to FITS file

   if (sxpar(flathdr,'quality') EQ 'excellent') then begin
      mwrfits, flux, tsetfile, /create
      mwrfits, fluxivar, tsetfile
      mwrfits, tset, tsetfile
      mwrfits, plugmap, tsetfile
      mwrfits, fibermask, tsetfile
   endif else begin
      splog, 'Quality is not excellent - do not write tsetfile'
   endelse

   ;----------
   ; Construct the returned structure

   ; Compute the X position between the central 2 fibers
   ; at their central 2 pixels.
   traceset2xy, tset, xx, yy
   xmid = mean( yy[nrow/2-1:nrow/2,ntrace/2-1:ntrace/2] )

   rstruct = create_struct('TSETFILE', fileandpath(tsetfile), $
                           'NGOODFIBER', float(ngfiber), $
                           'XMID', double(xmid), $
                           'XSIGMA_QUADRANT', float(medwidth), $
                           'XSIGMA', float(max(medwidth)) )

   obj_destroy, configuration

   return, rstruct
end
;------------------------------------------------------------------------------
