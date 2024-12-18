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
;   sos_checklimits()
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
function quicktrace, filename, tsetfile, plugmapfile=plugmapfile, nbin=nbin, $
 do_lock=do_lock, fps=fps, plugdir=plugdir, noreject=noreject

   if (NOT keyword_set(nbin)) then nbin = 16

   ;----------
   ; Dispose of any pre-existing tsetfile for this exposure.

   if (keyword_set(findfile(tsetfile))) then rmfile, tsetfile

   ;----------
   ; Read in image
   ;print,filename
   sdssproc, filename, flatimg, flativar, hdr=flathdr, $
    nsatrow=nsatrow, fbadpix=fbadpix, $
    spectrographid=spectrographid, camname=camname, do_lock=do_lock

   if strmatch(string(sxpar(flathdr,'CARTID')), '*FPS-S*', /fold_case) then obs='LCO' else obs='APO'
   
   ; load in configuration parameters from object
   configuration=obj_new('configuration', sxpar(flathdr, 'MJD'), obs)

   ;-----
   ; Decide if this flat is bad

   qbadflat = reject_flat(flatimg, flathdr, nsatrow=nsatrow, fbadpix=fbadpix, noreject=noreject, $
    percent80thresh=configuration->spcalib_reject_calib_percent80thresh())

   if (qbadflat) then begin
      splog, 'ABORT: Unable to reduce flat'
      return, 0
   endif

   ;----------
   ; Read in the plug map file, and sort it
   if (NOT keyword_set(fps)) then begin
     plugmap = readplugmap(plugmapfile, spectrographid, /deredden, /sostags, $
                           hdr=hdrplug, fibermask=fibermask, /plates, mjd=sxpar(flathdr, 'MJD'))
     cartid = long(yanny_par_fc(hdrplug, 'cartridgeId'))
   endif else begin
     if keyword_set(plugmapfile)then begin
        plugmap = readplugmap(plugmapfile, spectrographid, /deredden, /sostags,$
                              hdr=hdrplug, fibermask=fibermask,ccd=camname, $
                              savdir=plugdir, mjd=sxpar(flathdr, 'MJD'))
     endif else fibermask = lonarr(500)
     fibermask[where(fibermask ne 0)] = 0
   endelse

   ;----------
   ; Compute the trace set, but binning every NBIN rows for speed
   ; This is not necessary any more, and it doesn't account for bad columns

   dims = size(flatimg, /dimens)
   ncol = dims[0]
   nrow = dims[1]

; Needed to increase MAXDEV from 0.15 to not drop end fibers on bundles ???
   if (strmid(camname,0,1) EQ 'b') then color = 'blue' $
    else color = 'red'
   ; Set the maxdev to twice what it would be for optimal extraction...
   if (NOT keyword_set(fps)) then begin
        xsol = tracefibercrude(flatimg, flativar, yset=ycen, maxdev=0.30, $
            fibermask=fibermask, cartid=cartid, xerr=xerr, flathdr=flathdr, $
            padding=configuration->spcalib_trace320crude_padding(),/plates, $
            bundlefibers=bundlefibers, nbundle=nbundle)
   endif else begin
        cartid=sxpar(flathdr,'CARTID')
        xsol = tracefibercrude(flatimg, flativar, yset=ycen, maxdev=0.30, $
            fibermask=fibermask, xerr=xerr, flathdr=flathdr,cartid=cartid, $
            padding=configuration->spcalib_trace320crude_padding(), $
            bundlefibers=bundlefibers, nbundle=nbundle)
   endelse

   splog, 'nbundle', nbundle
   nbun = nbundle
   ; Consider a fiber bad only if any of the following mask bits are set,
   ; but specifically not if BADTRACE is set.
   badbits = sdss_flagval('SPPIXMASK','NOPLUG') $
    OR sdss_flagval('SPPIXMASK','BADFLAT')
   ngfiber = total((fibermask AND badbits) EQ 0)
splog, size(xsol, /dimensions)
splog, transpose(xsol[2056,*])
; The following with XERR takes 10X longer than the above, but won't crash ???
     outmask = 0
     ; Ignore values whose central point falls on a bad pixel
     inmask = flativar[xsol,ycen] GT 0
     if (strmid(camname,0,1) EQ 'b') then color = 'blue' else color = 'red'
     xy2traceset, ycen, xsol, tset, $
      ncoeff=configuration->spcalib_xy2traceset_ncoeff(color), $
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
                            ncoeff=5, sigma=sigma, medwidth=medwidth,/quick, $
                            bundlefibers=bundlefibers, nbundle=nbundle)

   if (sos_checklimits('flat', 'XSIGMA', camname, max(medwidth)) $ 
    EQ 'red') then $
    splog, 'WARNING: Median spatial widths = ' $
    + string(medwidth,format='(4f5.2)') + ' pix (LL LR UL UR)'

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
      mwrfits_named, flux, tsetfile, name='FLUX', /create
      mwrfits_named, fluxivar, tsetfile, name='IVAR'
      mwrfits_named, tset, tsetfile, name='TSET'
      mwrfits_named, plugmap, tsetfile,name='PLUGMAP'
      mwrfits_named, fibermask, tsetfile, name='FIBERMASK'
      mwrfits_named, [nbun], tsetfile, name='NUM_BUNDLES'
      mwrfits_named, bundlefibers, tsetfile,name='BUNDLE_FIBERS'
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
