;+
; NAME:
;   quicktrace
;
; PURPOSE:
;   Trace and boxcar extract an SDSS flat field image
;
; CALLING SEQUENCE:
;   rstruct = quicktrace (filename, tsetfile, plugmapfile, nbin=nbin)
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
;   extract_image
;   fileandpath()
;   fitflatwidth()
;   mwrfits
;   quickboxcar()
;   readplugmap()
;   reject_flat()
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
;-
;------------------------------------------------------------------------------
function quicktrace, filename, tsetfile, plugmapfile, nbin=nbin

   if (NOT keyword_set(nbin)) then nbin = 16

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

   plugmap = readplugmap(plugmapfile, /deredden)
   plugsort = sortplugmap(plugmap, spectrographid, fibermask=fibermask)

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

   xsol = trace320crude(flatimg, flativar, yset=ycen, maxdev=0.15, $
                        fibermask=fibermask)
   ngfiber = total(fibermask EQ 0)
   xy2traceset, ycen, xsol, tset, ncoeff=7, maxdev=0.1

   ;----------
   ; Boxcar extract

   flux = quickboxcar(flatimg, flativar, tset=tset, fluxivar=fluxivar)

   ;----------
   ; Optimal-extraction of a sparse number of rows simply to measure the
   ; profile widths, and trigger warnings if the spectrographs look
   ; out-of-focus.

   traceset2xy, tset, rownums, xcen
   yrow = lindgen(long(nrow/nbin)) * nbin
   sigma = 1.0

   extract_image, flatimg, flativar, xcen, sigma, $
    tempflux, tempfluxivar, proftype=1, wfixed=[1,1], yrow=yrow, $
    highrej=5, lowrej=5, npoly=10, ansimage=ansimage, relative=1

   widthset = fitflatwidth(tempflux, tempfluxivar, ansimage, fibermask, $
    ncoeff=5, sigma=sigma)

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

   if (argonsn GT 15.0) then $
    splog,'WARNING: Emission lines (Argon?) in flats at significance=', argonsn

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
   rstruct = create_struct('TSETFILE', fileandpath(tsetfile), $
                           'NGOODFIBER', ngfiber, $
                           'XMIN', xmin, $
                           'XMAX', xmax )

   return, rstruct
end
;------------------------------------------------------------------------------
