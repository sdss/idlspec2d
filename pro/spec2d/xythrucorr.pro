;+
; NAME:
;   xythrucorr
;
; PURPOSE:
;   Calculate throughput corrections based upon xy offsets of holes
;   relative to calibration reference wavelength
;
; CALLING SEQUENCE:
;   xythrucorr, framefile, [outfilename, outdir]
;
; INPUTS:
;   framefile - spFrame filename
;
; OPTIONAL INPUTS: 
;   outfilename - output file name [default spXYthrucorr-{cam}-{expid}.fits]
;   outdir - output directory to prepend to output file name
;           
; OUTPUTS:
;   Writes spXYfluxcorr-{cam}-{expid}.fits with multiplicative xy input
;       flux corretion factor, matched to the wavelength grid of the
;       input spFrame file.
;
; COMMENTS:
;   This has a dependency upon platedesign.
;   If SEEING50 is 0.0, then the seeing information has not yet been
;       propagated from the guider and the correction will be 1.0 for
;       everything.
;           
; TODO:
;   Some of the calculations are common to a plate and don't need to be
;   redone for each exposure.  This could be refactored to do all exposures
;   for a plate at once.
;   
; BUGS:
;   This will gzip the output file whether you want that or not
;
; PROCEDURES CALLED:
;   plate_center
;   plate_ad2xy
;   ha_fit
;   ha_apply
;   yanny_readone
;   djs_filepath
;   splog
;           
; REVISION HISTORY:
;   22-Oct-2014  Written by Daniel Margala (dmargala@uci.edu), UC Irvine.
;   17-Dec-2014  Adapted by Stephen Bailey, LBL for per-exposure files
;-

pro xythrucorr, framefile, outfilename=outfilename, outdir=outdir, $
    debug=debug

;; Make sure platedesign is setup first
platedesign_dir = getenv('PLATEDESIGN_DIR')
if not keyword_set(platedesign_dir) then begin
    message, 'ERROR - xythrucorr requires platedesign to be setup first'
endif

spframe_read, framefile, hdr=hdr, loglam=loglam
plateid = sxpar(hdr, 'PLATEID')
mjd = sxpar(hdr, 'MJD')
camname = strtrim(sxpar(hdr, 'CAMERAS'))
fwhm = sxpar(hdr, 'SEEING50')

;; Determine output filename
if not keyword_set(outfilename) then begin
    expnum = sxpar(hdr, 'EXPOSURE')
    outfilename = djs_filepath(string(camname, expnum, $
     format='("spXYthrucorr-", a2, "-", i8.8, ".fits")'), root_dir=outdir )
endif else begin
    outfilename = djs_filepath(outfilename, root_dir=outdir)
endelse

;; If guider data hasn't been propagated yet, fwhm will be 0.
;; Don't crash but also don't write a crazy corretion.
if (fwhm EQ 0.0) then begin
    splog, "WARNING - Seeing info not yet in "+framefile
    splog, "WARNING - not generating XY throughput corrections"
    outputcorr = loglam * 0.0 + 1.0
    mwrfits, outputcorr, outfilename, /create
    spawn, ['gzip','-f',outfilename], /noshell
    return
endif

; Calculate the hour angle (ha) at the mid-point of the exposure
ra = sxpar(hdr, 'RADEG')
dec = sxpar(hdr, 'DECDEG')
taimid = 0.5*(sxpar(hdr, 'TAI-BEG') + sxpar(hdr, 'TAI-END'))
jd = 2400000.5D + taimid / (24.D*3600.D)
eq2hor, ra, dec, jd, alt, az, ha, obsname='apo'

; Assume guiding for 5400 Angstroms light
guideon=5400.
; Assume this is pointing number one and no offset
pointing=1L
offset=0L

; Set path to directory for the specified plate number
platedir= plate_dir(plateid)

; Construct path to input plateHoles file, contains fiber positions and other relevant information
holefile= platedir+'/plateHolesSorted-'+strtrim(string(f='(i6.6)',plateid),2)+'.par'
;; splog, 'Opening plateHoles file: '+holefile
check_file_exists, holefile, plateid=plateid

; Parse contents of input file, 
; the input file has a global header and a data entry with multiple fields for each fiber
plateholes = yanny_readone(holefile, hdr=phdr, /anon)
definition= lines2struct(phdr)
default= definition

; trim to BOSS targets (not guide, not light trap, not center hole...)
plateholes = plateholes[where(plateholes.HOLETYPE eq 'BOSS')]

; Sort plateholes into the same fiber order as plugmap
plugmapname = getenv('SPECLOG_DIR') + '/' + strtrim(mjd,2) + '/' $
    + sxpar(hdr, 'PLUGFILE')
plugmap = readplugmap(plugmapname)

spherematch, plugmap.ra, plugmap.dec, $
             plateholes.target_ra, plateholes.target_dec, $
             1./3600, i1, i2, d12
nfiber = n_elements(plateholes.target_ra)
isort = lonarr(nfiber)
isort[i1] = lindgen(nfiber)

plateholes = plateholes[i2[isort]]

; Shorthand: Collect ra/dec/lambda/x/y info for each target
ra= plateholes.target_ra
dec= plateholes.target_dec
lambda= plateholes.lambda_eff   ;; e.g., 5400 for LRGs, 4000 for QSOs
xforig= plateholes.xfocal
yforig= plateholes.yfocal

; Read design hour angle, temperature
; Temperature is set per plate in the platePlans.par file
design_ha=float(strsplit(definition.ha, /extr))
temp=float(definition.temp)

design_platescale_alt=float(definition.design_platescale_alt)
mm_to_arcsec = 3600./design_platescale_alt
fiber_diameter = 2.0 

; Calculate ra and dec for center of the plate
; The values are accessible through 'racen' and 'deccen' after this
; This is necessary for handling multiple pointing plates (and 
; we keep it general for multi-offset plates).
plate_center, definition, default, pointing, offset, $
              racen=racen, deccen=deccen

; Calculate xfocal and yfocal for this pointing (should be similar 
; to xforig/yforig up to round-off)
plate_ad2xy, definition, default, pointing, offset, ra, dec, $
             lambda, xf=xfocal, yf=yfocal, lst=racen+design_ha[pointing-1L], $
             airtemp=temp

; Grab all fibers with with lambda_eff at 5400. These targets are used to fit 
; for guiding parameters: rotation, scale, xshift, yshift.
ifit= where(plateholes.lambda_eff eq guideon, nfit)

;; Calculate xtmp, ytmp at the exposure hour angle
plate_ad2xy, definition, default, pointing, offset, ra, dec, $
            lambda, xf=xtmp, yf=ytmp, lst=racen+ha, $
            airtemp=temp

;; Fit rotation, scale, shift parameters in guide targets
ha_fit, xfocal[ifit], yfocal[ifit], xtmp[ifit], ytmp[ifit], $
       xnew=xtmp2, ynew=ytmp2, rot=rot, scale=scale, $
       xshift=xshift, yshift=yshift

splog, format='(%"rot, scale, xshift, yshift: %f, %f, %f, %f",$)', rot, scale, xshift, yshift

;; Trim to science fibers for this camera
if strmid(camname,1,1) EQ 1 then begin
    ii = where( (plateholes.fiberid GT 0) AND (plateholes.fiberid LE 500), nspec)
endif else begin
    ii = where( (plateholes.fiberid GT 500) AND (plateholes.fiberid LE 1000), nspec)    
endelse
xfocal = xfocal[ii]
yfocal = yfocal[ii]
ra = ra[ii]
dec = dec[ii]

;; Apply rotation, scale, shift adjustments (i4000 targets)
ha_apply, xfocal, yfocal, xnew=xfocal_exp, ynew=yfocal_exp, $
          rot=rot, scale=scale, xshift=xshift, yshift=yshift

;; Calculate xfocal and yfocal at the reference wavelength, where the targets
lambda_ref= replicate(5400., nspec) 

plate_ad2xy, definition, default, pointing, offset, ra, dec, $
             lambda_ref, xf=xfocalref, yf=yfocalref, $
             lst=racen+design_ha[pointing-1L], airtemp=temp

;; Apply rotation, scale, shift adjustments (i4000 targets)
ha_apply, xfocalref, yfocalref, xnew=xfocalref_exp, ynew=yfocalref_exp, $
          rot=rot, scale=scale, xshift=xshift, yshift=yshift

; set up wavelengths to calculate offsets at
nlambda = 71L
lambda_min = 3500.
lambda_max = 10500.
lambda_array = lambda_min+(lambda_max-lambda_min)*(findgen(nlambda)/float(nlambda-1L))

;; Create empty arrays to store guiding corrections
;; and position offsets at each wavelength
tpcorr= fltarr(nspec, nlambda) + 1.0

for i=0L, nlambda-1L do begin
  new_lambda = replicate(lambda_array[i], nspec) 
  ;; Calculate xtmp, ytmp at this wavelength
  plate_ad2xy, definition, default, pointing, offset, ra, dec, $
               new_lambda, xf=xtmp, yf=ytmp, $
               lst=racen+ha, airtemp=temp

  ;; Apply rotation, scale, shift adjustments
  ha_apply, xtmp, ytmp, xnew=xnew, ynew=ynew, rot=rot, scale=scale, $
            xshift=xshift, yshift=yshift
  ;; Calculate throughput for this wavelength at the reference hole position
  ref_offsets = mm_to_arcsec*sqrt((xnew-xfocalref_exp)^2+(ynew-yfocalref_exp)^2)
  tpref = fiberfraction(fwhm, ref_offsets, fiber_diameter)
  ;; Calculate throughput for this wavelength at the actual hole position
  offsets = mm_to_arcsec*sqrt((xnew-xfocal_exp)^2+(ynew-yfocal_exp)^2)
  tp = fiberfraction(fwhm, offsets, fiber_diameter)
  ;; Save throughput correction for this wavelength
  tpcorr[*,i] = tpref/tp
endfor

;; Interpolate to sample wavelengths specific to this exposure
wave = 10^loglam
nwave = (size(loglam, /dim))[0]
outputcorr= fltarr(nwave, nspec)
for i=0L,nspec-1L do begin
    outputcorr[*,i] = interpol(tpcorr[i,*], lambda_array, wave[*,i])
endfor

if keyword_set(debug) then stop

;; Write it
splog, 'Writing '+outfilename
mwrfits, outputcorr, outfilename, /create
spawn, ['gzip','-f',outfilename], /noshell

end
;------------------------------------------------------------------------------
