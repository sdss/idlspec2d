;+
; NAME:
;   quickextract
;
; PURPOSE:
;   Science frame extraction with scattered light removal and sky subtraction.
;   S/N is estimated and output for HTML generation.
;
; CALLING SEQUENCE:
;   rstruct = quickextract(tsetfile, wsetfile, fflatfile, rawfile, outsci, $
;    [ radius=, filtsz=, /do_lock ])
;
; INPUTS:
;   tsetfile   - Name of fits file which contains matched trace
;   wsetfile   - Name of fits file which contains matched wavelengths
;   fflatfile  - Name of fits file which containes flat field vectors
;   rawfile    - Name of SDSS science frame to extract
;   outsci     - Name of fits file to store workings of quickextract
;
; OPTIONAL INPUTS:
;   radius     - Radius for boxcar extraction (default 3)
;   filtsz     - Median filter size to apply before average S/N is calculated;
;                default to 25 pixels.
;   do_lock    - Keyword for SDSSPROC
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
;   calcscatimage
;   divideflat
;   djs_mean()
;   djs_median()
;   extract_image
;   extract_boxcar()
;   fibermask_bits()
;   fileandpath()
;   find_whopping()
;   fitflatwidth()
;   fitsn()
;   match_trace()
;   quickboxcar()
;   mrdfits()
;   mwrfits
;   sdssproc
;   skysubtract()
;   splog
;   sxpar()
;   traceset2xy
;
; REVISION HISTORY:
;   3-Apr-2000  Written by S. Burles & D. Schlegel, APO
;-  6-Dec-2014  Included 'splitsky' by Vivek M.
;------------------------------------------------------------------------------
Function quickextract, tsetfile, wsetfile, fflatfile, rawfile, outsci, fullplugfile, outdir, mjd, $
 radius=radius, filtsz=filtsz, splitsky=splitsky, do_lock=do_lock, threshold=threshold,$
 sdssv_sn2=sdssv_sn2, sn2_15=sn2_15, arc2trace=arc2trace, forcea2t=forcea2t
print,'quickextract:',splitsky
   if (n_params() LT 4) then begin
      print, 'Syntax - rstruct = quickextract(tsetfile, wsetfile, fflatfile, $'
      print, ' rawfile, outsci, radius=, filtsz= ])'
      return, 0
   endif

   if (n_elements(tsetfile) NE 1) then return, 0
   if (n_elements(wsetfile) NE 1) then return, 0
   if (n_elements(rawfile) NE 1) then return, 0
   if (NOT keyword_set(radius)) then radius = 3.0
   if (NOT keyword_set(filtsz)) then filtsz = 25

   ;----------
   ; Read in the raw science image

   sdssproc, rawfile, image, invvar, hdr=hdr, camname=camname, $
    nsatrow=nsatrow, fbadpix=fbadpix, do_lock=do_lock
   colorband = strmid(camname,0,1)
   spectroid = strmid(camname,1,1)
   exptime = sxpar(hdr, 'EXPTIME')
   confid = strtrim(sxpar(hdr, 'CONFID'),2)
   fieldid = strtrim(sxpar(hdr,'FIELDID'),2)
   ;----------
   ; Decide if this science exposure is bad
   if not keyword_set(threshold) then threshold=1000.

   qbadsci = reject_science(image, hdr, nsatrow=nsatrow, fbadpix=fbadpix, threshold=threshold)
   if (qbadsci) then begin
      splog, 'ABORT: Unable to reduce science exposure'
      return, 0
   endif

   ;----------
   ; Read in the reduced data from the flat and arc

   tset = mrdfits(tsetfile,2)
   fits_info, tsetfile, n_ext=next
   if next eq 6 then begin
        nbundle = (mrdfits(tsetfile,5))[0]
        bundlefibers = mrdfits(tsetfile,6)
   endif else begin
        nbundle = 25
        bundlefibers = make_array(25,/integer,value=20)
   endelse
   obs = getenv('OBSERVATORY')
   if strmatch(obs, 'LCO', /fold_case) then sp = 2 else sp = 1
   plugsort = readplugmap(fullplugfile ,sp,/deredden, /sostags, fibermask=fibermask, $
                          hdr=pmhdr, savdir=outdir, ccd=camname)

   fflat = mrdfits(fflatfile,0)
   fflatmask = mrdfits(fflatfile,1)
   fibermask = fibermask OR fflatmask
   wset = mrdfits(wsetfile,1)
   traceset2xy, wset, ytemp, logwave

   ;---------------------------------------------------------------------------
   ; Extract and sky-subtract the science image
   ;---------------------------------------------------------------------------

   ;----------
   ; First let's try scattered light fit

   nrow = (size(image,/dimens))[1]
   ncol = (size(image,/dimens))[0]
   skiprow = 8
   yrow = lindgen(nrow/skiprow) * skiprow + skiprow/2

   nfirst = n_elements(yrow)
   wfixed = [1,1] ; Fit gaussian height + width (fixed center position)

   sigma = 1.0
   proftype = 1 ; Gaussian
   npoly=8
   nterms=2

   traceset2xy, tset, ytemp, xcen
   tset_field = (strsplit(tsetfile,'-',/extract))[2]
   
   if keyword_set(forcea2t) then begin
       if (long(tset_field) eq long(fieldid)) then begin
           splog,'Flat for Field Exists but forcing Arc2Trace'
           tset_field = -1
       endif
   endif
   if keyword_set(arc2trace) then begin
       if (long(tset_field) eq long(fieldid)) then begin
            splog, 'Using flat associated with fieldid instead of arc2trace'
       endif else begin
            arcexpid = FILE_BASENAME(wsetfile)
            arcexpid = (strsplit(arcexpid,'-',/extract))[3]
            traceflat = filepath('spTraceTab-'+camname+'-'+arcexpid+'.fits',$
                                root_dir=FILE_DIRNAME(outsci),$
                                subdirectory=['trace',strtrim(mjd,2)])
            traceflat = file_search(traceflat, /fold_case, count=ct)
            if ct gt 0 then begin
                traceflat = traceflat[0]
                xcen = mrdfits(traceflat,0)
                splog, 'Using arc2trace: '+traceflat
            endif
       endelse
   endif
   ;----------
   ; Calculate the shift of the traces between the flat and science exposures

   xnew = match_trace(image, invvar, xcen, /sos)
   if not keyword_set(xnew) then return, 0
   bestlag = median(xnew-xcen)
   minlag = min(xnew-xcen)
   maxlag = max(xnew-xcen)
   splog, 'Match_trace range: ', minlag, bestlag, maxlag
   if (abs(bestlag) GT 1.50) then $
    splog, 'WARNING: Large flexure flat<->science ' $
    + string(bestlag,format='(f5.2)') + ' pix (Post-calibs recommended!)'

   ;----------
   ; Do an optimal extraction for the purpose of measuring scattered
   ; light terms, and for checking the spatial profile widths to see
   ; that the spectrographs are indeed in focus.

   extract_image, image, invvar, xnew, sigma, tempflux, tempfluxivar, $
    proftype=proftype, wfixed=wfixed, yrow=yrow, highrej=5, lowrej=5, $
    npoly=npoly, ansimage=ansimage, relative=1

   ntrace = (size(tempflux,/dimens))[1]
   dims = size(image, /dimens)
   scatfit = calcscatimage(ansimage[ntrace*nterms:*,*], yrow, $
    nscatbkpts=npoly, nx=dims[0], ny=dims[1])
   scatflux = extract_boxcar(scatfit, xcen, radius=radius)

   exptime_factor = (exptime/900.0) > 1.0
   if (colorband EQ 'b') then scatlimit = 20 *exptime_factor $
    else scatlimit = 30 * exptime_factor

   scatmed = fix(median(scatfit))
   scatmax = fix(max(scatfit))

   if (scatmed GT scatlimit) then $
     splog, 'WARNING: Scattered light median = ', scatmed, ' electrons' $
      + ' (Warm CCD or twi?)' $
    else $
     splog, 'Scattered light median = ', scatmed, ' electrons'
   if (scatmax GT 2*scatlimit) then $
     splog, 'WARNING: Scattered light max = ', scatmax, ' electrons' $
      + ' (Warm CCD or twi?)' $
    else $
     splog, 'Scattered light max = ', scatmax, ' electrons'

   ;----------
   ; Check the spatial profile widths, and trigger warning messages
   ; if the spectrographs appear out of focus.

    widthset = fitflatwidth(tempflux, tempfluxivar, ansimage, fibermask, $
                            ncoeff=5, sigma=sigma, medwidth=medwidth, /double, $
                            /quick, bundlefibers=bundlefibers, nbundle=nbundle)
                            
   ; Use the limits as set for the flats, since we don't set limits
   ; for the science exposure widths.  Do not issue a warning message
   ; for smear exposures (which have low S/N), but only science exposures.
   if (sos_checklimits('flat', 'XSIGMA', camname, max(medwidth)) $
    EQ 'red' AND strtrim(sxpar(hdr,'FLAVOR'),2) NE 'smear') then $
    splog, 'WARNING: Median spatial widths = ' $
    + string(medwidth,format='(4f5.2)') + ' pix (LL LR UL UR)'

   ;----------
   ; Boxcar extract - no scattered light correction!

   flux = quickboxcar(image, invvar, tset=tset, fluxivar=fluxivar)
   fluxsub = flux - scatflux
   nfiber = (size(flux,/dimens))[1]

   ;----------
   ; Flat-field

   divideflat, fluxsub, invvar=fluxivar, fflat, /quiet

   ;----------
   ; Check for whopping fibers in SOS reductions, which is
   ; especially useful for flagging affected sky fibers.

   scrunch = djs_median(flux * (fluxivar GT 0), 1) 
   whopping = find_whopping(scrunch, 10000.0, whopct)
   if (whopct GT 0) then begin
      ; Only print the fiber numbers for the first 5 whopping fibers
      wstring = strcompress(string(whopping+1+(spectroid EQ '2')*320))
      splog, 'WARNING: Whopping fiber #' $
       + string(wstring[0:4<(whopct-1)],format='(5a)') $
       + ((whopct GT 5) ? $
        ' (+ ' + strtrim(string(whopct-5),2) + ' more)' : '')
      wp = [whopping - 2 , whopping -1, whopping, whopping+1 , whopping+2]
      wp = (wp > 0) < (nfiber - 1)
      fibermask[wp] = fibermask[wp] OR pixelmask_bits('NEARWHOPPER')
   endif
   iskies = where(strtrim(plugsort.objtype,2) EQ 'SKY' $
      AND (plugsort.fiberid GT 0) AND (fibermask EQ 0), nskies)
   tskys = where(strtrim(plugsort.objtype,2) EQ 'SKY' AND (plugsort.fiberid GT 0), nts)
   splog, 'nplug ', n_elements(plugsort)
   splog, 'tskys ', nts
   tskys = where((fibermask EQ 0), nts)
   splog, 'tskys ', nts
   ;print,nskies
   splog, 'Nskys ',nskies
   ;fibermask[iskies] = 0
   if nskies GT 10 then begin
      skylevel = djs_median(fluxsub[*,iskies], 1)
      outlier = (sort(skylevel))[[0,1,nskies-2,nskies-1]]
      fibermask[iskies[outlier]] = fibermask[iskies[outlier]] $
                                OR fibermask_bits('BRIGHTSKY')
      splog, 'Warning: Rejecting Bright Sky Fibers ', iskies[outlier]
   endif

   ;----------
   ; If too many sky fibers then only choose the first+last per bundle.
   ; Use 20 fibers per bundle, although it really doesn't matter

   iskies = where(strtrim(plugsort.objtype,2) EQ 'SKY' $
    AND (plugsort.fiberid GT 0) AND (fibermask EQ 0), nskies)
splog, 'Nskys ',nskies
;??? modify???
   if (nskies GT 50) then begin
      for i=0, nbundle-1 do begin
         iposs = where(fix(iskies/20) EQ i, nposs)
         if (nposs GT 2) then begin
            ; Trim all but the first + last sky fiber per bundle
            itrim = iskies[iposs[1:nposs-2]]
            splog, 'Trimming extra sky fibers ', itrim
            fibermask[itrim] = fibermask[itrim] $
             OR fibermask_bits('BADSKYFIBER')
         endif
      endfor
      ; Re-select the sky fibers
      iskies = where(strtrim(plugsort.objtype,2) EQ 'SKY' $
       AND (plugsort.fiberid GT 0) AND (fibermask EQ 0), nskies)
splog, 'Nskys ',nskies
   endif
   ;----------
   ; Sky-subtract

   get_tai, hdr, tai_beg, tai_mid, tai_end

   ;skystruct = skysubtract(fluxsub, fluxivar, plugsort, wset, $
   ; objsub, objsubivar, iskies=iskies, fibermask=fibermask, tai=tai_mid, $
   ; sset=sset, npoly=3)

   ;----------
 ;----------
   ; Sky-subtract one final time, this time with dispset (PSF subtraction)
   ; (rejected sky fibers from above remain rejected).
   ; Modify pixelmask in this call.

   nskypoly = 3L
   ;if ~keyword_set(skipsky) then begin
      if keyword_set(splitsky) then begin
; Split sky subtraction between halves of the CCD if requested
; (bolton@utah 2011).
         splog, 'Splitting sky model across spatial CCD halves.'
; Dial down the sptial polynomial order in this case:
         nskypoly = 2L
; Determine where the traces cross the amp break:
         nxfull = (size(image))[1]
         nyfull = (size(image))[2]
         yhw = nyfull / 2
         xhw = nxfull / 2
         isplit = max(where(reform(xnew[yhw,*]) le (xhw + 0.5)))
         splog,'split: ',isplit
	 wset0={FUNC:wset.FUNC,XMIN:wset.XMIN,XMAX:wset.XMAX,COEFF:wset.COEFF(*,0:isplit)}
	 wset1={FUNC:wset.FUNC,XMIN:wset.XMIN,XMAX:wset.XMAX,COEFF:wset.COEFF(*,isplit+1:*)}
	print,'------------------------------------------------------------------'
; Sky subtract both halves:
obs = getenv('OBSERVATORY')
	skystruct0 = skysubtract(fluxsub[*,0:isplit], fluxivar[*,0:isplit], plugsort[0:isplit],wset0, $
             skysub0, skysubivar0, iskies=iskies0,fibermask=fibermask[0:isplit], tai=tai_mid, $
             sset=sset,npoly=nskypoly,obs=obs)
obs = getenv('OBSERVATORY')
         skystruct1 = skysubtract(fluxsub[*,isplit+1:*], fluxivar[*,isplit+1:*], plugsort[isplit+1:*],wset1,  $
             skysub1, skysubivar1, iskies=iskies1,fibermask=fibermask[isplit+1:*], tai=tai_mid, $
             sset=sset,npoly=nskypoly, obs=obs)
; Reassemble outputs for use further below:
         objsub = [[skysub0], [skysub1]]
         objsubivar = [[skysubivar0], [skysubivar1]]
         iskies = [iskies0, iskies1+isplit+1]
; These are needed for QA, although we'll only get the
; low-fiber-number picture in this case:
         skystruct = skystruct0
      endif else begin
obs = getenv('OBSERVATORY')
         skystruct = skysubtract(fluxsub, fluxivar, plugsort,wset,  $
             objsub, objsubivar, iskies=iskies,fibermask=fibermask, tai=tai_mid, $
             sset=sset,npoly=3L, obs=obs)
      endelse
   ;endif else begin
   ;   splog, 'Skipping sky subtraction'
   ;endelse


   ; Issue warnings about very large sky-subtraction chi^2

;   if (keyword_set(skystruct)) then begin
;      thiswave = logwave[*,0]
;      rchi2 = bspline_valu(thiswave, relchi2set)
;      ; Ignore wavelengths near 5577 Ang`
;      i5577 = where(thiswave GT alog10(5577.-10.) $
;       AND thiswave LT alog10(5577.+10.))
;      if (i5577[0] NE -1) then rchi2[i5577] = 0
;      medval = median(rchi2)
;      maxval = max(rchi2, imax)
;      maxwave = 10.^thiswave[imax]
;      if (medval GT 2.) then $
;       splog, 'Warning: Median sky-residual chi2 = ', medval
;      if (maxval GT 60.) then $
;       splog, 'Warning: Max sky-residual chi2 = ', maxval, $
;        ' at' , maxwave, ' Ang (ignoring 5577)'
;   endif

   ;-----------------------------------------------------------------
   ; Analyze spectra for the sky level and signal-to-noise
   ;-----------------------------------------------------------------

   ;----------
   ; Select wavelength range to analyze

   if (colorband EQ 'b') then begin
      icolor = 1
      snfilter = 'g'
      wrange = [4000,5500] ; coverage of g-band
   endif else begin
      icolor = 3
      snfilter = 'i'
      wrange = [6910,8500] ; coverage of i-band
      
      
   endelse
  
  
   ;----------
   ; Compute average (but median-filtered) flux and signal-to-noise
  
   meanflux = fltarr(nfiber)
   meansn = fltarr(nfiber)
   meanobjsub = fltarr(nfiber)
   meanobjsubivar = fltarr(nfiber)
   for ifib=0, nfiber-1 do begin
      ; Select unmasked pixels in the wavelength range for this fiber
      iwave = where(logwave[*,ifib] GT alog10(wrange[0]) $
                AND logwave[*,ifib] LT alog10(wrange[1]) $
                AND objsubivar[*,ifib] GT 0, nwave)
      if (nwave GT 0) then begin
         meanfluxvec = djs_median( flux[iwave,ifib], width=filtsz<nwave, $
          boundary='reflect' )
         meansnvec = djs_median( objsub[iwave,ifib] $
          * sqrt(objsubivar[iwave,ifib]), width=filtsz<nwave, $
          boundary='reflect' )
         meanobjsubvec = djs_median( objsub[iwave,ifib], width=filtsz<nwave, $
          boundary='reflect' )
         meanobjsubivarvec = djs_median( objsubivar[iwave, ifib], width=filtsz<nwave, $
          boundary='reject' ) 
         meanflux[ifib] = djs_mean(meanfluxvec)
         meansn[ifib] = djs_mean(meansnvec)
         meanobjsub[ifib] = djs_mean(meanobjsubvec)
         meanobjsubivar[ifib] = djs_mean(meanobjsubivarvec)
      endif
   endfor

   ;----------
   ; Compute median of the sky fibers

   if (iskies[0] NE -1) then begin
      skylevel = median( meanflux[iskies] )
      if (exptime GT 0) then skylevel = skylevel / float(exptime)
   endif else begin
      skylevel = 0.
   endelse
   ;print,skylevel
 ;----------
   ; Find which fibers are sky fibers + object fibers

   iobj = where(strtrim(plugsort.objtype,2) NE 'SKY' AND plugsort.fiberid GT 0 AND meansn GT 0.2 AND strtrim(plugsort.objtype,2) NE 'NA' AND $
                plugsort.fibermask EQ 0 AND plugsort.FIBER_OFFSET EQ 0 AND plugsort.mag[icolor] GT 0 AND $
                not strmatch(plugsort.FIRSTCARTON, 'bhm_gua_*', /fold_case))
   if (iobj[0] NE -1) then begin
	splog,'#########################     CHECK       ##############################'
        splog,'N-elements--quickextract-fitsn (mag) ', n_elements(plugsort[iobj].mag[icolor]),icolor
       splog,'--quickextract-fitsn (mag) ', plugsort[iobj].mag[icolor]
        splog,'--quickextract-fitsn (sn) ', meansn[iobj]
	splog,'#########################     CHECK       ##############################'
      coeffs = fitsn(plugsort[iobj].mag[icolor], meansn[iobj], $
                    sncode='sos', filter=snfilter, sn2=sn2)
      if keyword_set(sdssv_sn2) then begin
        coeffs2 = fitsn(plugsort[iobj].mag[icolor], meansn[iobj], $
                        sncode='sos2', filter=snfilter, sn2=sn2_v2)
      endif
      if keyword_set(sn2_15) then begin
        ; custom bright SN2 test
        coeffs3 = fitsn(plugsort[iobj].mag[icolor], meansn[iobj], $
                        sncode='sos15', filter=snfilter, sn2=sn2_m_15)
      endif
splog, prelog = ''
splog, sn2
	; Modification by Vivek for RM plates to have original depth of b=10 and r=22
	;Hardcoding plateids
	;if (plugsort[0].plateid eq 7338 or plugsort[0].plateid eq 7339 or plugsort[0].plateid eq 7340) then begin
	;sn2 = 0.8 * sn2
	;endif
   endif else begin
      sn2 = 0
      sn2_v2 = 0
   endelse

   if (keyword_set(rchi2)) then skychi2 = mean(rchi2) $
    else skychi2 = 0.0

   sn2_v2_def = !Values.F_NAN
   sn2_15_def = !Values.F_NAN
   if keyword_set(sdssv_sn2) then begin
    if keyword_set(sn2_v2) then sn2_v2_def = float(sn2_v2)
   endif
   if keyword_set(sn2_15) then begin
    if keyword_set(sn2_m_15) then sn2_15_def = float(sn2_m_15)
   endif

   rstruct = create_struct('SCIFILE', fileandpath(outsci), $
                                'SKYPERSEC', float(skylevel), $
                                'XSIGMA_QUADRANT', float(medwidth), $
                                'XSIGMA', float(max(medwidth)), $
                                'SKYCHI2', float(skychi2), $
                                'FIBERMAG', plugsort.mag[icolor], $
                                'RAWFLUX', float(meanobjsub),$
                                'RAWFLUX_IVAR',float(meanobjsubivar),$
                                'SN2VECTOR', float(meansn^2), $
                                'SN2', float(sn2),$
                                'SEEING', float(sxpar(hdr, 'GSEEING')),$
                                'SN2_v2', sn2_v2_def, $
                                'SN2_15', sn2_15_def)

   ;----------
   ; Write out the extracted spectra

   sxaddpar, hdr, 'FRAMESN2', sn2
   if keyword_set(sdssv_sn2) then sxaddpar, hdr, 'FSN2_v2', sn2_v2_def
   if keyword_set(sn2_15) then sxaddpar, hdr, 'FSN2_15', sn2_15_def
   
   mwrfits_named, objsub, outsci, hdr=hdr, name='OBJSUB', /create
   mwrfits_named, objsubivar, outsci, name='OBJSUBIVAR'
   mwrfits_named, meansn, outsci, name='MEANSN'
   mwrfits_named, sset, outsci, name='SSET'
;   mwrfits_named, relchi2set, outsci, name='RELCHI2'

   return, rstruct
end
;------------------------------------------------------------------------------
