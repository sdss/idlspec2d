;+
; NAME:
;   spcalib
;
; PURPOSE:
;   Extract calibration frames.
;
; CALLING SEQUENCE:
;   spcalib, flatname, arcname, fibermask=, cartid=, $
;    lampfile=, indir=, timesep=, ecalibfile=, plottitle=, $
;    minflat=, maxflat=, arcinfoname=, flatinfoname=, $
;    arcstruct=, flatstruct=, writeflatmodel=, /bbspec]
;
; INPUTS:
;   flatname   - Name(s) of flat-field SDSS image(s)
;   arcname    - Name(s) of arc SDSS image(s)
;   cartid     - Cartridge ID from plugmap
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER].
;                Note this is not modified, but modified copies appear
;                in the returned structures ARCSTRUCT and FLATSTRUCT.
;   lampfile   - Name of file describing arc lamp lines, which would
;                over-ride the default file read by FITARCIMAGE.
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to '.'
;   timesep    - Maximum time separation between flats and arcs to pair them;
;                set to zero to disable this test; default to 7200 sec.
;   ecalibfile - opECalib file to pass to SDSSPROC
;   plottitle  - Prefix for titles in QA plots.
;   minflat    - Parameter for SDSSPROC for pixel flats; default to 0.8
;   maxflat    - Parameter for SDSSPROC for pixel flats; default to 1.2
;   arcinfoname- File name (with path) to output arc extraction and fitting
;                information
;   flatinfoname-File name (with path) to output flat field extraction and
;                fitting information
; writeflatmodel-Set this keyword to write flat data image, ivar, and
;                final extraction model image to a file.  Will only
;                work if "flatinfoname" is present also (ASB).
; writearcmodel- Set this keyword to write arc data image, ivar, and
;                final extraction model image to a file.  Will only
;                work if "arcinfoname" is present also (ASB).
;   bbspec         - use bbspec extraction code
;
; OUTPUTS:
;   arcstruct  - Structure array with extracted arc calibration information
;   flatstruct - Structure array with extracted flat calibration information
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Always pair arcs to the nearest good flat, and flats to the nearest good arc
;   (nearest in time, as defined by the TAI-BEG keyword in the FITS headers).
;
;   Also store SUPERFLATSET from fiberflat, since we need this to remove
;   small scale features present in all spectra

; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   extract_image
;   fiberflat()
;   fitarcimage
;   fitdispersion
;   fitflatwidth()
;   get_tai
;   reject_arc()
;   reject_flat()
;   sdssproc
;   splog
;   trace320crude()
;   traceset2xy
;   xy2traceset
;
; INTERNAL SUPPORT ROUTINES:
;   create_arcstruct()
;   create_flatstruct()
;
; REVISION HISTORY:
;   24-Jan-2000  Written by D. Schlegel, Princeton
;   27-Nov-2000  Changed to proftype 3, added minflat, maxflat keywords
;    8-Jan-2001  And now back to proftype 1, more robust against bad columns
;   26-Jan-2001  And now let's check both 1&3, and use the better fit
;      Apr-2010  Added "write[flat,arc]model" option (A. Bolton, Utah)
;   25-Jan-2011  Added "twophase" test and switching, A. Bolton, Utah
;   29-Mar-2011  Switched to bundle-wise pure IDL extraction, A. Bolton, Utah
;
;-
;------------------------------------------------------------------------------
function create_arcstruct, narc

  ftemp = create_struct( name='ARC_STRUCT', $
    'NAME', '', $
    'TAI', 0D, $
    'TSEP', 0D, $
    'QBAD', 0B, $
    'IFLAT', -1, $
    'BESTCORR', 0.0, $
    'NMATCH', 0L, $
    'MEDWIDTH', fltarr(4), $
    'LAMBDA', ptr_new(), $
    'REJLINE', ptr_new(), $
    'XPEAK', ptr_new(), $
    'XDIF_TSET', ptr_new(), $
    'WSET', ptr_new(), $
    'DISPSET', ptr_new(), $
    'FIBERMASK', ptr_new(), $ 
    'RESLSET', ptr_new(), $
    'MEDRESOL', fltarr(4), $
    'TRACEFLAT', 0,$
    'HDR', ptr_new())

  arcstruct = replicate(ftemp, narc)
  
  return, arcstruct
end
;------------------------------------------------------------------------------
function create_flatstruct, nflat

  ftemp = create_struct( name='FLAT_STRUCT', $
    'NAME', '', $
    'TAI', 0D, $
    'TSEP', 0D, $
    'QBAD', 0, $
    'IARC', -1, $
    'PROFTYPE', 0, $
    'MEDWIDTH', fltarr(4), $
    'FIBERMASK', ptr_new(), $
    'TSET', ptr_new(), $
    'XSOL', ptr_new(), $
    'WIDTHSET', ptr_new(), $
    'FFLAT', ptr_new(), $
    'SUPERFLATSET', ptr_new(), $
    'NBRIGHT', 0, $
    'YMODEL', ptr_new(),$
    'SCATTER', ptr_new(),$
    'FIELDID', 0,$
    'HDR', ptr_new())

  flatstruct = replicate(ftemp, nflat)
  
  return, flatstruct
end
;------------------------------------------------------------------------------

pro spcalib, flatname, arcname, fibermask=fibermask, cartid=cartid, $
             lampfile=lampfile, indir=indir, timesep=timesep, $
             ecalibfile=ecalibfile, plottitle=plottitle, $
             arcinfoname=arcinfoname, flatinfoname=flatinfoname, $
             arcstruct=arcstruct, flatstruct=flatstruct, $
             minflat=minflat, maxflat=maxflat, debug=debug,$
             writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, $
             bbspec=bbspec,plates=plates,legacy=legacy, noreject=noreject, $
             nbundles=nbundles, bundlefibers=bundlefibers, saveraw=saveraw, $
             noarc=noarc, nowrite=nowrite, traceflat=traceflat, $
             force_arc2trace=force_arc2trace, outdir=outdir
    
  if (NOT keyword_set(indir)) then indir = '.'
  if (NOT isa(timesep)) then timesep = 50400
  if (NOT keyword_set(minflat)) then minflat = 0.8
  if (NOT keyword_set(maxflat)) then maxflat = 1.2
  ;timesep = 28800; note coment this line for the final version
  stime1 = systime(1)
 

    if keyword_set(fibermask) then begin
       fibermask_bkup=fibermask
       nt=where(fibermask EQ -100)
       fibermask[*] = 0
       fibermask[nt] = -100
    endif
 
  ;---------------------------------------------------------------------------
  ; Determine spectrograph ID and color from first flat file
  ;---------------------------------------------------------------------------
  sdssproc, flatname[0], indir=indir, $
    spectrographid=spectrographid, color=color
    
  ;---------------------------------------------------------------------------
  ; LOOP THROUGH FLATS + TRACE
  ;---------------------------------------------------------------------------
    
  nflat = N_elements(flatname)
  
  flatstruct = create_flatstruct(nflat)
  
  for iflat=0, nflat-1 do begin
  
    splog, iflat+1, nflat, format='("Tracing flat #",I3," of",I3)'
    
    ;---------------------------------------------------------------------
    ; Read flat-field image
    ;---------------------------------------------------------------------
    
    splog, 'Reading flat ', flatname[iflat]
    sdssproc, flatname[iflat], flatimg, flativar, rdnoiseimg=flatrdnoise, $
      indir=indir, hdr=flathdr, camname=fcamname, outfile=saveraw,$
      nsatrow=nsatrow, fbadpix=fbadpix,$
      ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat, /applycrosstalk
     
    if strmatch(string(sxpar(flathdr,'CARTID')), '*FPS-S*', /fold_case) then obs='LCO' else obs='APO'
    configuration=obj_new('configuration', sxpar(flathdr, 'MJD'), obs)
    
    ;-----
    ; Decide if this flat is bad
    
    qbadflat = reject_flat(flatimg, flathdr, nsatrow=nsatrow, fbadpix=fbadpix, noreject=noreject, $
      percent80thresh=configuration->spcalib_reject_calib_percent80thresh())
    
    if keyword_set(plates) or keyword_set(legacy) then begin
      if (NOT keyword_set(fibermask)) then tmp_fibmask = 0 $
        else tmp_fibmask = fibermask
    endif else begin
      if (NOT keyword_set(fibermask)) then begin
        tmp_fibmask = 0 
      endif else begin
        nt=where(fibermask EQ -100)
        i = 0
        tmp_fibmask = fibermask[nt[i]+1:nt[i+1]-1]
      endelse
    endelse
    
    if (NOT qbadflat) then begin
      ;------------------------------------------------------------------
      ; Create spatial tracing from flat-field image
      ;------------------------------------------------------------------
    
      splog, 'Tracing fibers in ', flatname[iflat]
      xsol = tracefibercrude(flatimg, flativar, yset=ycen, maxdev=1.0, $ ;0.15, $
                             fibermask=tmp_fibmask, cartid=cartid, xerr=xerr, $
                             flathdr=flathdr, plates=plates, $
                             padding=configuration->spcalib_trace320crude_padding(), $
                             plottitle=plottitle+' Traces '+flatname[iflat], $
                             nbundle=nbundles, bundlefibers = bundlefibers, $
                             flatname=flatname[iflat])
 
  
      splog, 'Fitting traces in ', flatname[iflat]
      ntrace = (size(xsol, /dimens))[1]
      outmask = 0
      ; Ignore values whose central point falls on a bad pixel
      ;inmask = flativar[xsol,ycen] GT 0
      ; ASB: New recipe for inmask, just masking fully useless rows,
      ;      since trace320crude has already done clever fill-ins:
      inmask = (total(flativar gt 0., 1) gt 0.) # replicate(1B, ntrace)
      xy2traceset, ycen, xsol, tset, $
       ncoeff=configuration->spcalib_xy2traceset_ncoeff(color), $
       maxdev=0.5, outmask=outmask, /double, xerr=xerr, inmask=inmask

      junk = where(outmask EQ 0, totalreject)
      if (totalreject GT configuration->spcalib_rejecttheshold()) then begin
        splog, 'Reject flat ' + flatname[iflat] + $
          ': ' + string(format='(i8)', totalreject) + ' rejected pixels'
        qbadflat = 1
      endif

      traceset2xy, tset, ycen, xsol

      flatstruct[iflat].tset = ptr_new(tset)
      flatstruct[iflat].xsol = ptr_new(xsol)
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)
      flatstruct[iflat].fieldid = Long(sxpar(flathdr, 'FIELDID'))
      flatstruct[iflat].hdr = ptr_new(flathdr)
    endif else begin
      xsol = 0
      flatstruct[iflat].qbad = 1
      flatstruct[iflat].hdr = ptr_new(flathdr)
    endelse
    
    ;----------
    ; Verify that traces are separated by > 3 pixels
    
    if (qbadflat EQ 0) then begin
      ;sep = xsol[*,1:ntrace-1] - xsol[*,0:ntrace-2]
      ;;- JEB checking trace quality on valid pixels only  (2017-03-03)
      wvalid = where( total(inmask, 2) gt 0.) ;;-- JEB
      sep = xsol[wvalid, 1:ntrace-1] - xsol[wvalid, 0:ntrace-2]
      tooclose = where(sep LT 3)
      if (tooclose[0] NE -1) then begin
        splog, 'Reject flat ' + flatname[iflat] + $
          ': Traces not separated by more than 3 pixels'
; ??? Should reject here!!!
;        qbadflat = 1
      endif
    endif
    
    if (NOT qbadflat) then begin
      ;---------------------------------------------------------------------
      ; Extract the flat-field image to obtain width and flux
      ;---------------------------------------------------------------------

      sigma = configuration->spcalib_sigmaguess() ; Initial guess for gaussian width
      highrej = 15
      lowrej = 15
      npoly = 10 ; Fit 1 terms to background
      wfixed = [1,1] ; Fit the first gaussian term + gaussian width
 
      ; SDSS-I better with proftype=3 for exponential cubic
      ; BOSS better with proftype=1
      proftype = configuration->spcalib_extract_image_proftype()    ; |x|^3
      splog, 'Extracting flat with proftype=', proftype
        outname = string(format='(a,i8.8,a)',flatinfoname, sxpar(flathdr, 'EXPOSURE'), '.fits')
	       outname = repstr(outname, 'spFlat', 'spFrame')

      extract_image, flatimg, flativar, xsol, sigma, flux, fluxivar, $
       proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
       npoly=npoly, relative=1, ansimage=ansimage, reject=[0.1, 0.6, 0.6], $
       chisq=chisq3, outname =outname, plottitle=plottitle+' Flat Extraction profile for '+flatname[iflat], $
       debug=debug

      widthset3 = fitflatwidth(flux, fluxivar, ansimage, tmp_fibmask, $
                               ncoeff=configuration->spcalib_fitflatwidth_ncoeff(), $
                               sigma=sigma, medwidth=medwidth, $
                               mask=configuration->spcalib_fitflatwidth_mask(flux,fluxivar), $
                               inmask=configuration->spcalib_fitflatwidth_inmask(flux,fluxivar), $
                               /double, nbundles=nbundles, bundlefibers = bundlefibers)
      ansimage = 0

      widthset = widthset3
      splog, 'Using proftype=', proftype
      
      junk = where(flux GT 1.0e5, nbright)
      splog, 'Found ', nbright, ' bright pixels in extracted flat ', $
        flatname[iflat], format='(a,i7,a,a)'
        
      flatstruct[iflat].proftype  = proftype
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)
      flatstruct[iflat].widthset = ptr_new(widthset)
      flatstruct[iflat].medwidth  = medwidth
      
    endif
    
    flatstruct[iflat].name = flatname[iflat]
    get_tai, flathdr, tai_beg, tai_mid, tai_end
    flatstruct[iflat].tai = tai_mid
    flatstruct[iflat].qbad = qbadflat
    obj_destroy,configuration
  endfor
  
  ;---------------------------------------------------------------------------
  ; LOOP THROUGH ARCS + FIND WAVELENGTH SOLUTIONS
  ;---------------------------------------------------------------------------
  
  narc = N_elements(arcname)
  
  arcstruct = create_arcstruct(narc)
  for iarc=0, narc-1 do begin
    ccd = strtrim(sxpar(flathdr, 'CAMERAS'),2)
    arcid = FILE_BASENAME(arcname[iarc])
    arcid = (strsplit((strsplit(arcid,'-',/extract))[2],'.',/extract))[0]
    traceflat = filepath('spTraceTab-'+ccd+'-'+arcid+'.fits',root_dir='.',$
                         subdirectory=['..','trace',strtrim(sxpar(flathdr, 'MJD'),2)])
    traceflat = file_search(traceflat, /fold_case, count=ct)
    if ct gt 0 then begin
        traceflat = traceflat[0]
        traceflat_xsol = ptr_new(mrdfits(traceflat,0))
    endif else traceflat = 0


    noarc=0
    splog, iarc+1, narc, format='("Extracting arc #",I3," of",I3)'
    
    ;---------------------------------------------------------------------
    ; Read the arc
    ;---------------------------------------------------------------------
    
    splog, 'Reading arc ', arcname[iarc]
    
    sdssproc, arcname[iarc], arcimg, arcivar, rdnoiseimg=arcrdnoise, $
      indir=indir, hdr=archdr, outfile=saveraw,$
      nsatrow=nsatrow, fbadpix=fbadpix, $
      ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat,/applycrosstalk
    ny = (size(arcimg,/dimens))[1]
     
    cart = sxpar(archdr,'CARTID')
    if strmatch(cart, '*FPS-S*', /fold_case) then lco=1 else lco = 0
    configuration=obj_new('configuration', sxpar(archdr, 'MJD'), obs)
    
    splog, 'Fraction of bad pixels in arc = ', fbadpix
    
    ;----------
    ; Decide if this arc is bad
    
    qbadarc = reject_arc(arcimg, archdr, nsatrow=nsatrow, fbadpix=fbadpix, noreject=noreject)
    
    ;----------
    ; Identify the nearest flat-field for this arc, which must be
    ; within TIMESEP seconds and be a good flat.
    
    get_tai, archdr, tai_beg, tai_mid, tai_end
    tai = tai_mid
    
    iflat = -1
    igood = where(flatstruct.qbad EQ 0)
    if (igood[0] NE -1) then begin
      tsep = min( abs(tai - flatstruct[igood].tai), ii )
      if (tsep LE timesep AND timesep NE 0) then iflat = igood[ii] $
      else begin
        if timesep eq 0 then iflat = igood[ii]
      endelse
    endif
    
    if (iflat GE 0) then begin
      splog, 'Arc ' + arcname[iarc] + ' paired with flat ' + flatname[iflat]
    endif else begin
      splog, 'Arc ' + arcname[iarc] + ' paired with no flat'
      qbadarc = 1
      noarc=0
    endelse
    
    if (NOT qbadarc) then begin
      if Long(sxpar(archdr, 'FIELDID')) eq flatstruct[iflat].fieldid then begin
          if not keyword_set(force_arc2trace) then begin
            traceflat = 0
          endif else begin
            if keyword_set(traceflat) then $
              splog, 'Matching Flat and spTraceTab exists, overriding Flat with spTraceTab'
          endelse
      endif
      if keyword_set(traceflat) then begin
        splog, 'Using adjusted xsol from ',traceflat
        xsol = *(traceflat_xsol)
      endif else begin
        xsol = *(flatstruct[iflat].xsol)
      endelse
      widthset = *(flatstruct[iflat].widthset)
      tmp_fibmask = *(flatstruct[iflat].fibermask)
      proftype = flatstruct[iflat].proftype
      
      ;----------
      ; Calculate possible shift between arc and flat
      
      xcor = match_trace(arcimg, arcivar, xsol)
      
      bestlag = median(xcor-xsol)
      if (abs(bestlag) GT 2.0) then begin
        qbadarc = 1
        splog, 'Reject arc: pixel shift is larger than 2 pixel'
        splog, 'Reject arc ' + arcname[iarc] + ': Pixel shift = ', bestlag
      endif
    endif
    
    if (NOT qbadarc) then begin
      splog, 'Shifting traces with match_trace', bestlag
      ;         splog, 'Shifting traces to fit arc by pixel shift of ', bestlag
      
      ;---------------------------------------------------------------------
      ; Extract the arc image
      ;---------------------------------------------------------------------
      
      traceset2xy, widthset, xx, sigma2
      
      highrej = 100 ; JG (some trouble with bad trace match at earlier step)
      lowrej  = 100 ; JG

      wfixed = [1,0] ; ASB: Don't fit for width terms.
      ;;;;;;; debug;;;;;;;;;
      print, '-------------------------'
      print, nbundles
      print, '-------------------------'
      print, bundlefibers
      print, '-------------------------'
      print, max(xcor-xsol), min(xcor-xsol)
      print, max(sigma2), min(sigma2), median(sigma2)
      print, '-------------------------'

      splog, 'Extracting arc'
      pixelmask=lonarr(size(flux,/dimens)) ; JG : add a mask
      extract_bundle_image, arcimg, arcivar, arcrdnoise, xcor, sigma2, $
                            flux, fluxivar, proftype=proftype, wfixed=wfixed, $
                            highrej=highrej, lowrej=lowrej, npoly=0L, relative=1, $
                            reject=[0.1, 0.6, 0.6], ymodel=ymodel, $
                            buffsize=8L, pixelmask=pixelmask, debug=debug, $
                            use_image_ivar=1, $ ; JG more robust to trace offsets
                            nbundles=nbundles, bundlefibers=bundlefibers

      if keyword_set(debug) then begin
        flatextfile = string(format='(a,i8.8,a)',flatinfoname, sxpar(flathdr, 'EXPOSURE'), '.fits')
        arcextfile = repstr(repstr(arcname[iarc], 'sdR', 'spArcFlux'),'.fit','.fits')

        mwrfits_named, flux,     arcextfile, name='FLUX',/create
        mwrfits_named, fluxivar, arcextfile, name='IVAR'
      endif


      ;JG debug
      ;outname="arc.fits"
      ;sxaddpar, bighdr, 'BUNIT', 'electrons/row'
      ;mwrfits_named, flux, outname, hdr=bighdr, name='FLUX','/create
      ;sxaddpar, hdrfloat, 'BUNIT', 'electrons/row'
      ;mwrfits_named, fluxivar, outname, hdr=hdrfloat, name='IVAR', desc=' Inverse variance'
      ;mwrfits_named, pixelmask, outname, hdr=hdrfloat, name='MASK', desc=' MASK'
      ;STOP


      ; flag to determine whether or not to do 2-phase arc solution:
      twophase = sxpar(archdr, 'TWOPHASE')
      if keyword_set(twophase) then splog, 'Setting 2-phase readout flag'

      ;---------------------------------------------------------------------
      ; Compute correlation coefficient for this arc image
      ;---------------------------------------------------------------------

      splog, 'Searching for wavelength solution'
      aset = 0
      fitarcimage, flux, fluxivar, aset=aset, color=color, lco=lco,$
       lampfile=lampfile, fibermask=tmp_fibmask, bestcorr=bestcorr, $
       acoeff=configuration->spcalib_arcfitguess_acoeff(color), $
       dcoeff=configuration->spcalib_arcfitguess_dcoeff(color), $
       wrange=configuration->spcalib_fitarcimage_wrange(color), $
       twophase=twophase, nbundle=nbundles, bundlefibers = bundlefibers

      arcstruct[iarc].bestcorr = bestcorr
      
      if ((color EQ 'blue' AND bestcorr LT 0.5) $
        OR (color EQ 'red'  AND bestcorr LT 0.5) ) then begin
        qbadarc = 1
        splog, 'Reject arc ' + arcname[iarc] + $
          ': correlation is only = ' + string(format='(i4)', bestcorr)
      endif
    endif
    
    if (NOT qbadarc) then begin
    
      ;---------------------------------------------------------------------
      ; Compute wavelength calibration
      ;---------------------------------------------------------------------
    
      arccoeff = configuration->spcalib_arccoeff()
      
      splog, 'Searching for wavelength solution'
      fitarcimage, flux, fluxivar, xpeak, ypeak, wset, ncoeff=arccoeff, $
       aset=aset, color=color, lampfile=lampfile, fibermask=tmp_fibmask, $
       lambda=lambda, rejline=rejline, xdif_tset=xdif_tset, lco=lco, $
       acoeff=configuration->spcalib_arcfitguess_acoeff(color), $
       dcoeff=configuration->spcalib_arcfitguess_dcoeff(color), $
       wrange=configuration->spcalib_fitarcimage_wrange(color), $
       twophase=twophase, nbundle=nbundles, bundlefibers = bundlefibers

      if (NOT keyword_set(wset)) then begin
        splog, 'Wavelength solution failed'
        qbadarc = 1
      endif else begin
        nfitcoeff = configuration->spcalib_ncoeff(color)
        ilamp = where(rejline EQ '')
        if keyword_set(debug) then $
            arc_test_file = repstr(string(format='(a,i8.8,a)',arcinfoname, $
                            sxpar(archdr, 'EXPOSURE'), '.fits'),'spArc','fitdispersion')

        dispset = fitdispersion(flux, fluxivar, xpeak[*,ilamp], $
          sigma=configuration->spcalib_sigmaguess(), ncoeff=nfitcoeff, $
          xmin=0.0, xmax=ny-1, bundlefibers = bundlefibers,$
          medwidth=wsigarr, numbundles=nbundles, width_final=width_final,$
          arc_test_file = arc_test_file)
          
        if keyword_set(debug) then $
            arc_test_file = repstr(string(format='(a,i8.8,a)',arcinfoname, $
                            sxpar(archdr, 'EXPOSURE'), '.fits'),'spArc','fitspectraresol')
        reslset = fitspectraresol(flux,fluxivar, xpeak[*,ilamp], wset,  $
          ncoeff=nfitcoeff, xmin=0.0, xmax=ny-1, bundlefibers = bundlefibers, $
          medresol=sresarr, numbundles=nbundles, resol_final=resol_final, $
          arc_test_file=arc_test_file, waves = lambda[ilamp])

        arcstruct[iarc].dispset = ptr_new(dispset)
        arcstruct[iarc].wset = ptr_new(wset)
        arcstruct[iarc].nmatch = N_elements(lambda)
        arcstruct[iarc].lambda = ptr_new(lambda)
        arcstruct[iarc].rejline = ptr_new(rejline)
        arcstruct[iarc].tsep = tsep
        arcstruct[iarc].xpeak = ptr_new(xpeak)
        arcstruct[iarc].xdif_tset = ptr_new(xdif_tset)
        arcstruct[iarc].fibermask = ptr_new(tmp_fibmask)
        arcstruct[iarc].medwidth = wsigarr
        arcstruct[iarc].medresol = sresarr
        arcstruct[iarc].traceflat = traceflat
        arcstruct[iarc].reslset = ptr_new(reslset)
        ;------------------------------------------------------------------
        ; Write information on arc lamp processing
        
        if (keyword_set(arcinfoname)) then begin
           write_sparc, arcinfoname, iarc, arcstruct, archdr, $
             flatname[iflat], arcname, fbadpix, bestcorr, tai, $
             lambda, xpeak, ntrace, width_final, ilamp, $
             arcimg, arcivar, ymodel, nowrite=nowrite, $
             writearcmodel=writearcmodel
        
          ymodel = 0
        endif
        
      endelse
      
    endif
    
    arcstruct[iarc].name = arcname[iarc]
    arcstruct[iarc].tai = tai
    arcstruct[iarc].iflat = iflat
    arcstruct[iarc].qbad = qbadarc
    
    obj_destroy,configuration
  endfor
  
  arcimg = 0
  arcivar = 0
  
  ;---------------------------------------------------------------------------
  ; LOOP THROUGH FLATS + CREATE FIBERFLATS
  ;---------------------------------------------------------------------------
  
  for iflat=0, nflat-1 do begin
  
    splog, iflat+1, nflat, $
      format='("Create fiberflats for flat #",I3," of",I3)'
      
    ;----------
    ; Identify the nearest arc for each flat-field, which must be
    ; within TIMESEP seconds and be good.
      
    iarc = -1
    igood = where(arcstruct.qbad EQ 0)
    if (igood[0] NE -1) then begin
      tsep = min( abs(flatstruct[iflat].tai - arcstruct[igood].tai), ii )
      if (tsep LE timesep AND timesep NE 0) then iarc = igood[ii] $
      else begin
        if timesep eq 0 then iarc=igood[ii]
      endelse
      flatstruct[iflat].tsep = tsep
    endif
    
    if (iarc GE 0) then begin
      splog, 'Flat ' + flatname[iflat] + ' paired with arc ' + arcname[iarc]
    endif else begin
      splog, 'Flat ' + flatname[iflat] + ' paired with no arc'
      flatstruct[iflat].qbad = 1 ; Flat is bad if no companion arc exists
      noarc = 1
    endelse
    
    flatstruct[iflat].iarc = iarc
    
    if (NOT flatstruct[iflat].qbad) then begin
    
      widthset = *(flatstruct[iflat].widthset)
      wset = *(arcstruct[iarc].wset)
      xsol = *(flatstruct[iflat].xsol)
      traceflat = arcstruct[iarc].traceflat
      tmp_fibmask = *(flatstruct[iflat].fibermask)
      proftype = flatstruct[iflat].proftype
      
      ;---------------------------------------------------------------------
      ; Read flat-field image (again)
      ;---------------------------------------------------------------------
      
      ; If there is only 1 flat image, then it's still in memory
      if (nflat GT 1) then begin
        splog, 'Reading flat ', flatname[iflat]
        sdssproc, flatname[iflat], flatimg, flativar, $
          indir=indir, hdr=flathdr, $
          ecalibfile=ecalibfile, $
          minflat=minflat, maxflat=maxflat,/applycrosstalk
      endif
      configuration=obj_new('configuration',sxpar(flathdr, 'MJD'), obs)

      ;---------------------------------------------------------------------
      ; Extract the flat-field image
      ;---------------------------------------------------------------------
      
      traceset2xy, widthset, xx, sigma2   ; sigma2 is real width
      highrej = 15
      lowrej = 15
      npoly = 5 ; Fit 5 terms to background, just get best model
;      wfixed = [1,1] ; Fit gaussian plus both derivatives
      wfixed = [1,0] ; Do not refit for Gaussian widths, only flux ???

        outname = string(format='(a,i8.8,a)',flatinfoname, sxpar(flathdr, 'EXPOSURE'), '.fits')
       outname = repstr(outname, 'spFlat', 'spFrame')
      extract_bundle_image, flatimg, flativar, flatrdnoise, xsol, sigma2, $
                            flux, fluxivar, proftype=proftype, wfixed=wfixed, $
                            highrej=highrej, lowrej=lowrej, npoly=0L, $
                            relative=1, chisq=schisq, ansimage=ansimage2, $
                            reject=[0.1, 0.6, 0.6], ymodel=ymodel, outname=outname,$
                            buffsize=8L, nbundles=nbundles, bundlefibers=bundlefibers, $
                            debug=debug

      if keyword_set(debug) then begin
        flatextfile = string(format='(a,i8.8,a)',flatinfoname, sxpar(flathdr, 'EXPOSURE'), '.fits')
        flatextfile = repstr(flatextfile, 'spFlat', 'spFlatFlux')

        mwrfits_named, flux,     flatextfile, hdr=flathdr, name='FLUX' /create
        mwrfits_named, fluxivar, flatextfile, name='IVAR'
      endif

      if (keyword_set(bbspec)) then begin
         basisfile = 'spBasisPSF-*-'+strmid(arcstruct[iarc].name,4,11)+'.fits'
         tmproot = 'tmp-'+strmid(flatstruct[iflat].name,4,11)
         bbspec_extract, flatimg, flativar, bbflux, bbfluxivar, $
              basisfile=basisfile, ximg=xsol, ymodel=bb_ymodel, $
              tmproot=tmproot, nbundles=nbundles, bundlefibers=bundlefibers, $
              /batch ; ??? set batch

         ; Deal with case of only the first few spectra being re-extracted...
         dims = size(bbflux,/dimens)
         flux[0:dims[0]-1,0:dims[1]-1] = bbflux
         fluxivar[0:dims[0]-1,0:dims[1]-1] = bbfluxivar $
          * (fluxivar[0:dims[0]-1,0:dims[1]-1] GT 0) ; <- Retain old rejection

         outfile = 'ymodel-'+strmid(flatstruct[iflat].name,4,11)+'.fits'
         mwrfits_named, bb_ymodel, outfile, name='BB_YMODEL', /create
         mwrfits_named, ymodel, outfile, name='YMODEL'
      endif

;x      splog, 'First  extraction chi^2 ', minmax(fchisq)
      splog, 'Second extraction chi^2 ', minmax(schisq)
      
      xaxis = lindgen(n_elements(schisq)) + 1
      djs_plot, xaxis, schisq, $
        xrange=[0,N_elements(schisq)], xstyle=1, $
;x        yrange=[0,max([max(fchisq), max(schisq)])], $
        yrange=[0,max(schisq)], $
        xtitle='Row number',  ytitle = '\chi^2', $
        title=plottitle+' flat extraction chi^2 for '+flatname[iflat]
        
      djs_oplot, !x.crange, [1,1]
;x      djs_oplot, xaxis, fchisq, color='green'
      
      xyouts, 100, 0.05*!y.crange[0]+0.95*!y.crange[1], $
        'BLACK = Final chisq extraction'
;x      xyouts, 100, 0.08*!y.crange[0]+0.89*!y.crange[1], $
;x        'GREEN = Initial chisq extraction'
        
      ;---------------------------------------------------------------------
      ; Compute fiber-to-fiber flat-field variations
      ;---------------------------------------------------------------------
        
      sigma2 = 0
      xsol = 0
      
      fflat = fiberflat(flux, fluxivar, wset, fibermask=tmp_fibmask, $
        /dospline, pixspace=5, $
        plottitle=plottitle+' Superflat '+flatstruct[iflat].name, $
        superflatset=superflatset, $
        badflatfracthresh=configuration->spcalib_fiberflat_badflatfracthresh(),$
        minval=configuration->spcalib_fiberflat_minval(flux))
        
      if (n_elements(fflat) EQ 1) then begin
        flatstruct[iflat].qbad  = 1
        splog, 'Reject flat ' + flatname[iflat] + ': No good traces'
      endif
      
      flatstruct[iflat].fflat = ptr_new(fflat)
      flatstruct[iflat].superflatset = ptr_new(superflatset)
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)
      
      flatstruct[iflat].nbright = nbright
      flatstruct[iflat].ymodel = ptr_new(ymodel)
      flatstruct[iflat].scatter = ptr_new(scatter)

      ;------------------------------------------------------------------
      ; Write information on flat field processing
      
      if (keyword_set(flatinfoname)) then begin
            
            write_spflat, flatinfoname, iflat, flatstruct, flathdr, $
                  arcname, nbright, ymodel, scatter, outdir=outdir,$
                  nowrite=nowrite, writeflatmodel=writeflatmodel
      
        ymodel = 0
      endif
      
      if keyword_set(traceflat) then begin
          flatstruct[iflat].xsol = traceflat_xsol
      endif

      obj_destroy,configuration
    endif
  endfor
  
  if (not keyword_set(plates)) and (not keyword_set(legacy)) and keyword_set(fibermask_bkup)  then fibermask=fibermask_bkup
  splog, 'Elapsed time = ', systime(1)-stime1, ' seconds', format='(a,f6.0,a)'
  return
end
;------------------------------------------------------------------------------
