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
    'FIBERMASK', ptr_new() )

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
    'SUPERFLATSET', ptr_new() )
    
  flatstruct = replicate(ftemp, nflat)
  
  return, flatstruct
end
;------------------------------------------------------------------------------

pro spcalib, flatname, arcname, fibermask=fibermask, cartid=cartid, $
    lampfile=lampfile, indir=indir, timesep=timesep, $
    ecalibfile=ecalibfile, plottitle=plottitle, $
    arcinfoname=arcinfoname, flatinfoname=flatinfoname, $
    arcstruct=arcstruct, flatstruct=flatstruct, $
    minflat=minflat, maxflat=maxflat, $
    writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, bbspec=bbspec
    
  if (NOT keyword_set(indir)) then indir = '.'
  if (NOT keyword_set(timesep)) then timesep = 7200
  if (NOT keyword_set(minflat)) then minflat = 0.8
  if (NOT keyword_set(maxflat)) then maxflat = 1.2
  
  stime1 = systime(1)
  
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
    sdssproc, flatname[iflat], flatimg, flativar, indir=indir, hdr=flathdr, $
      /applybias, /applypixflat, nsatrow=nsatrow, fbadpix=fbadpix,$
      ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat,/applycrosstalk
      
    configuration=obj_new('configuration', sxpar(flathdr, 'MJD'))
    
    ;-----
    ; Decide if this flat is bad
    
    qbadflat = reject_flat(flatimg, flathdr, nsatrow=nsatrow, fbadpix=fbadpix, $
      percent80thresh=configuration->spcalib_reject_calib_percent80thresh())

    if (NOT keyword_set(fibermask)) then tmp_fibmask = 0 $
     else tmp_fibmask = fibermask

    if (NOT qbadflat) then begin
      ;------------------------------------------------------------------
      ; Create spatial tracing from flat-field image
      ;------------------------------------------------------------------
    
      splog, 'Tracing fibers in ', flatname[iflat]
      xsol = trace320crude(flatimg, flativar, yset=ycen, maxdev=1.0, $ ;0.15, $
       fibermask=tmp_fibmask, cartid=cartid, xerr=xerr, $
       flathdr=flathdr, $
       padding=configuration->spcalib_trace320crude_padding(), $
       plottitle=plottitle+' Traces '+flatname[iflat])
        
      splog, 'Fitting traces in ', flatname[iflat]
      ntrace = (size(xsol, /dimens))[1]
      outmask = 0
      ; Ignore values whose central point falls on a bad pixel
      ;inmask = flativar[xsol,ycen] GT 0
      ; ASB: New recipe for inmask, just masking fully useless rows,
      ;      since trace320crude has already done clever fill-ins:
      inmask = (total(flativar gt 0., 1) gt 0.) # replicate(1B, ntrace)
      xy2traceset, ycen, xsol, tset, $
       ncoeff=configuration->spcalib_xy2traceset_ncoeff(), $
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
    endif else begin
      xsol = 0
      flatstruct[iflat].qbad = 1
    endelse
    
    ;----------
    ; Verify that traces are separated by > 3 pixels
    
    if (qbadflat EQ 0) then begin
      sep = xsol[*,1:ntrace-1] - xsol[*,0:ntrace-2]
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
      extract_image, flatimg, flativar, xsol, sigma, flux, fluxivar, $
       proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
       npoly=npoly, relative=1, ansimage=ansimage, reject=[0.1, 0.6, 0.6], $
       chisq=chisq3

      widthset3 = fitflatwidth(flux, fluxivar, ansimage, tmp_fibmask, $
       ncoeff=configuration->spcalib_fitflatwidth_ncoeff(), sigma=sigma, $
       medwidth=medwidth, $
       mask=configuration->spcalib_fitflatwidth_mask(flux,fluxivar), $
       inmask=configuration->spcalib_fitflatwidth_inmask(flux,fluxivar), $
       /double)
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
  
    splog, iarc+1, narc, format='("Extracting arc #",I3," of",I3)'
    
    ;---------------------------------------------------------------------
    ; Read the arc
    ;---------------------------------------------------------------------
    
    splog, 'Reading arc ', arcname[iarc]
    
    sdssproc, arcname[iarc], arcimg, arcivar, indir=indir, hdr=archdr, $
      /applybias, /applypixflat, nsatrow=nsatrow, fbadpix=fbadpix, $
      ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat,/applycrosstalk
    ny = (size(arcimg,/dimens))[1]
      
    configuration=obj_new('configuration', sxpar(archdr, 'MJD'))
    
    splog, 'Fraction of bad pixels in arc = ', fbadpix
    
    ;----------
    ; Decide if this arc is bad
    
    qbadarc = reject_arc(arcimg, archdr, nsatrow=nsatrow, fbadpix=fbadpix)
    
    ;----------
    ; Identify the nearest flat-field for this arc, which must be
    ; within TIMESEP seconds and be a good flat.
    
    get_tai, archdr, tai_beg, tai_mid, tai_end
    tai = tai_mid
    
    iflat = -1
    igood = where(flatstruct.qbad EQ 0)
    if (igood[0] NE -1) then begin
      tsep = min( abs(tai - flatstruct[igood].tai), ii )
      if (tsep LE timesep AND timesep NE 0) then iflat = igood[ii]
    endif
    
    if (iflat GE 0) then begin
      splog, 'Arc ' + arcname[iarc] + ' paired with flat ' + flatname[iflat]
    endif else begin
      splog, 'Arc ' + arcname[iarc] + ' paired with no flat'
      qbadarc = 1
    endelse
    
    if (NOT qbadarc) then begin
      xsol = *(flatstruct[iflat].xsol)
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
        splog, 'Reject arc ' + arcname[iarc] + $
          ': Pixel shift = ', bestlag
      endif
    endif
    
    if (NOT qbadarc) then begin
      splog, 'Shifting traces with match_trace', bestlag
      ;         splog, 'Shifting traces to fit arc by pixel shift of ', bestlag
      
      ;---------------------------------------------------------------------
      ; Extract the arc image
      ;---------------------------------------------------------------------
      
      traceset2xy, widthset, xx, sigma2
      
      highrej = 15
      lowrej = 15

      wfixed = [1,0] ; ASB: Don't fit for width terms.

      splog, 'Extracting arc'
      extract_bundle_image, arcimg, arcivar, xcor, sigma2, $
        flux, fluxivar, proftype=proftype, wfixed=wfixed, $
        highrej=highrej, lowrej=lowrej, npoly=2L, relative=1, $
        reject=[0.1, 0.6, 0.6], ymodel=ymodel, nperbun=20L, buffsize=8L
        
      ; flag to determine whether or not to do 2-phase arc solution:
      twophase = sxpar(archdr, 'TWOPHASE')
      if keyword_set(twophase) then splog, 'Setting 2-phase readout flag'

      ;---------------------------------------------------------------------
      ; Compute correlation coefficient for this arc image
      ;---------------------------------------------------------------------

      splog, 'Searching for wavelength solution'
      aset = 0
      fitarcimage, flux, fluxivar, aset=aset, color=color, $
       lampfile=lampfile, fibermask=tmp_fibmask, bestcorr=bestcorr, $
       acoeff=configuration->spcalib_arcfitguess_acoeff(color), $
       dcoeff=configuration->spcalib_arcfitguess_dcoeff(color), $
       wrange=configuration->spcalib_fitarcimage_wrange(color), $
       twophase=twophase

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
       lambda=lambda, rejline=rejline, xdif_tset=xdif_tset, $
       acoeff=configuration->spcalib_arcfitguess_acoeff(color), $
       dcoeff=configuration->spcalib_arcfitguess_dcoeff(color), $
       wrange=configuration->spcalib_fitarcimage_wrange(color), $
       twophase=twophase

      if (NOT keyword_set(wset)) then begin
        splog, 'Wavelength solution failed'
        qbadarc = 1
      endif else begin

        nfitcoeff = configuration->spcalib_ncoeff(color)
        ilamp = where(rejline EQ '')
        dispset = fitdispersion(flux, fluxivar, xpeak[*,ilamp], $
          sigma=configuration->spcalib_sigmaguess(), ncoeff=nfitcoeff, $
          xmin=0.0, xmax=ny-1, $
          medwidth=wsigarr, numbundles=ntrace/20) ; Hard-wires 20 fibers/bundle???

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
        
        ;------------------------------------------------------------------
        ; Write information on arc lamp processing
        
        if (keyword_set(arcinfoname)) then begin
          sxaddpar, archdr, 'FBADPIX', fbadpix, $
            'Fraction of bad pixels in raw image'
          sxaddpar, archdr, 'BESTCORR', bestcorr, $
            'Best Correlation coefficient'
            
          arcinfofile = string(format='(a,i8.8,a)',arcinfoname, $
            sxpar(archdr, 'EXPOSURE'), '.fits')
            
          mwrfits, flux, arcinfofile, archdr, /create
          mwrfits, [transpose(lambda), xpeak], arcinfofile
          mwrfits, *arcstruct[iarc].wset, arcinfofile
          mwrfits, *arcstruct[iarc].fibermask, arcinfofile
          mwrfits, *arcstruct[iarc].dispset, arcinfofile
          spawn, ['gzip', '-f', arcinfofile], /noshell
         ; ASB: write arc image model info if requested:
          if keyword_set(writearcmodel) then begin
             arcmodelfile = string(format='(a,i8.8,a)',arcinfoname + $
               'MODELIMG-', sxpar(archdr, 'EXPOSURE'), '.fits')
             mwrfits, arcimg, arcmodelfile, /create
             mwrfits, arcivar, arcmodelfile
             mwrfits, ymodel, arcmodelfile
             spawn, ['gzip', '-f', arcmodelfile], /noshell
          endif
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
      if (tsep LE timesep AND timesep NE 0) then iarc = igood[ii]
      flatstruct[iflat].tsep = tsep
    endif
    
    if (iarc GE 0) then begin
      splog, 'Flat ' + flatname[iflat] + ' paired with arc ' + arcname[iarc]
    endif else begin
      splog, 'Flat ' + flatname[iflat] + ' paired with no arc'
      flatstruct[iflat].qbad = 1 ; Flat is bad if no companion arc exists
    endelse
    
    flatstruct[iflat].iarc = iarc
    
    if (NOT flatstruct[iflat].qbad) then begin
    
      widthset = *(flatstruct[iflat].widthset)
      wset = *(arcstruct[iarc].wset)
      xsol = *(flatstruct[iflat].xsol)
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
          /applybias, /applypixflat, ecalibfile=ecalibfile, $
          minflat=minflat, maxflat=maxflat,/applycrosstalk
      endif
      configuration=obj_new('configuration',sxpar(flathdr, 'MJD'))

      ;---------------------------------------------------------------------
      ; Extract the flat-field image
      ;---------------------------------------------------------------------
      
      traceset2xy, widthset, xx, sigma2   ; sigma2 is real width
      highrej = 15
      lowrej = 15
      npoly = 5 ; Fit 5 terms to background, just get best model
;      wfixed = [1,1] ; Fit gaussian plus both derivatives
      wfixed = [1,0] ; Do not refit for Gaussian widths, only flux ???

      extract_bundle_image, flatimg, flativar, xsol, sigma2, flux, fluxivar, $
        proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
        npoly=2L, relative=1, chisq=schisq, ansimage=ansimage2, $
        reject=[0.1, 0.6, 0.6], ymodel=ymodel, nperbun=20L, buffsize=8L

      if (keyword_set(bbspec)) then begin
         basisfile = 'spBasisPSF-*-'+strmid(arcstruct[iarc].name,4,11)+'.fits'
         tmproot = 'tmp-'+strmid(flatstruct[iflat].name,4,11)
         bbspec_extract, flatimg, flativar, bbflux, bbfluxivar, $
          basisfile=basisfile, ximg=xsol, ymodel=bb_ymodel, $
          tmproot=tmproot, /batch ; ??? set batch

         ; Deal with case of only the first few spectra being re-extracted...
         dims = size(bbflux,/dimens)
         flux[0:dims[0]-1,0:dims[1]-1] = bbflux
         fluxivar[0:dims[0]-1,0:dims[1]-1] = bbfluxivar $
          * (fluxivar[0:dims[0]-1,0:dims[1]-1] GT 0) ; <- Retain old rejection

         outfile = 'ymodel-'+strmid(flatstruct[iflat].name,4,11)+'.fits'
         mwrfits, bb_ymodel, outfile, /create
         mwrfits, ymodel, outfile
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
      
      ;------------------------------------------------------------------
      ; Write information on flat field processing
      
      if (keyword_set(flatinfoname)) then begin
      
        sxaddpar, flathdr, 'NBRIGHT', nbright, $
          'Number of bright pixels (>10^5) in extracted flat-field'
          
        flatinfofile = string(format='(a,i8.8,a)',flatinfoname, $
          sxpar(flathdr, 'EXPOSURE'), '.fits')
          
        mwrfits, *flatstruct[iflat].fflat, flatinfofile, flathdr, /create
        mwrfits, *flatstruct[iflat].tset, flatinfofile
        mwrfits, *flatstruct[iflat].fibermask, flatinfofile
        mwrfits, *flatstruct[iflat].widthset, flatinfofile
        mwrfits, *flatstruct[iflat].superflatset, flatinfofile
        spawn, ['gzip', '-f', flatinfofile], /noshell
        ; ASB: write flat image model info if requested:
        if keyword_set(writeflatmodel) then begin
           flatmodelfile = string(format='(a,i8.8,a)',flatinfoname + $
             'MODELIMG-', sxpar(flathdr, 'EXPOSURE'), '.fits')
           mwrfits, flatimg, flatmodelfile, /create
           mwrfits, flativar, flatmodelfile
           mwrfits, ymodel + scatter, flatmodelfile
           spawn, ['gzip', '-f', flatmodelfile], /noshell
        endif
        ymodel = 0
      endif
      
      obj_destroy,configuration
    endif
  endfor
  
  splog, 'Elapsed time = ', systime(1)-stime1, ' seconds', format='(a,f6.0,a)'
  return
end
;------------------------------------------------------------------------------
