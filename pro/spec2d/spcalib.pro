;+
; NAME:
;   spcalib
;
; PURPOSE:
;   Extract calibration frames.
;
; CALLING SEQUENCE:
;   spcalib, flatname, arcname, [pixflatname=, fibermask=, $
;    lampfile=, indir=, timesep=, ecalibfile=, plottitle=, $
;    arcinfoname=, flatinfoname=, arcstruct=, flatstruct=]
;
; INPUTS:
;   flatname   - Name(s) of flat-field SDSS image(s)
;   arcname    - Name(s) of arc SDSS image(s)
;
; OPTIONAL KEYWORDS:
;   pixflatname- Name of pixel-to-pixel flat, produced with SPFLATTEN.
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER].
;                Note this is not modified, but modified copies appear
;                in the returned structures ARCSTRUCT and FLATSTRUCT.
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in the IDL path.
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to '.'
;   timesep    - Maximum time separation between flats and arcs to pair them;
;                set to zero to disable this test; default to 7200 sec.
;   ecalibfile - opECalib file to pass to SDSSPROC
;   plottitle  - Prefix for titles in QA plots.
;   arcinfoname- File name (with path) to output arc extraction and fitting
;                information
;   flatinfoname-File name (with path) to output flat field extraction and
;                fitting information
;
; OUTPUTS:
;   arcstruct  - Structure array with extracted arc calibration information
;   flatstruct - Structure array with extracted flat calibration information
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Always pair arcs to the nearest good flat, and flats to the nearest good
;   arc (nearest in time, as defined by the TAI keyword in the FITS headers).
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   extract_image
;   fiberflat()
;   fitarcimage
;   fitarcimage_old
;   sdssproc
;   shift_trace()
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
    'LAMBDA', ptr_new(), $
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
    'FIBERMASK', ptr_new(), $
    'XSOL', ptr_new(), $
    'WIDTHSET', ptr_new(), $
    'FFLAT', ptr_new() )

   flatstruct = replicate(ftemp, nflat)

   return, flatstruct
end
;------------------------------------------------------------------------------

pro spcalib, flatname, arcname, pixflatname=pixflatname, fibermask=fibermask, $
 lampfile=lampfile, indir=indir, timesep=timesep, $
 ecalibfile=ecalibfile, plottitle=plottitle, $
 arcinfoname=arcinfoname, flatinfoname=flatinfoname, arcstruct, flatstruct

   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(timesep)) then timesep = 7200

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
      sdssproc, flatname[iflat], flatimg, flativar, indir=indir, $
       hdr=flathdr, pixflatname=pixflatname, nsatrow=nsatrow, fbadpix=fbadpix,$
       ecalibfile=ecalibfile

      ;-----
      ; Decide if this flat is bad:
      ;   Reject if more than 1% of the pixels are marked as bad.
      ;   Reject if more than 20 rows are saturated.

      qbadflat = 0
      if (fbadpix GT 0.01) then begin
         qbadflat = 1
         splog, 'Reject flat ' + flatname[iflat] + $
          ' (' + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels)'
      endif
      if (nsatrow GT 20) then begin
         qbadflat = 1
         splog, 'Reject flat ' + flatname[iflat] + $
          ' (' + string(format='(i4)', nsatrow) + ' saturated rows)'
      endif

      if (NOT keyword_set(fibermask)) then tmp_fibmask = bytarr(nfiber) $
       else tmp_fibmask = fibermask

      if (NOT qbadflat) then begin
         ;------------------------------------------------------------------
         ; Create spatial tracing from flat-field image
         ;------------------------------------------------------------------

         splog, 'Tracing 320 fibers in ',  flatname[iflat]
         xsol = trace320crude(flatimg, flativar, yset=ycen, maxdev=0.15, $
          fibermask=tmp_fibmask)

         splog, 'Fitting traces in ',  flatname[iflat]
         xy2traceset, ycen, xsol, tset, ncoeff=5, maxdev=0.1, $
             totalreject=totalreject

         if (totalreject GT 10000) then begin
           splog, 'WARNING: over 10000 pixels rejected!'
           qbadflat = 1
         endif
         traceset2xy, tset, ycen, xsol
      endif else begin
         xsol = 0
      endelse

      flatstruct[iflat].name = flatname[iflat]
      flatstruct[iflat].tai = sxpar(flathdr, 'TAI')
      flatstruct[iflat].qbad = qbadflat
      flatstruct[iflat].xsol = ptr_new(xsol)
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)

      ;---------------------------------------------------------
      ; Check to see if traces are separated by > 3 pixels
      ;

      if (qbadflat EQ 0) then begin 
        ntrace = (size(xsol, /dimens))[1]
        sep = xsol[*,1:ntrace-1] - xsol[*,0:ntrace-2]
        tooclose = where(sep LT 3)
        if (tooclose[0] NE -1) then begin
          splog, 'WARNING: Traces are not separated more than 3 pixels'
          flatstruct[iflat].qbad = 1
        endif
      endif

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
      sdssproc, arcname[iarc], arcimg, arcivar, indir=indir, $
       hdr=archdr, pixflatname=pixflatname, nsatrow=nsatrow, fbadpix=fbadpix, $
       ecalibfile=ecalibfile

splog,'Arc fbadpix ', fbadpix ; ???

      ;-----
      ; Decide if this arc is bad:
      ;   Reject if more than 1% of the pixels are marked as bad.
      ;   Reject if more than 100 rows are saturated.

      qbadarc = 0
      if (fbadpix GT 0.01) then begin
         qbadarc = 1
         splog, 'Reject arc ' + arcname[iarc] + $
          ' (' + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels)'
      endif
      if (nsatrow GT 100) then begin
         qbadarc = 1
         splog, 'Reject arc ' + arcname[iarc] + $
          ' (' + string(format='(i4)', nsatrow) + ' saturated rows)'
      endif

      tai = sxpar(archdr, 'TAI')

      ;-----
      ; Identify the nearest flat-field for this arc, which must be
      ; within TIMESEP seconds and be a good flat.

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
        tmp_fibmask = *(flatstruct[iflat].fibermask)



;
;	Calculate possible shift between arc and flat
;

         skiptrace = 20L
         nrow = (size(xsol))[1]
         nfiber = (size(xsol))[2]
         ysample = lindgen(nrow) # replicate(1,nfiber - 2*skiptrace)
         xsample = xsol[*,skiptrace:nfiber - skiptrace - 1]

         bestlag = shift_trace(arcimg, xsample, ysample, $
              lagrange=1.0, lagstep=0.1)

         if (abs(bestlag) GT 2.0) then begin
            qbadarc = 1
            splog, 'Reject arc: pixel shift is larger than 2 pixels'
         endif 
       endif

       if (NOT qbadarc) then begin
           splog, 'Shifting traces to fit arc by pixel shift of ', bestlag

           ;------------------------------------------------------------------
           ; Extract the arc image
           ;------------------------------------------------------------------

           splog, 'Extracting arc image with simple gaussian'
           sigma = 1.0
           proftype = 1 ; Gaussian
           highrej = 15
           lowrej = 15
           npoly = 1 ; maybe more structure
           wfixed = [1,1] ; Just fit the first gaussian term

           extract_image, arcimg, arcivar, xsol + bestlag, $
             sigma, flux, fluxivar, proftype=proftype, wfixed=wfixed, $
             highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1

           ;-------------------------------------------------------------------
           ; Compute correlation coefficient for this arc image
           ;-------------------------------------------------------------------

           splog, 'Searching for wavelength solution'
           aset = 0
; FOR NOW, REVERT TO THE OLD CODE! ???
          fitarcimage, flux, fluxivar, aset=aset, color=color, $
            lampfile=lampfile, fibermask=tmp_fibmask, bestcorr=bestcorr
;         fitarcimage_old, flux, fluxivar, aset=aset, color=color, $
;          lampfile=lampfile, fibermask=tmp_fibmask, bestcorr=bestcorr

          arcstruct[iarc].bestcorr = bestcorr

          if ((color EQ 'blue' AND bestcorr LT 0.5) $
            OR (color EQ 'red'  AND bestcorr LT 0.5) ) then begin
              qbadarc = 1
              splog, 'Reject arc ' + arcname[iarc] + $
               ' with correlation = ' + string(format='(i4)', bestcorr)
          endif

      endif

      if (NOT qbadarc) then begin

         ;---------------------------------------------------------------------
         ; Compute wavelength calibration
         ;---------------------------------------------------------------------

         arccoeff = 5

         splog, 'Searching for wavelength solution'
;       FOR NOW, REVERT TO THE OLD CODE! ???
         fitarcimage, flux, fluxivar, xpeak, ypeak, wset, $
          ncoeff=arccoeff, aset=aset, $
          color=color, lampfile=lampfile, fibermask=tmp_fibmask, $
          lambda=lambda, xdif_tset=xdif_tset
;         fitarcimage_old, flux, fluxivar, xpeak, ypeak, wset, $
;          ncoeff=arccoeff, aset=aset, $
;          color=color, lampfile=lampfile, fibermask=tmp_fibmask, $
;          lambda=lambda, xdif_tset=xdif_tset

         if (NOT keyword_set(wset)) then begin
            splog, 'Wavelength solution failed'
            qbadarc = 1
         endif else begin

            nfitcoeff = color EQ 'red' ? 4 : 3
            dispset = fitdispersion(flux, fluxivar, xpeak, $
               sigma=1.0, ncoeff=nfitcoeff, xmin=0.0, xmax=2047.0)

            arcstruct[iarc].dispset = ptr_new(dispset)
            arcstruct[iarc].wset = ptr_new(wset)
            arcstruct[iarc].nmatch = N_elements(lambda)
            arcstruct[iarc].lambda = ptr_new(lambda)
            arcstruct[iarc].tsep = tsep
            arcstruct[iarc].xpeak = ptr_new(xpeak)
            arcstruct[iarc].xdif_tset = ptr_new(xdif_tset)
            arcstruct[iarc].fibermask = ptr_new(tmp_fibmask) 

            ;------------------------------------------------------------------
            ; Write information on arc lamp processing
            ;

            if (keyword_set(arcinfoname)) then begin

              sxaddpar, archdr, 'FBADPIX', fbadpix, $
                  'Fraction of bad pixels in raw image'
              sxaddpar, archdr, 'BESTCORR', bestcorr, $
                  'Best Correlation coefficient'

              arcinfofile = string(format='(a,i8.8,a)',arcinfoname, $
                 sxpar(archdr, 'EXPOSURE'), '.fits')

              mwrfits, flux, arcinfofile, archdr, /create
              mwrfits, [transpose(lambda), xpeak], arcinfofile
              mwrfits, wset, arcinfofile
              mwrfits, fibermask, arcinfofile 
            endif

         endelse

;         qaplot_arcline, xdif_tset, lambda, filename=arcname[iarc], color=color

      endif

      arcstruct[iarc].name = arcname[iarc]
      arcstruct[iarc].tai = tai
      arcstruct[iarc].iflat = iflat
      arcstruct[iarc].qbad = qbadarc

   endfor

   ;---------------------------------------------------------------------------
   ; LOOP THROUGH FLATS + CREATE FIBERFLATS
   ;---------------------------------------------------------------------------

   for iflat=0, nflat-1 do begin

      splog, iflat+1, nflat, $
       format='("Create fiberflats for flat #",I3," of",I3)'

      ;-----
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

         wset = *(arcstruct[iarc].wset)
         xsol = *(flatstruct[iflat].xsol)
         tmp_fibmask = *(flatstruct[iflat].fibermask)

         ;---------------------------------------------------------------------
         ; Read flat-field image (again)
         ;---------------------------------------------------------------------

         ; If there is only 1 flat image, then it's still in memory
         if (nflat GT 1) then begin
            splog, 'Reading flat ', flatname[iflat]
            sdssproc, flatname[iflat], flatimg, flativar, indir=indir, $
             hdr=flathdr, pixflatname=pixflatname, ecalibfile=ecalibfile
         endif

         ;---------------------------------------------------------------------
         ; Extract the flat-field image
         ;---------------------------------------------------------------------

         splog, 'Extracting flat-field image with simple gaussian'
         sigma = 1.0
         proftype = 1 ; Gaussian
         highrej = 15
         lowrej = 15
         npoly = 1  ; just fit flat background to each row
         wfixed = [1,1] ; Just fit the first gaussian term

         extract_image, flatimg, flativar, xsol, sigma, $
          flux, fluxivar, proftype=proftype, wfixed=wfixed, $
          highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1, $
          ansimage=ansimage



         junk = where(flux GT 1.0e5, nbright)
         splog, 'Found ', nbright, ' bright pixels in extracted flat ', $
          flatname[iflat], format='(a,i7,a,a)'

         ;---------------------------------------------------------------------
         ; Compute fiber-to-fiber flat-field variations
         ;---------------------------------------------------------------------

         ntrace = (size(flux, /dimens))[1]

;         fflat = fiberflat(flux, fluxivar, wset, fibermask=tmp_fibmask, $
;          plottitle=plottitle)
         fflat = fiberflat(flux, fluxivar, wset, fibermask=tmp_fibmask, $
          /dospline, plottitle=plottitle)

         if (n_elements(fflat) EQ 1) then begin
            flatstruct[iflat].qbad  = 1
            splog, 'Reject flat ' + flatname[iflat] + $
               ':  No good traces?!?'
         endif

         widthset = fitflatwidth(flux, fluxivar, ansimage, tmp_fibmask, $
               ncoeff=5, sigma=sigma)

         flatstruct[iflat].fflat = ptr_new(fflat)
         flatstruct[iflat].widthset = ptr_new(widthset)
         flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)

         ;------------------------------------------------------------------
         ; Write information on flat field processing

         if (keyword_set(flatinfoname)) then begin

           sxaddpar, flathdr, 'NBRIGHT', nbright, $
               'Number of bright pixels (>10^5) in extracted flat-field'

           flatinfofile = string(format='(a,i8.8,a)',flatinfoname, $
                 sxpar(flathdr, 'EXPOSURE'), '.fits')

           mwrfits, fflat, flatinfofile, flathdr, /create
           mwrfits, tset, flatinfofile
           mwrfits, fibermask, flatinfofile
         endif

      endif

   endfor

   splog, 'Elapsed time = ', systime(1)-stime1, ' seconds', format='(a,f6.0,a)'

   return
end
;------------------------------------------------------------------------------
