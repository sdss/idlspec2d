;+
; NAME:
;   spcalib
;
; PURPOSE:
;   Extract calibration frames.
;
; CALLING SEQUENCE:
;   spcalib, flatname, arcname, pixflatname=, fibermask=, $
;    lampfile=, indir=, timesep=, arcstruct, flatstruct
;
; INPUTS:
;   flatname   - Name(s) of flat-field SDSS image(s)
;   arcname    - Name(s) of arc SDSS image(s)
;
; OPTIONAL KEYWORDS:
;   pixflatname- Name of pixel-to-pixel flat, produced with SPFLATTEN.
;   fibermask  - Mask of 0 for bad fibers and 1 for good fibers [NFIBER]
;   lampfile   - Name of file describing arc lamp lines;
;                default to the file 'lamphgcdne.dat' in the IDL path.
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to '.'
;   timesep    - Maximum time separation between flats and arcs to pair them;
;                default to 3600 sec.
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
    'QBAD', 0B, $
    'IFLAT', -1, $
    'BESTCORR', 0.0, $
    'NMATCH', 0L, $
    'LAMBDA', ptr_new(), $
    'XPEAK', ptr_new(), $
    'XDIF_TSET', ptr_new(), $
    'WSET', ptr_new() )

   arcstruct = replicate(ftemp, narc)

   return, arcstruct
end
;------------------------------------------------------------------------------
function create_flatstruct, nflat

   ftemp = create_struct( name='FLAT_STRUCT', $
    'NAME', '', $
    'TAI', 0D, $
    'QBAD', 0, $
    'IARC', -1, $
    'FIBERMASK', ptr_new(), $
    'XSOL', ptr_new(), $
    'FFLAT', ptr_new() )

   flatstruct = replicate(ftemp, nflat)

   return, flatstruct
end
;------------------------------------------------------------------------------

pro spcalib, flatname, arcname, pixflatname=pixflatname, fibermask=fibermask, $
 lampfile=lampfile, indir=indir, timesep=timesep, arcstruct, flatstruct

   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(timesep)) then timesep = 3600

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
       hdr=flathdr, pixflatname=pixflatname, nsatrow=nsatrow, fbadpix=fbadpix

      ;-----
      ; Decide if this flat is bad:
      ;   Reject if more than 1% of the pixels are marked as bad.
      ;   Reject if more than 10 rows are saturated.

      qbadflat = 0
      if (fbadpix GT 0.01) then begin
         qbadflat = 1
         splog, 'Reject flat ' + flatname[iflat] + $
          ' (' + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels)'
      endif
      if (nsatrow GT 10) then begin
         qbadflat = 1
         splog, 'Reject flat ' + flatname[iflat] + $
          ' (' + string(format='(i4)', nsatrow) + ' saturated rows)'
      endif

      if (NOT qbadflat) then begin
         ;------------------------------------------------------------------
         ; Create spatial tracing from flat-field image
         ;------------------------------------------------------------------

         splog, 'Tracing 320 fibers in ',  flatname[iflat]
         xsol = trace320crude(flatimg, flativar, yset=ycen, maxdev=0.15)

         splog, 'Fitting traces in ',  flatname[iflat]
         xy2traceset, ycen, xsol, tset, ncoeff=5, maxdev=0.1
         traceset2xy, tset, ycen, xsol
      endif else begin
         xsol = 0
      endelse

      flatstruct[iflat].name = flatname[iflat]
      flatstruct[iflat].tai = sxpar(flathdr, 'TAI')
      flatstruct[iflat].qbad = qbadflat
      flatstruct[iflat].xsol = ptr_new(xsol)

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
       hdr=archdr, pixflatname=pixflatname, nsatrow=nsatrow, fbadpix=fbadpix

splog,'Arc fbadpix ', fbadpix ; ???

      ;-----
      ; Decide if this arc is bad:
      ;   Reject if more than 1% of the pixels are marked as bad.
      ;   Reject if more than 40 rows are saturated.

      qbadarc = 0
      if (fbadpix GT 0.01) then begin
         qbadarc = 1
         splog, 'Reject arc ' + flatname[iflat] + $
          ' (' + string(format='(i3)', fix(fbadpix*100)) + '% bad pixels)'
      endif
      if (nsatrow GT 40) then begin
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
         if (tsep LE timesep) then iflat = igood[ii]
      endif

      if (iflat GE 0) then begin
         splog, 'Arc ' + arcname[iarc] + ' paired with flat ' + flatname[iflat]
      endif else begin
         splog, 'Arc ' + arcname[iarc] + ' paired with no flat'
         qbadarc = 1
      endelse

      if (NOT qbadarc) then begin

         xsol = *(flatstruct[iflat].xsol)

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

         extract_image, arcimg, arcivar, xsol, sigma, flux, fluxivar, $
          proftype=proftype, wfixed=wfixed, $
          highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1

         ;-------------------------------------------------------------------
         ; Compute correlation coefficient for this arc image
         ;-------------------------------------------------------------------

         splog, 'Searching for wavelength solution'
         aset = 0
; FOR NOW, REVERT TO THE OLD CODE! ???
         fitarcimage, flux, fluxivar, aset=aset, $
          color=color, lampfile=lampfile, bestcorr=bestcorr
;         fitarcimage_old, flux, fluxivar, aset=aset, $
;          color=color, lampfile=lampfile, bestcorr=bestcorr

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
          color=color, lampfile=lampfile, lambda=lambda, xdif_tset=xdif_tset
;         fitarcimage_old, flux, fluxivar, xpeak, ypeak, wset, $
;          ncoeff=arccoeff, aset=aset, $
;          color=color, lampfile=lampfile, lambda=lambda, xdif_tset=xdif_tset

         if (NOT keyword_set(wset)) then begin
            splog, 'Wavelength solution failed'
            qbadarc = 1
         endif else begin
            arcstruct[iarc].wset = ptr_new(wset)
            arcstruct[iarc].nmatch = N_elements(lambda)
            arcstruct[iarc].lambda = ptr_new(lambda)
            arcstruct[iarc].xpeak = ptr_new(xpeak)
            arcstruct[iarc].xdif_tset = ptr_new(xdif_tset)
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
         if (tsep LE timesep) then iarc = igood[ii]
      endif

      if (iarc GE 0) then $
       splog, 'Flat ' + flatname[iflat] + ' paired with arc ' + arcname[iarc] $
      else $
       splog, 'Flat ' + flatname[iflat] + ' paired with no arc'

      flatstruct[iflat].iarc = iarc

      if (NOT flatstruct[iflat].qbad AND iarc NE -1) then begin

         wset = *(arcstruct[iarc].wset)
         xsol = *(flatstruct[iflat].xsol)

         ;---------------------------------------------------------------------
         ; Read flat-field image (again)
         ;---------------------------------------------------------------------

         ; If there is only 1 flat image, then it's still in memory
         if (nflat GT 1) then begin
            splog, 'Reading flat ', flatname[iflat]
            sdssproc, flatname[iflat], flatimg, flativar, indir=indir, $
             hdr=flathdr, pixflatname=pixflatname
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
          flux, fluxivar, $
          proftype=proftype, wfixed=wfixed, $
          highrej=highrej, lowrej=lowrej, npoly=npoly, relative=1

         junk = where(flux GT 1.0e5, nbright)
         splog, 'Found ', nbright, ' bright pixels in extracted flat ', $
          flatname[iflat], format='(a,i7,a,a)'

         ;---------------------------------------------------------------------
         ; Compute fiber-to-fiber flat-field variations
         ;---------------------------------------------------------------------

         ntrace = (size(flux, /dimens))[1]
         if (keyword_set(fibermask)) then fmask = fibermask $
          else fmask = bytarr(ntrace)

;         fflat = fiberflat(flux, fluxivar, wset, fibermask=fmask)
         fflat = fiberflat(flux, fluxivar, wset, fibermask=fmask, $
          /dospline)

;         qaplot_fflat, fflat, wset, filename=flatname[iflat]

         flatstruct[iflat].fflat = ptr_new(fflat)
         flatstruct[iflat].fibermask = ptr_new(fmask)

      endif

   endfor

   return
end
;------------------------------------------------------------------------------
