;+
; NAME:
;   sdssproc
;
; PURPOSE:
;   Read in Raw SDSS files, and process with opConfig, opECalib, opBC par files.
;
; CALLING SEQUENCE:
;   sdssproc, infile, [image, invvar, indir=, $
;    outfile=, varfile=, nsatrow=, fbadpix=, $
;    hdr=hdr, configfile=, ecalibfile=, bcfile=, $
;    pixflatname=, minflat=, spectrographid=, color= ]
;
; INPUTS:
;   infile     - Raw SDSS file name
;
; OPTIONAL KEYWORDS:
;   indir      - Input directory for INFILE
;   outfile    - Calibrated 2d frame, after processing
;   varfile    - Inverse variance frame after processing
;   nsatrow    - Number of saturated rows, assuming that a row is saturated
;                if at least 20 of its pixels are above saturation level
;   fbadpix    - Fraction of bad pixels, not including bad columns
;   hdr        - Header returned in memory
;   configfile - Default to "opConfig.par"
;   ecalibfile - Default to "opECalib.par"
;   bcfile     - Default to "opBC.par"
;   minflat    - Minimum values allowed for pixflat
;                   (lower values of pixflat are set to 0 invvar)
;   pixflatname- Name of pixel-to-pixel flat, produced with SPFLATTEN.
;
; OUTPUTS:
;   image      - Processed 2d image
;   invvar     - Associated inverse variance
;   spectrographid - Return spectrograph ID (1 or 2)
;   color      - Return spectrograph color ('red' or 'blue')
;
; COMMENTS:
;   Only the header is read from the image if IMAGE, INVVAR, OUTFILE and
;   VARFILE are all not set.
;
;   Required header keywords: EXPTIME.
;
; BUGS:
;   The open-shutter correction SMEARIMG will include smeared data from
;   any cosmic rays, which is wrong.  At the minimum, I could interpolate
;   over A/D saturations (in ADMASK) before constructing SMEARIMG.
;
; PROCEDURES CALLED:
;   djs_iterstat
;   headfits()
;   idlspec2d_version()
;   idlutils_version()
;   rdss_fits()
;   readfits()
;   splog
;   sxaddpar
;   sxpar()
;   writefits
;   yanny_free
;   yanny_read
;
; REVISION HISTORY:
;   13-May-1999  Written by Scott Burles & David Schlegel, Apache Point.
;   08-Sep-1999  Modified to read Yanny param files instead of FITS
;                versions of the same (DJS).
;   01-Dec-1999  Added version stamping (DJS).
;   07-Dec-1999  Mask neighbors of pixels that saturated the A/D converter.
;                Identify blead trails and mask from that row up (DJS).
;   10-Dec-1999  Test if the shutter was open during readout, and try
;                to correct the light for that (DJS).
;   04-Feb-2000  Declare that the shutter was open if it is a >640 sec
;                exposure taken before MJD=51570 (DJS).
;-
;------------------------------------------------------------------------------

function findopfile, expres, mjd, indir

   cd, indir, current=olddir
   files = findfile(expres)
   nfiles = n_elements(files)
   if (nfiles EQ 0) then begin
     cd, olddir
     message, 'ABORT: Cannot find opFile '+expres
   endif

   mjdlist = lonarr(nfiles)
   for i=0,nfiles-1 do begin
     yanny_read, files[i], info, hdr=hdr
     mjdlist[i] = yanny_par(hdr, 'mjd')
     yanny_free, info
   endfor 

   diff = mjd - mjdlist 
   score = min(diff  + 10000000L*(diff LT 0), bestmjd)

   cd, olddir
   return, files[bestmjd]
      
end

;------------------------------------------------------------------------------
pro sdssproc, infile, image, invvar, indir=indir, $
 outfile=outfile, varfile=varfile, nsatrow=nsatrow, fbadpix=fbadpix, $
 hdr=hdr, configfile=configfile, ecalibfile=ecalibfile, bcfile=bcfile, $
 pixflatname=pixflatname, spectrographid=spectrographid, color=color, $
 minflat=minflat

   if (N_params() LT 1) then begin
      print, 'Syntax - sdssproc, infile, [image, invvar, indir=, $'
      print, ' outfile=, varfile=, nsatrow=, fbadpix=, $' 
      print, ' hdr=, configfile=, ecalibfile=, bcfile=, $'
      print, ' pixflatname=, spectrographid=, color= ]'
      return
   endif


   readimg = arg_present(image) OR keyword_set(outfile)
   readivar = arg_present(invvar) OR keyword_set(varfile) $
    OR arg_present(nsatrow) OR arg_present(fbadpix)


   if (keyword_set(indir)) then fullname = filepath(infile, root_dir=indir) $
    else fullname = infile
   fullname = (findfile(fullname, count=ct))[0]
   if (ct NE 1) then $
    message, 'Cannot find image ' + infile

   if (readimg OR readivar) then $
    rawdata = rdss_fits(fullname, hdr, /nofloat) $
   else $
    hdr = headfits(fullname)

   mjd = sxpar(hdr, 'MJD')

   pp = filepath('', root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')

   if (NOT keyword_set(configfile)) then $
       configfile = findopfile('opConfig*par',mjd,pp)
   if (NOT keyword_set(ecalibfile)) then $
       ecalibfile = findopfile('opECalib*par',mjd,pp)
   if (NOT keyword_set(bcfile)) then $
       bcfile = findopfile('opBC*par',mjd,pp)

   naxis = sxpar(hdr,'NAXIS*')
   if (naxis[0] NE 2128 OR naxis[1] NE 2069) then $
    splog, 'WARNING: Expecting 2128x2069, found '+string(naxis[0])+'x'$
     +string(naxis[1])

   ;------
   ; Determine which CCD from the file name itself, using either the
   ; numbering scheme (01,02,03,04) or naming scheme (b1,r2,b2,r1).
   ; Very bad form, but this information is not in the header since
   ; the CAMERAS keyword is sometimes wrong.

   i = rstrpos(infile, '-')
   if (i[0] EQ -1 OR i-2 LT 0) then $
    message, 'Cannot determine CCD number from file name ' + infile

   camnames = ['b1', 'r2', 'b2', 'r1']
   camnums = ['01', '02', '03', '04']

   indx = where(strmid(infile, i-2, 2) EQ camnames, ct)
   if (ct NE 1) then $
     indx = where(strmid(infile, i-2, 2) EQ camnums, ct)
   if (ct NE 1) then $
    message, 'Cannot determine CCD number from file name ' + infile

;   cameras = strtrim( sxpar(hdr, 'CAMERAS'), 2 )
   case camnames[indx[0]] of
    'b1': begin
          spectrographid = 1
          color = 'blue'
          end
    'r2': begin
;
;;	r2 means red on spectrograph 2
;;  and r1 means red on spectrograph 1
;;  I don't think any red reductions have been working since this
;;  was changed.
;;
;
          spectrographid = 2
          color = 'red'
          end
    'b2': begin
          spectrographid = 2
          color = 'blue'
          end
    'r1': begin
          spectrographid = 1
          color = 'red'
          end
   endcase
   camcol = indx[0] + 1
   camrow = 0

   sxaddpar, hdr, 'CAMROW', camrow
   sxaddpar, hdr, 'CAMCOL', camcol
   sxaddpar, hdr, 'TELESCOP', 'SDSS 2.5-M', ' Sloan Digital Sky Survey'
   sxaddpar, hdr, 'AUTHOR', 'Scott Burles & David Schlegel'
   sxaddpar, hdr, 'SPEC2D_V', idlspec2d_version(), ' Version of idlspec2d'
   sxaddpar, hdr, 'UTILS_V', idlutils_version(), ' Version of idlutils'
 
   ; Read in opConfig.par file
   ; Take the first entry for the configuration of each CCD in the event
   ; that there are several.

   yanny_read, filepath(configfile, root_dir=pp), pdata
   config = *pdata[0]
   yanny_free, pdata
   i = where(config.camrow EQ camrow AND config.camcol EQ camcol)
   config = config[i[0]]

   if (naxis[0] NE config.ncols OR naxis[1] NE config.nrows) then $
      splog, 'WARNING! Config file dimensions do not match raw image'
   

   qexist = [config.amp0, config.amp1, config.amp2, config.amp3]

   ; Define the "overscan" regions
   sover = [config.soverscan0, config.soverscan1, config.soverscan2, $
    config.soverscan3]
   nover = [config.noverscan0, config.noverscan1, config.noverscan2, $
    config.noverscan3]

   ; Define the "mapped overscan" regions
   smapover = [config.smapoverscan0, config.smapoverscan1, $
    config.smapoverscan2, config.smapoverscan3]
   nmapover = [config.nmapoverscan0, config.nmapoverscan1, $
    config.nmapoverscan2, config.nmapoverscan3]

   ; Define the "overscan rows" (at the bottom of the CCD)
   soverrow = [config.soverscanrows0, config.soverscanrows1, $
    config.soverscanrows2, config.soverscanrows3]
   noverrow = [config.noverscanrows0, config.noverscanrows1, $
    config.noverscanrows2, config.noverscanrows3]

   ; Data position in the original image
   sdatarow = [config.sdatarow0, config.sdatarow1, $
    config.sdatarow2, config.sdatarow3]
   sdatacol = [config.sdatasec0, config.sdatasec1, config.sdatasec2, $
    config.sdatasec3]
   nrow = [config.ndatarow0, config.ndatarow1, $
    config.ndatarow2, config.ndatarow3]
   ncol = [config.ndatasec0, config.ndatasec1, config.ndatasec2, $
    config.ndatasec3]

   ; Data position in the final (trimmed) image
   srow = [config.sccdrowsec0, config.sccdrowsec1, $
    config.sccdrowsec2, config.sccdrowsec3]
   scol = [config.sccdcolsec0, config.sccdcolsec1, $
    config.sccdcolsec2, config.sccdcolsec3]

   if (naxis[1] EQ 2049) then begin
     splog, 'WARNING: NROWS is 2049, adjusting config entries' 
     sdatarow = sdatarow - 20
     noverrow = 1
   endif
    

   ; Read in ECalib File
   yanny_read, filepath(ecalibfile, root_dir=pp), pdata
   ecalib = *pdata[0]
   yanny_free, pdata
   ecalib = ecalib[ where(ecalib.camrow EQ camrow AND ecalib.camcol EQ camcol) ]

   gain = [ecalib.gain0, ecalib.gain1, ecalib.gain2, ecalib.gain3]
   readnoiseDN = [ecalib.readnoiseDN0, ecalib.readnoiseDN1, $
    ecalib.readnoiseDN2, ecalib.readnoiseDN3]
   fullWellDN = [ecalib.fullWellDN0, ecalib.fullWellDN1, $
    ecalib.fullWellDN2, ecalib.fullWellDN3]

   ; Construct the final image
   igood = where(qexist)
   nr = max((srow+nrow)[igood])
   nc = max((scol+ncol)[igood])
   if ((size(image))[0] NE 2) then image = fltarr(nc, nr) $
   else if ((size(image))[1] NE nc OR (size(image))[2] NE nr OR $
            (size(image))[3] NE 4) then image = fltarr(nc, nr) 

   yanny_read, filepath(bcfile, root_dir=pp), pdata
   bc = *pdata[0]
   yanny_free, pdata
   
   bchere = where(bc.camrow EQ camrow AND bc.camcol EQ camcol, nbc)
   if (nbc GT 0) then bc = bc[ bchere ]

   ;------
   ; Test to see if the shutter was open during readout if the exposure
   ; was longer than 640 seconds.

   exptime = sxpar(hdr, 'EXPTIME')
   flavor = sxpar(hdr, 'FLAVOR')
   qshutter = 0

   ; Toggle the variable QSHUTTER if the observation was taken before
   ; MJD=51570 and this was not a bias or dark exposure.

   if (exptime GT 640 AND (readimg OR readivar) $
    AND mjd GT 0 AND mjd LT 51570 $
    AND flavor NE 'bias' AND flavor NE 'dark') then qshutter = 1

   ; Look at the signal in the overscan rows (at the bottom of the CCD).
   ; Toggle the variable QSHUTTER if this appears to be true in any
   ; of the amplifiers.

;   if (exptime GT 640 AND (readimg OR readivar)) then begin
;      nskip = 2  ; Ignore the first and last NSKIP rows of these overscan rows
;      for iamp=0, 3 do begin
;         if (qexist[iamp] EQ 1) then begin
;            biasreg = rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
;                soverrow[iamp]+nskip:soverrow[iamp]+noverrow[iamp]-1-nskip]
;            biasvec = djs_median(biasreg, 2)
;            ; Count the number of "hot" overscan columns, hotter than 3-sigma
;            ; above the median
;            junk = where(biasvec GT median(biasvec) + 4, nhot)
;            if (nhot GE 15) then qshutter = 1 ; Flag the shutter as being open
;            splog, 'Number of hot overscan columns for amp', iamp, ' = ', nhot
;         endif
;      endfor
;   endif

   ;------
   ; Construct IMAGE

   for iamp=0, 3 do begin
      if (qexist[iamp] EQ 1) then begin
         if (readimg OR readivar) then begin
            if (nover[iamp] NE 0) then begin
               ; Use the "overscan" region
               biasreg = rawdata[sover[iamp]:sover[iamp]+nover[iamp]-1, $
                   sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1]
            endif else if (nmapover[iamp] NE 0) then begin
               ; Use the "mapped overscan" region
               biasreg = rawdata[smapover[iamp]:smapover[iamp]+nmapover[iamp]-1, $
                sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] 
            endif

            biasval = median( biasreg )
            djs_iterstat, biasreg, sigma=readoutDN

            splog, 'Measured read-noise in DN for amp#', iamp, ' = ', readoutDN

            ; Copy the data for this amplifier into the final image
            ; Subtract the bias (in DN), and then multiply by the gain
            ; Now image is in electrons

            image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                   srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
             (rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
                   sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] - biasval) $
             * gain[iamp]

            ; Add to the header
            sxaddpar, hdr, 'GAIN'+string(iamp,format='(i1)'), $
             gain[iamp], ' Gain in electrons per ADU'
            sxaddpar, hdr, 'RDNOISE'+string(iamp,format='(i1)'), $
             gain[iamp]*readnoiseDN[iamp], ' Readout noise in electrons'
         endif
      endif
   endfor

   ;------
   ; If the shutter was open during readout, try to correct for that.
   ; Construct an image of the smeared light on the CCD during readout.
   ; Note that we work from IMAGE, which already has the bias removed
   ; and is multiplied by the gain.

   if (qshutter) then begin
      splog, 'WARNING: Correcting for open shutter during readout '

      t1 = exptime ; Read time for entire frame
      t2 = 0.026976 ; Read time for one row of data (from Connie Rockosi)

      smearimg = 0 * image
      smearimg[*,0] = image[*,0]
      ny = (size(image,/dimens))[1]
      for i=1, ny-1 do begin
         ; Burles counter of row number...
         ;print, format='($, ".",i4.4,a5)', i, string([8b,8b,8b,8b,8b])

         smearimg[*,i] = smearimg[*,i-1] + image[*,i]
      endfor
      smearimg = (t2/t1) * smearimg
      image = image - smearimg

      splog, 'Median value of open-shutter contamination = ', median(smearimg)
      splog, 'Max value of open-shutter contamination = ', max(smearimg)
   endif

   ;------
   ; Construct INVVAR

   if (readivar) then begin
      if ((size(invvar))[0] NE 2) then $
       invvar = fltarr(nc, nr) $
      else if ((size(invvar))[1] NE nc OR (size(invvar))[2] NE nr OR $
       (size(invvar))[3] NE 4) then $
       invvar = fltarr(nc, nr) 

      ;------
      ; SATMASK = Mask for saturated the detector, 0=bad
      ; ADMASK = Mask for saturating the A/D converter (at 65535), 1=bad
      ; BCMASK = Mask for bad columns, 1=bad

      satmask = bytarr(nc, nr)
      admask = bytarr(nc, nr)
      bcmask = bytarr(nc, nr)
 
      for iamp=0, 3 do begin
         if (qexist[iamp] EQ 1) then begin

            satmask[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                     srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
              rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
                     sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] LT $
                     fullWellDN[iamp]

            admask[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                     srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
              rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
                     sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] EQ 65535

            if (qshutter) then begin
               ; Add to the variance image from the open shutter
               ; by adding 1% of that signal^2 to the variance.
               ; This says that the uncertainty in this subtracted
               ; quantity is about 10%.
               invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                        srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
                 1.0/(abs(image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                        srow[iamp]:srow[iamp]+nrow[iamp]-1]) + $
                      abs(smearimg[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                        srow[iamp]:srow[iamp]+nrow[iamp]-1]) + $
                      0.01 * abs(smearimg[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                        srow[iamp]:srow[iamp]+nrow[iamp]-1])^2 + $
                        (readnoiseDN[iamp]*gain[iamp])^2) + $ 
                          1.0e-4  ; floor to variance
            endif else begin
               invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                        srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
                 1.0/(abs(image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                        srow[iamp]:srow[iamp]+nrow[iamp]-1]) + $
                        (readnoiseDN[iamp]*gain[iamp])^2) + $
                        1.0e-4  ; floor to variance
            endelse

         endif
      endfor

      ;------
      ; Look for blead trails and mask them.
      ; At present, SATMASK is set to 0 for saturated pixels.  Look for any
      ; column with >=11 saturated pixels in a row.  In that case, mask
      ; all pixels from that row number up to the top of the CCD.

      kern = transpose(fltarr(11) + 1.0)
      mask1 = fix( convol(satmask+0.0, kern, /center, /edge_truncate) ) EQ 0
         ; 1=bad
      qblead = total(mask1, 2) GT 0 ; =1 for each column with a blead trail
      iblead = where(qblead, nblead)
      sxaddpar, hdr, 'NBLEAD', nblead, ' Number of columns with blead trails'
      if (nblead GT 0) then begin
         splog, 'Number of bleading columns = ', nblead
         for i=0, nblead-1 do begin
            icol = iblead[i] ; Column number for this blead trail
            irow = (where(mask1[icol,*]))[0] ; First bad row in this column
            satmask[icol,irow:nr-1] = 0
         endfor
      endif

      ;------
      ; Count the number of rows (wavelengths) that we think are saturated.
      ; A row is considered to be saturated if at least 20 of its pixels are.
      ; Note that this counting is done before masking out bad columns.
      ; What we really should do is ignore bad columns when computing this.

      if (arg_present(nsatrow)) then begin
         totsat = total((1-satmask), 1)
         junk = where(totsat GE 20, nsatrow)
         splog, 'Number of saturated rows = ', nsatrow
      endif

      ;------
      ; Mask out pixels that saturated the A/D converter, plus mask
      ; all neighbors within 1 pixel

      ngrow = 1
      width = 2*ngrow + 1
      admask = smooth(admask * width^2, width) GT 0 ; 1=bad

      ;------
      ; Mask out bad columns

      if (nbc GT 0) then begin
         bcsc = (bc.dfcol0 > 0) < (nc-1)
         bcec = (bc.dfcol0 + bc.dfncol - 1 < (nc-1)) > bcsc
         bcsr = (bc.dfrow0 > 0) < (nr-1)
         bcer = (bc.dfrow0 + bc.dfnrow - 1 < (nr-1)) > bcsr

         for i=0, nbc-1 do bcmask[bcsc[i]:bcec[i],bcsr[i]:bcer[i]] = 1
      endif

      ;------
      ; For masked pixels, set INVVAR=0

      invvar = invvar * satmask * (1-admask) * (1-bcmask)

      ;------
      ; Count the fraction of bad pixels, not including bad columns

      if (arg_present(fbadpix)) then begin
         junk = where((satmask EQ 0 OR admask EQ 1) AND (bcmask EQ 0), njunk)
         fbadpix = float(njunk) / (float(nc) * float(nr))
      endif

satmask = 0
admask = 0
bcmask = 0
   endif

   ;---------------------------------------------------------------------------
   ; Read pixel-to-pixel flat-field
   ;---------------------------------------------------------------------------

   if (keyword_set(pixflatname) AND (readimg OR readivar)) then begin
      splog, 'Correcting with pixel flat ' + pixflatname
      fullname = findfile(pixflatname, count=ct)
      if (ct EQ 0) then $
       message, 'Cannot find pixflat image ' + pixflatname
      pixflat = readfits(fullname[0])

      if (readimg) then image = image / pixflat
      if (NOT keyword_set(minflat)) then minflat = 0.0
      if (readivar) then invvar = invvar * pixflat^2 * (pixflat GT minflat)
pixflat = 0
   endif

   ;---------------------------------------------------------------------------
   ; Check for NaN's
   ;---------------------------------------------------------------------------

   ; This should never happen, but just in case...

   if (readimg OR readivar) then begin
      inan = where(finite(image) EQ 0, nnan)
      if (nnan GT 0) then begin
         splog, 'WARNING: Replacing ', nnan, ' NaN values'
         image[inan] = 0
         invvar[inan] = 0
      endif
   endif

   ;---------------------------------------------------------------------------
   ; Write output files
   ;---------------------------------------------------------------------------

   if (keyword_set(outfile)) then begin
      if (keyword_set(varfile)) then $
       sxaddpar, hdr, 'VARFILE', varfile, ' Corresponding inverse var file'
      writefits, outfile, image, hdr
   endif

   if (readivar) then begin
      varhdr = hdr
      if (keyword_set(outfile)) then $
       sxaddpar, hdr, 'IMGFILE', outfile, ' Corresponding image file'
      if (keyword_set(varfile)) then $
       writefits, varfile, invvar, varhdr
   endif
 
   return
end
;------------------------------------------------------------------------------
