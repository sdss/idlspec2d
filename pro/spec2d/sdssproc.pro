;+
; NAME:
;   sdssproc
;
; PURPOSE:
;   Read in raw SDSS files, and process with opConfig, opECalib, opBC par files.
;
; CALLING SEQUENCE:
;   sdssproc, infile, [image, invvar, indir=, $
;    outfile=, varfile=, nsatrow=, fbadpix=, $
;    hdr=hdr, configfile=, ecalibfile=, bcfile=, $
;    /applybias, /applypixflat, minflat=, maxflat=, $
;    spectrographid=, color=, camname= ]
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
;   configfile - Default to "opConfig*par", selecting the file with the
;                appropriate MJD.
;   ecalibfile - Default to "opECalib*par", selecting the file with the
;                appropriate MJD.
;   bcfile     - Default to "opBC*par", selecting the file with the
;                appropriate MJD.
;   applybias  - Apply 2-D bias image.
;   applypixflat- Apply 2-D pixel-to-pixel flat (after subtracting bias).
;   minflat    - Minimum values allowed for pixflat; pixels with the
;                flat out of range are masked; default to 0.
;   maxflat    - Maximum values allowed for pixflat; pixels with the
;                flat out of range are masked; default to 1e10.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   image      - Processed 2d image
;   invvar     - Associated inverse variance
;   hdr        - Processed FITS header
;   spectrographid - Return spectrograph ID (1 or 2)
;   color      - Return spectrograph color ('red' or 'blue')
;   camname    - Return camera name: 'b1', 'r1', 'b2', or 'r2'
;
; COMMENTS:
;   Only the header is read from the image if IMAGE, INVVAR, OUTFILE and
;   VARFILE are all not set.
;
;   Required header keywords: EXPTIME.
;
;   The returned image is in electrons, not ADU.
;
;   The signal-to-noise is limited to never exceed 100, by adding 1.e-4
;   times the flux to the variance term.
;
;   Change the CAMERAS keyword to the camera as specified by the file name.
;
;   Rename 'target' to 'science', and 'calibration' to 'arc' in the
;   header keyword FLAVOR.
;
;   Determine the exposure number from the file name itself.
;
; BUGS:
;   The open-shutter correction SMEARIMG will include smeared data from
;   any cosmic rays, which is wrong.  At the minimum, I could interpolate
;   over A/D saturations (in ADMASK) before constructing SMEARIMG.
;
; PROCEDURES CALLED:
;   djs_filepath()
;   djs_iterstat
;   fileandpath()
;   fits_wait()
;   headfits()
;   idlspec2d_version()
;   idlutils_version()
;   lookforgzip()
;   rdss_fits()
;   readfits()
;   sphdrfix
;   splog
;   sxaddpar
;   sxpar()
;   writefits
;   yanny_free
;   yanny_par()
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;   findopfile()
;   make_badcolumn_mask()
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/opConfig*par
;   $IDLSPEC2D_DIR/examples/opECalib*par
;   $IDLSPEC2D_DIR/examples/opBC*par
;   $SPECFLAT_DIR/biases/pixbias*.fits
;   $SPECFLAT_DIR/flats/pixflat*.fits
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
;   26-Jul-2000  Added fix for "dropped pixel" problem for data on or after
;                MJD=51688 (23 May 2000).  Should disable this code for later
;                MJD's once this problem is fixed in the electronics.
;   26-Jul-2000  Added fix for more severe "shifted row" electronics problem
;                for data taken on MJD=51578 to 51580 (3 to 5 Feb 2000).
;   26-Aug-2000  Horrible kludge for this night MJD=51779 (22/23 Aug 2000).
;                Add a noise term of 100 ADU to the left amplifier of r2.
;                That amplifier had random bits being set incorrectly,
;                in particular the 32 bit and 256 bit.
;   04-Nov-2000  Measure the bias values using DJS_ITERSTAT instead of MEDIAN,
;                since the median is always only an integer value.
;   31-Jan-2001  Determine the exposure number from the file name itself,
;                since the counting got off by one on MJD=51882.
;-
;------------------------------------------------------------------------------
function findopfile, expres, mjd, indir

   files = findfile(filepath(expres, root_dir=indir), count=nfile)
   if (nfile EQ 0) then $
    message, 'Cannot find opFile '+expres

   mjdlist = lonarr(nfile)
   for i=0,nfile-1 do begin
      thisfile = fileandpath(files[i])
      if (strmid(thisfile, strpos(thisfile, '.')) EQ '.par') then begin
         ; Get the MJD from the header of the Yanny file
         yanny_read, files[i], hdr=hdr
         mjdlist[i] = long(yanny_par(hdr, 'mjd'))
      endif else begin
         ; Get the MJD from the file name; use the 5 digits after the first minus sign
         mjdlist[i] = strmid(thisfile, strpos(thisfile,'-')+1,5)
      endelse
   endfor

   diff = mjd - mjdlist
   score = min(diff  + 10000000L*(diff LT 0), bestmjd)

   return, fileandpath(files[bestmjd])
end
;------------------------------------------------------------------------------
;  Create the bad column mask (1 for a masked pixel) with image size nc,nr
;  If the operation doesn't work just return 0 for no masked pixels
;
function make_badcolumn_mask, bcfile, camrow, camcol, nc=nc, nr=nr

   if NOT keyword_set(nc) then nc=2048L
   if NOT keyword_set(nr) then nr=2048L

   yanny_read, bcfile, pdata
   if (size(pdata,/tname)) EQ 'INT' then begin
     splog,' Could not read BC file '+bcfile
     return, 0
   endif

   bc = *pdata[0]
   yanny_free, pdata

   bchere = where(bc.camrow EQ camrow AND bc.camcol EQ camcol, nbc)
   if NOT keyword_set(nbc) then begin
     splog,' Could not parse BC file contents or find corresponding camera' $
               +bcfile
     return, 0
   endif

   bc = bc[ bchere ]
   bcmask = bytarr(nc, nr)

   ; Mask out bad columns

   if (nbc GT 0) then begin
     bcsc = (bc.dfcol0 > 0) < (nc-1)
     bcec = (bc.dfcol0 + bc.dfncol - 1 < (nc-1)) > bcsc
     bcsr = (bc.dfrow0 > 0) < (nr-1)
     bcer = (bc.dfrow0 + bc.dfnrow - 1 < (nr-1)) > bcsr

     for i=0, nbc-1 do bcmask[bcsc[i]:bcec[i],bcsr[i]:bcer[i]] = 1
   endif

   return, bcmask
end

;------------------------------------------------------------------------------
pro sdssproc, infile, image, invvar, indir=indir, $
 outfile=outfile, varfile=varfile, nsatrow=nsatrow, fbadpix=fbadpix, $
 hdr=hdr, configfile=configfile, ecalibfile=ecalibfile, bcfile=bcfile, $
 applybias=applybias, applypixflat=applypixflat, $
 minflat=minflat, maxflat=maxflat, $
 spectrographid=spectrographid, color=color, camname=camname

   if (N_params() LT 1) then begin
      print, 'Syntax - sdssproc, infile, [image, invvar, indir=, $'
      print, ' outfile=, varfile=, nsatrow=, fbadpix=, $' 
      print, ' hdr=, configfile=, ecalibfile=, bcfile=, $'
      print, ' /applybias, /applypixflat, minflat=, maxflat=, $'
      print, ' spectrographid=, color=, camname= ]'
      return
   endif

   readimg = arg_present(image) OR keyword_set(outfile)
   readivar = arg_present(invvar) OR keyword_set(varfile) $
    OR arg_present(nsatrow) OR arg_present(fbadpix)

   fullname = djs_filepath(infile, root_dir=indir)
   fullname = (lookforgzip(fullname, count=ct))[0]
;   fullname = (findfile(fullname, count=ct))[0]
   if (ct NE 1) then $
    message, 'Cannot find image ' + infile

   if (readimg OR readivar) then $
    rawdata = rdss_fits(fullname, hdr, /nofloat) $
   else $
    hdr = headfits(fullname)

   ;-----------
   ; Fix the headers with any hand-edits that we have determined.

   if (!version.release LT '5.3') then $
    splog, 'Warning: Unable to fix headers with this version of IDL' $
   else $
    sphdrfix, infile, hdr

   ;-----------
   ; Determine the exposure number from the file name itself.
   ; Very bad form, but this information is sometimes wrong in the header.
   ; In particular, the counting got off by 1 on MJD=51882.

   i = rstrpos(infile, '-')
   expnum = long( strmid(infile, i+1, 8) )
   if (NOT keyword_set(expnum)) then $
    message, 'Cannot determine exposure number from file name ' + infile
   hdrexp = sxpar(hdr, 'EXPOSURE')
   if (expnum NE hdrexp) then begin
      splog, 'WARNING: Exposure number in header (' $
       + strtrim(string(hdrexp),2) + ') disagrees w/filename (' $
       + strtrim(string(expnum),2) + ') !!'
      sxaddpar, hdr, 'EXPOSURE', expnum
   endif

   ;-----------
   ; Determine which CCD from the file name itself, using either the
   ; numbering scheme (01,02,03,04) or naming scheme (b1,r2,b2,r1).
   ; Very bad form, but this information is not in the header since
   ; the CAMERAS keyword is sometimes wrong.

   i = rstrpos(infile, '-')
   if (i[0] EQ -1 OR i-2 LT 0) then $
    message, 'Cannot determine CCD number from file name ' + infile

   camnames = ['b1', 'r2', 'b2', 'r1']
   camnums = ['01', '02', '03', '04']

   ; First try to match a camera name (e.g., 'b1'), then try to match
   ; a camera number (e.g., '01').  If both fail, then abort.
   indx = where(strmid(infile, i-2, 2) EQ camnames, ct)
   if (ct NE 1) then $
     indx = where(strmid(infile, i-2, 2) EQ camnums, ct)
   if (ct NE 1) then $
    message, 'Cannot determine CCD number from file name ' + infile

; Do not read the camera from the CAMERAS keyword, since this was often
; wrong in the early days!
;   cameras = strtrim( sxpar(hdr, 'CAMERAS'), 2 )
   camname = camnames[indx[0]]
   case camname of
    'b1': begin
          spectrographid = 1
          color = 'blue'
          end
    'r1': begin
          spectrographid = 1
          color = 'red'
          end
    'b2': begin
          spectrographid = 2
          color = 'blue'
          end
    'r2': begin
          spectrographid = 2
          color = 'red'
          end
   endcase
   camcol = indx[0] + 1
   camrow = 0

   if (!version.release LT '5.4') then $
    spawn, 'speclog_version', verslog $
   else $
    spawn, 'speclog_version', verslog, errstring
   if (NOT keyword_set(verslog)) then verslog = 'Unknown'

   if (!version.release LT '5.4') then $
    spawn, 'specflat_version', versflat $
   else $
    spawn, 'specflat_version', versflat, errstring
   if (NOT keyword_set(versflat)) then versflat = 'Unknown'

   sxaddpar, hdr, 'CAMROW', camrow
   sxaddpar, hdr, 'CAMCOL', camcol
   sxaddpar, hdr, 'TELESCOP', 'SDSS 2.5-M', ' Sloan Digital Sky Survey'
   sxaddpar, hdr, 'AUTHOR', 'Scott Burles & David Schlegel'
   sxaddpar, hdr, 'VERSIDL', !version.release, ' Version of IDL'
   sxaddpar, hdr, 'VERSUTIL', idlutils_version(), ' Version of idlutils'
   sxaddpar, hdr, 'VERSREAD', idlspec2d_version(), $
    ' Version of idlspec2d for pre-processing raw data', after='VERSUTIL'
   sxaddpar, hdr, 'VERSLOG', verslog[0], ' Version of SPECLOG product', $
    after='VERSREAD'
   sxaddpar, hdr, 'VERSFLAT', versflat[0], ' Version of SPECFLAT product', $
    after='VERSFLAT'
 
   ;-----------
   ; Rename 'target' -> 'science', and 'calibration' -> 'arc'

   mjd = sxpar(hdr, 'MJD')
   flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)
   if (flavor EQ 'target') then flavor = 'science'
   if (flavor EQ 'calibration') then flavor = 'arc'
   if (mjd GT 51576) then begin
      if (sxpar(hdr, 'COLBIN') NE 1 OR $
          sxpar(hdr, 'ROWBIN') NE 1) then flavor = 'unknown'
   endif
   sxaddpar, hdr, 'FLAVOR', flavor
   sxaddpar, hdr, 'CAMERAS', camname

   if (NOT readimg AND NOT readivar) then return

   ;-----------
   ; Fix the shifted-row problem in the electronics that appeared
   ; on MJD=51578 to 51580 (3-5 Feb 2000) for spectrograph-2.
   ; Note that the bad rows need to be identified on the red frame,
   ; so if we are reducing b2, we need to read in the image for r2
   ; to identify the bad rows.

   if ((mjd GE 51578 AND mjd LE 51580) $
    AND (camname EQ 'b2' OR camname EQ 'r2') $
    AND (readimg OR readivar)) then begin
      if (camname EQ 'b2') then begin
;         i1 = strpos(infile,'b2',/reverse_search) ; IDL 5.3 command
         i1 = rstrpos(infile,'b2')
         if (i1 EQ -1) then $
          message, 'Unable to parse corresponding red file for '+infile
         redfile = infile
         strput, redfile, 'r2', i1
         redfile = (lookforgzip(djs_filepath(redfile, root_dir=indir)))[0]

         if (fits_wait(redfile, deltat=1, tmax=1)) then $
          reddata = rdss_fits(redfile, /nofloat) $
         else $
          splog, 'Warning: Could not read corresponding red file ' + redfile
      endif else begin
         reddata = rawdata
      endelse

      if (keyword_set(reddata)) then $
       ibad = where( reddata[20,*] LT median(reddata[20,*]) - 100 , nbad ) $
      else $
       nbad = 0

      if (nbad GT 0) then begin
         splog, 'WARNING: Fixing ', nbad, ' shifted rows (from electronics)'

         ; For unrecoverable data, set KILLDATA=0
         killdata = byte(0*rawdata) + 1b

         medval = median(rawdata[22:39,*])
         for ii=0, nbad-1 do begin
            nshift = $
             (reverse(where(rawdata[20:39,ibad[ii]] GT medval + 100)))[0]

            ; If nshift is -1, we will also zero the whole row

            if (nshift GE 0 AND nshift LT 10 ) then begin
                rawdata[0:1063-nshift,ibad[ii]] = rawdata[nshift:1063,ibad[ii]]
                killdata[1063-nshift+1:1063,ibad[ii]] = 0
            endif else begin
                killdata[*,ibad[ii]] = 0
            endelse

         endfor
      endif

      reddata = 0  ; Free memory
   endif

   ;-----------
   ; Fix the "dropped pixel" problem in the electronics that appeared
   ; on MJD=51688 (23 May 2000) for spectrograph-2.
   ; Note that the bad rows need to be identified on the red frame,
   ; so if we are reducing b2, we need to read in the image for r2
   ; to identify the bad rows.
   ; Reference e-mail discussion with JEG on 01-Jun-2000.

   if ((mjd GE 51688) $
    AND (camname EQ 'b2' OR camname EQ 'r2') $
    AND (readimg OR readivar)) then begin
      if (camname EQ 'b2') then begin
;         i1 = strpos(infile,'b2',/reverse_search) ; IDL 5.3 command
         i1 = rstrpos(infile,'b2') ; IDL 5.2 command
         if (i1 EQ -1) then $
          message, 'Unable to parse corresponding red file for '+infile
         redfile = infile
         strput, redfile, 'r2', i1
         redfile = (lookforgzip(djs_filepath(redfile, root_dir=indir)))[0]

         if (fits_wait(redfile, deltat=1, tmax=1)) then $
          reddata = rdss_fits(redfile, /nofloat) $
         else $
          splog, 'Warning: Could not read corresponding red file ' + redfile
      endif else begin
         reddata = rawdata
      endelse

      if (keyword_set(reddata)) then $
       ibad = where( reddata[20,*] LT median(reddata[20,*]) - 100 $
        AND reddata[2107,*] GT median(reddata[2107,*]) + 100, nbad ) $
      else $
       nbad = 0

      if (nbad GT 0) then begin
         splog, 'WARNING: Fixing ', nbad, ' dropped-pixel rows (from electronics)'
         rawdata[1:1063,ibad] = rawdata[0:1062,ibad]
         rawdata[20,ibad] = median( rawdata[20,ibad] )
         rawdata[2108:2127,ibad] = rawdata[2107:2126,ibad]
         rawdata[2107,ibad] = median( rawdata[2107,ibad] )
      endif

      reddata = 0  ; Free memory
   endif

   ;-----------
   ; Find names of the configurations files

   config_dir = filepath('', $
      root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')

   if (NOT keyword_set(configfile)) then $
       configfile = findopfile('opConfig*par',mjd,config_dir)
   if (NOT keyword_set(ecalibfile)) then $
       ecalibfile = findopfile('opECalib*par',mjd,config_dir)
   if (NOT keyword_set(bcfile)) then $
       bcfile = findopfile('opBC*par',mjd,config_dir)

   naxis1 = sxpar(hdr,'NAXIS1')
   naxis2 = sxpar(hdr,'NAXIS2')
   if (naxis1 NE 2128 OR naxis2 NE 2069) then $
    splog, 'WARNING: Expecting 2128x2069, found '+string(naxis1)+'x'$
     +string(naxis2)

   ;------
   ; Read in opConfig.par file
   ; Take the first entry for the configuration of each CCD in the event
   ; that there are several.

   yanny_read, filepath(configfile, root_dir=config_dir), pdata
   config = *pdata[0]
   yanny_free, pdata
   i = where(config.camrow EQ camrow AND config.camcol EQ camcol)
   config = config[i[0]]

   if (naxis1 NE config.ncols OR naxis2 NE config.nrows) then $
      splog, 'WARNING! Config file dimensions do not match raw image'

   qexist = [config.amp0, config.amp1, config.amp2, config.amp3]

   ;------
   ; Define the "overscan" regions

   sover = [config.soverscan0, config.soverscan1, config.soverscan2, $
    config.soverscan3]
   nover = [config.noverscan0, config.noverscan1, config.noverscan2, $
    config.noverscan3]

   ;------
   ; Define the "mapped overscan" regions

   smapover = [config.smapoverscan0, config.smapoverscan1, $
    config.smapoverscan2, config.smapoverscan3]
   nmapover = [config.nmapoverscan0, config.nmapoverscan1, $
    config.nmapoverscan2, config.nmapoverscan3]

   ;------
   ; Define the "overscan rows" (at the bottom of the CCD)

   soverrow = [config.soverscanrows0, config.soverscanrows1, $
    config.soverscanrows2, config.soverscanrows3]
   noverrow = [config.noverscanrows0, config.noverscanrows1, $
    config.noverscanrows2, config.noverscanrows3]

   ;------
   ; Data position in the original image

   sdatarow = [config.sdatarow0, config.sdatarow1, $
    config.sdatarow2, config.sdatarow3]
   sdatacol = [config.sdatasec0, config.sdatasec1, config.sdatasec2, $
    config.sdatasec3]
   nrow = [config.ndatarow0, config.ndatarow1, $
    config.ndatarow2, config.ndatarow3]
   ncol = [config.ndatasec0, config.ndatasec1, config.ndatasec2, $
    config.ndatasec3]

   ;------
   ; Data position in the final (trimmed) image

   srow = [config.sccdrowsec0, config.sccdrowsec1, $
    config.sccdrowsec2, config.sccdrowsec3]
   scol = [config.sccdcolsec0, config.sccdcolsec1, $
    config.sccdcolsec2, config.sccdcolsec3]

   if (naxis2 EQ 2049) then begin
     splog, 'WARNING: NROWS is 2049, adjusting config entries' 
     sdatarow = sdatarow - 20
     noverrow = 1
   endif

   ;------
   ; Read in ECalib File

   yanny_read, filepath(ecalibfile, root_dir=config_dir), pdata
   ecalib = *pdata[0]
   yanny_free, pdata
   ecalib = ecalib[ where(ecalib.camrow EQ camrow AND ecalib.camcol EQ camcol) ]

   gain = [ecalib.gain0, ecalib.gain1, ecalib.gain2, ecalib.gain3]
   readnoiseDN = [ecalib.readnoiseDN0, ecalib.readnoiseDN1, $
    ecalib.readnoiseDN2, ecalib.readnoiseDN3]
   fullWellDN = [ecalib.fullWellDN0, ecalib.fullWellDN1, $
    ecalib.fullWellDN2, ecalib.fullWellDN3]

   ;------
   ; Construct the final image

   igood = where(qexist)
   nr = max((srow+nrow)[igood])
   nc = max((scol+ncol)[igood])
   if ((size(image))[0] NE 2) then image = fltarr(nc, nr) $
   else if ((size(image))[1] NE nc OR (size(image))[2] NE nr OR $
            (size(image))[3] NE 4) then image = fltarr(nc, nr) 

   ;------
   ; Moved the Bad Column mask creation to the function above
   ;
   ; yanny_read, filepath(bcfile, root_dir=config_dir), pdata
   ;  bc = *pdata[0]
   ;  yanny_free, pdata
   ;  bchere = where(bc.camrow EQ camrow AND bc.camcol EQ camcol, nbc)
   ;  if (nbc GT 0) then bc = bc[ bchere ]

   ;------
   ; Test to see if the shutter was open during readout if the exposure
   ; was longer than 640 seconds.

   exptime = sxpar(hdr, 'EXPTIME')
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

            djs_iterstat, biasreg, sigrej=3.0, mean=biasval, sigma=readoutDN

            splog, 'Measured read-noise in DN for amp#', iamp, ' = ', readoutDN
            splog, 'Measured bias value in DN for amp#', iamp, ' = ', biasval
            splog, 'Applying gain for amp#', iamp, ' = ', gain[iamp]

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

            ; Flux term below
            expr1 = abs(image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                     srow[iamp]:srow[iamp]+nrow[iamp]-1])
            ; Add to the variance image from the open shutter
            ; by adding 1% of that signal^2 to the variance.
            ; This says that the uncertainty in this subtracted
            ; quantity is about 10%.
            if (qshutter) then $
               expr2 = abs(smearimg[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                        srow[iamp]:srow[iamp]+nrow[iamp]-1]) $
            else expr2 = 0
            ; Read noise term below
            expr3 = (readnoiseDN[iamp]*gain[iamp])^2
            ; Term below to limit best S/N to under 100
            expr4 = 1.e-4 * expr1^2

            ; Horrible kludge for this night MJD=51779 (22/23 Aug 2000).
            ; Add a noise term of 100 ADU to the left amplifier of r2
            if (mjd EQ 51779 AND camname EQ 'r2' AND iamp EQ 2) then $
             expr4 = expr4 + 100.^2

            invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                     srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
              1.0/(expr1 + expr2 + 0.01 * expr2^2 + expr3 + expr4)

            if (keyword_set(killdata)) then $
             invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                     srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
              invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                      srow[iamp]:srow[iamp]+nrow[iamp]-1] * $
              killdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
                     sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1]
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

      bcmask = make_badcolumn_mask($
          filepath(bcfile,root_dir=config_dir), camrow, camcol)

      ;------
      ; For masked pixels, set INVVAR=0

      invvar = invvar * satmask * (1-admask) * (1-bcmask)

      ;------
      ; Count the fraction of bad pixels, not including bad columns

      if (arg_present(fbadpix)) then begin
         junk = where((satmask EQ 0 OR admask EQ 1) AND (bcmask EQ 0), njunk)
         fbadpix = float(njunk) / (float(nc) * float(nr))
      endif

satmask = 0 ; clear memory
admask = 0 ; clear memory
; bcmask = 0 ; clear memory
   endif

   ;---------------------------------------------------------------------------
   ; Correct image with bias image
   ;---------------------------------------------------------------------------

   if (keyword_set(applybias) AND readimg) then begin
      pp = filepath('', root_dir=getenv('SPECFLAT_DIR'), subdirectory='biases')
      ; First search for files "pixbiasave-*.fits", and if not found then
      ; look for "pixbias-*.fits".
      pixbiasname = findopfile('pixbiasave-*-'+camname+'.fits', mjd, pp)
      if (NOT keyword_set(pixbiasname)) then $
       pixbiasname = findopfile('pixbias-*-'+camname+'.fits', mjd, pp)

      if (NOT keyword_set(pixbiasname)) then begin
         splog, 'WARNING: Bias image not found for camera ' + camname
      endif else begin
         splog, 'Correcting with bias image ' + pixbiasname
         pixbiasimg = mrdfits(djs_filepath(pixbiasname, root_dir=pp))

         image = image - pixbiasimg
pixbiasimg = 0 ; clear memory

         ; Add pixbiasname to header since it has just been applied
         sxaddpar, hdr, 'BIASIMG', pixbiasname
      endelse
   endif

   ;---------------------------------------------------------------------------
   ; Correct image with pixel-to-pixel flat-field
   ;---------------------------------------------------------------------------

   if (keyword_set(applypixflat) AND (readimg OR readivar)) then begin
      pp = filepath('', root_dir=getenv('SPECFLAT_DIR'), subdirectory='flats')
      ; First search for files "pixflatave-*.fits", and if not found then
      ; look for "pixflat-*.fits".
      pixflatname = findopfile('pixflatave-*-'+camname+'.fits', mjd, pp)
      if (NOT keyword_set(pixflatname)) then $
       pixflatname = findopfile('pixflat-*-'+camname+'.fits', mjd, pp)

      if (NOT keyword_set(pixflatname)) then begin
         splog, 'WARNING: Pixel flat not found for camera ' + camname
      endif else begin
         splog, 'Correcting with pixel flat ' + pixflatname

         pixflatimg = mrdfits(djs_filepath(pixflatname, root_dir=pp))

         if (readimg) then image = image / (pixflatimg + (pixflatimg LE 0))
         if (NOT keyword_set(minflat)) then minflat = 0.0
         if (NOT keyword_set(maxflat)) then maxflat = 1.0e10
         if (readivar) then $
          invvar = invvar * pixflatimg^2 * (pixflatimg GT minflat) $
           * (pixflatimg LT maxflat)
pixflatimg = 0 ; clear memory

         ; Add pixflatname to header since it has just been applied
         sxaddpar, hdr, 'PIXFLAT', pixflatname
      endelse
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

   if keyword_set(bcfile) then sxaddpar, hdr, 'OPBC', bcfile
   if keyword_set(configfile) then sxaddpar, hdr, 'OPCONFIG', configfile
   if keyword_set(ecalibfile) then sxaddpar, hdr, 'OPECALIB', ecalibfile
   sxdelpar, hdr, 'UNSIGNED'   

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
