;+
; NAME:
;   sdssproc
;
; PURPOSE:
;   Read in Raw SDSS files, and process with opECalib.par and opConfig.par
;
; CALLING SEQUENCE:
;   sdssproc, infile, [image, invvar, indir=indir, $
;    outfile=outfile, varfile=varfile, $
;    hdr=hdr, configfile=configfile, ecalibfile=ecalibfile, bcfile=bcfile, $
;    pixflatname=pixflatname, spectrographid=spectrographid, color=color ]
;
; INPUTS:
;   infile     - Raw SDSS frame
;
; OPTIONAL KEYWORDS:
;   indir      - Input directory for INFILE
;   outfile    - Calibrated 2d frame, after processing
;   varfile    - Inverse Variance Frame after processing
;   hdr        - Header returned in memory
;   configfile - Default to "opConfig.par"
;   ecalibfile - Default to "opECalib.par"
;   bcfile     - Default to "opBC.par"
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
; BUGS:
;
; PROCEDURES CALLED:
;   djs_iterstat
;   headfits()
;   idlspec2d_version()
;   idlutils_version()
;   rdss_fits()
;   readfits()
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
;                Identify blead trails and mask from that row up. (DJS)
;-
;------------------------------------------------------------------------------

pro sdssproc, infile, image, invvar, indir=indir, $
 outfile=outfile, varfile=varfile, $
 hdr=hdr, configfile=configfile, ecalibfile=ecalibfile, bcfile=bcfile, $
 pixflatname=pixflatname, spectrographid=spectrographid, color=color

   if (N_params() LT 1) then begin
      print, 'Syntax - sdssproc, infile, [image, invvar, indir=indir, '
      print, ' outfile=outfile, varfile=varfile, ' 
      print, ' hdr=hdr, configfile=configfile, ecalibfile=ecalibfile, bcfile=bcfile]'
      return
   endif
   if (NOT keyword_set(configfile)) then configfile = 'opConfig.par'
   if (NOT keyword_set(ecalibfile)) then ecalibfile = 'opECalib.par'
   if (NOT keyword_set(bcfile)) then bcfile = 'opBC.par'

   readimg = arg_present(image) OR keyword_set(outfile)
   readivar = arg_present(invvar) OR keyword_set(varfile)

   pp = getenv('IDLSPEC2D_DIR')+'/examples'

   junk = findfile(configfile, count=ct)
   if (ct NE 1) then begin
     tempname = findfile(filepath(configfile, root_dir=pp), count=ct)
   endif
   if (ct NE 1) then $
     message, 'No configuration file ' + string(configfile)

   realconfig = tempname[0]

   tempname = findfile(ecalibfile, count=ct)
   if (ct NE 1) then begin
     tempname = findfile(filepath(ecalibfile, root_dir=pp), count=ct)
   endif
   if (ct NE 1) then $
    message, 'No ECalib file ' + string(ecalibfile)

   realecalib = tempname[0]

   tempname = findfile(bcfile, count=ct)
   if (ct NE 1) then begin
     tempname = findfile(filepath(bcfile, root_dir=pp), count=ct)
   endif
   if (ct NE 1) then $
    message, 'No BC file ' + string(bcfile)

   realbc = tempname[0]

   if (keyword_set(indir)) then inpath = filepath(infile, root_dir=indir) $
    else inpath = infile
   fullname = (findfile(inpath, count=ct))[0]
   if (ct NE 1) then $
    message, 'Cannot find image ' + infile

   if (readimg OR readivar) then $
    rawdata = rdss_fits(fullname, hdr, /nofloat) $
   else $
    hdr = headfits(fullname)

   cards = sxpar(hdr,'NAXIS*')
;   if (cards[0] NE 2128 OR cards[1] NE 2069) then $
;      message, 'Expecting 2128x2069, found '+string(cards[0])+','$
;               +string(cards[1])

   ; Determine which CCD from the file name itself!
   ; Very bad form, but this information is not in the header.
   ; The CAMERAS keyword is sometimes for the wrong camera.
   i = rstrpos(infile, '-')
   if (i[0] EQ -1 OR i-2 LT 0) then $
    message, 'Cannot determine CCD number from file name ' + infile

   camnames = ['b1', 'r2', 'b2', 'r1']
   camnums = ['01', '02', '03', '04']

   ; They've changed filenames again, this works both ways

   camplace = where(strmid(infile, i-2, 2) EQ camnames, camct)
   if (camct NE 1) then $
     camplace = where(strmid(infile, i-2, 2) EQ camnums, camct)
   if (camct NE 1) then  message, 'do not know what camera this is'

   camcol = camplace[0] + 1
     
   cameras = strtrim( sxpar(hdr, 'CAMERAS'), 2 )
   case camcol of
     1: begin
        spectrographid = 1
        color = 'blue'
         end
     4: begin
        spectrographid = 1
        color = 'red'
         end
     3: begin
        spectrographid = 2
        color = 'blue'
         end
     2: begin
        spectrographid = 2
        color = 'red'
        end
     else: begin
        print, 'CAMERAS keyword not found, guessing b2'
        spectrographid = 2
        color = 'blue'
        camcol = 3
        sxaddpar, hdr, 'CAMERAS', 'b2       ', ' Guessed b2 by default'
;        message, 'Cannot determine CCD number from file name ' + infile
        end
   endcase
   camrow = 0

   sxaddpar, hdr, 'CAMROW', camrow
   sxaddpar, hdr, 'CAMCOL', camcol
   sxaddpar, hdr, 'TELESCOP', 'SDSS 2.5-M', ' Sloan Digital Sky Survey'
   sxaddpar, hdr, 'AUTHOR', 'Scott Burles & David Schlegel'
   sxaddpar, hdr, 'SPEC2D_V', idlspec2d_version(), ' Version of idlspec2d'
   sxaddpar, hdr, 'UTILS_V', idlutils_version(), ' Version of idlutils'
 
   ; Read in opConfig.par file

   yanny_read, realconfig, pdata
   config = *pdata[0]
   yanny_free, pdata
   config = config[ where(config.camrow EQ camrow AND config.camcol EQ camcol) ]

   if (cards[0] NE config.ncols OR cards[1] NE config.nrows) then $
      message, 'Config file dimensions do not match raw image'

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

   ; Read in ECalib File
   yanny_read, realecalib, pdata
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

   yanny_read, realbc, pdata
   bc = *pdata[0]
   yanny_free, pdata
   
   bchere = where(bc.camrow EQ camrow AND bc.camcol EQ camcol,nbc)
   if (nbc GT 0) then bc = bc[ bchere ]

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
           ; Now image is in electrons

           image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
            (rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
                  sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] - biasval) * gain[iamp]

         ; Add to the header
         sxaddpar, hdr, 'GAIN'+string(iamp,format='(i1)'), $
          gain[iamp], ' Gain in electrons per ADU'
         sxaddpar, hdr, 'RDNOISE'+string(iamp,format='(i1)'), $
          gain[iamp]*readnoiseDN[iamp], ' Readout noise in electrons'
        endif
      endif
   endfor

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

      satmask = bytarr(nc, nr)
      admask = bytarr(nc, nr)
 
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

         invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
           1.0/(abs(image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1]) + $
                  (readnoiseDN[iamp]*gain[iamp])^2)
         endif
      endfor

      ;------
      ; Look for blead trails and mask them.
      ; At present, SATMASK is set to 0 for saturated pixels.  Look for any
      ; column with >=11 saturated pixels in a row.

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
      ; Mask out bad columns

      if (nbc GT 0) then begin
         bcsc = (bc.dfcol0 > 0) < nc
         bcec = (bc.dfcol0 + bc.dfncol - 1 < nc) > bcsc
         bcsr = (bc.dfrow0 > 0) < nr
         bcer = (bc.dfrow0 + bc.dfnrow - 1 < nr) > bcsr

         for i=0,nbc-1 do satmask[bcsc[i]:bcec[i],bcsr[i]:bcer[i]] = 0
      endif

      ;------
      ; Mask out pixels that saturated the A/D converter, plus mask
      ; all neighbors within 1 pixel
      ngrow = 1
      width = 2*ngrow + 1
      admask = smooth(admask * width^2, width) GT 0 ; 1=bad

      ; For masked pixels, set INVVAR=0
      invvar = invvar * satmask * (1-admask)
   endif

   ;---------------------------------------------------------------------------
   ; Read pixel-to-pixel flat-field
   ;---------------------------------------------------------------------------

   if (keyword_set(pixflatname) AND (readimg OR readivar)) then begin
      fullname = findfile(pixflatname, count=ct)
      if (ct EQ 0) then $
       message, 'Cannot find pixflat image ' + pixflatname
      pixflat = readfits(fullname[0])

      if (readimg) then image = image / pixflat
      if (readivar) then invvar = invvar * pixflat^2
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
