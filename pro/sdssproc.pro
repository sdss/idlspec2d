;+
; NAME:
;   sdssproc
;
; PURPOSE:
;   Read in Raw SDSS files, and process with opECalib.par and opConfig.par
;
; CALLING SEQUENCE:
;   sdssproc, infile, [image, invvar, outfile=outfile, varfile=varfile, hdr=hdr,
;     configfile=configfile, ecalibfile=ecalibfile]
;
; INPUTS:
;   infile     - Raw SDSS frame
;
; OPTIONAL KEYWORDS:
;   outfile    - Calibrated 2d frame, after processing
;   varfile    - Inverse Variance Frame after processing
;   hdr        - Header returned in memory
;   configfile - Default to "opConfig.par"
;   ecalibfile - Default to "opECalib.par"
;   satvalue   - Saturation value; default to 1.0e6.
;                For all pixels with values above this, set INVVAR=0.
;
; OUTPUTS:
;   image      - Processed 2d image
;   invvar     - associated Inverse Variance
;
; COMMENTS:
;
; BUGS:
;
; PROCEDURES CALLED:
;   rdss_fits()
;   yanny_read
;
; REVISION HISTORY:
;   13-May-1999  Written by Scott Burles & David Schlegel, Apache Point.
;   08-Sep-1999  Modified to read Yanny param files instead of FITS
;                versions of the same (DJS).
;-
;------------------------------------------------------------------------------

pro sdssproc, infile, image, invvar, outfile=outfile, varfile=varfile, $
    hdr=hdr, configfile=configfile, ecalibfile=ecalibfile, satvalue=satvalue

   if (N_params() LT 1) then begin
      print, 'Syntax - sdssproc, infile, [image, invvar, outfile=outfile, varfile=varfile, ' 
      print, ' hdr=hdr, configfile=configfile, ecalibfile=ecalibfile]'
      return
   endif
   if (NOT keyword_set(configfile)) then configfile = 'opConfig.par'
   if (NOT keyword_set(ecalibfile)) then ecalibfile = 'opECalib.par'
   if (NOT keyword_set(satvalue)) then satvalue = 1.0e6
   
   junk = findfile(configfile, count=ct)
   if (ct NE 1) then begin
     pp = getenv('EVIL_PAR') 
     tempname = findfile(filepath(configfile, root_dir=pp), count=ct)
   endif
   if (ct NE 1) then $
     message, 'No configuration file ' + string(configfile)

   realconfig = tempname[0]

   tempname = findfile(ecalibfile, count=ct)
   if (ct NE 1) then begin
     pp = getenv('EVIL_PAR') 
     tempname = findfile(filepath(ecalibfile, root_dir=pp), count=ct)
   endif
   if (ct NE 1) then $
    message, 'No ECalib file ' + string(ecalibfile)

   realecalib = tempname[0]

   rawdata = rdss_fits(infile, hdr)

   cards = sxpar(hdr,'NAXIS*')
;   if (cards[0] NE 2128 OR cards[1] NE 2069) then $
;      message, 'Expecting 2128x2069, found '+string(cards[0])+','$
;               +string(cards[1])

   ; Test that the name of this file is consistent with an SDSS
   ; spectroscopic image.
   names = str_sep(infile,'/')
   camrow = long( strmid( names[N_elements(names) -1], 4,1 ))
   camcol = long( strmid( names[N_elements(names) -1], 5,1 ))

   if (camrow NE 0) then message, 'Not a spectroscopic image'
   if (camcol LT 1 OR camcol GT 4) then message, 'Not a spectroscopic image'
   sxaddpar, hdr, 'CAMROW', camrow
   sxaddpar, hdr, 'CAMCOL', camcol
 

   ; Read in opConfig.par file

   yanny_read, realconfig, pdata
   config = *pdata[0]
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
   image = fltarr(nc, nr)
   invvar = fltarr(nc, nr)
   mask = bytarr(nc, nr)

   for iamp=0, 3 do begin
      if (qexist[iamp] EQ 1) then begin
         if (nover[iamp] NE 0) then begin
            ; Use the "overscan" region
            biasval = median( $
             rawdata[sover[iamp]:sover[iamp]+nover[iamp]-1, $
             sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] )
         endif else if (nmapover[iamp] NE 0) then begin
            ; Use the "mapped overscan" region
            biasval = median( $
             rawdata[smapover[iamp]:smapover[iamp]+nmapover[iamp]-1, $
             sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] )
         endif

         ; Copy the data for this amplifier into the final image
         image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
         rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
                  sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] - biasval

         mask[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
           (rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
                  sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] LT $
                  fullWellDN[iamp])

         invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
           1.0/(abs(image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1]) /gain[iamp] + $
                  readnoiseDN[iamp]*readnoiseDN[iamp])

         ; Add to the header
         sxaddpar, hdr, 'BIAS'+string(iamp,format='(i1)'), biasval
         sxaddpar, hdr, 'GAIN'+string(iamp,format='(i1)'), gain[iamp]
         sxaddpar, hdr, 'RDNOISE'+string(iamp,format='(i1)'), $
          gain[iamp]*readnoiseDN[iamp]
      endif
   endfor

   if (keyword_set(outfile)) then $
    writefits, outfile, image, hdr

   if (keyword_set(varfile)) then begin
      varhdr = hdr
      sxaddpar, varhdr, 'LOOKHERE', 'INVERSE VARIANCE of ' + outfile
      writefits, varfile, invvar, varhdr
   endif

   ; For saturated pixels, set INVVAR=0
   invvar = invvar * mask

   return
end
