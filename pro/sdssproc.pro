;+
; Convert a raw SDSS spectroscopic image from unsigned 16-bit to float.
; Also, do a quick-and-dirty overscan-subtraction
; NAME:
;   sdssproc
;
; PURPOSE:
;   Read in Raw SDSS files, and process with opECalib and opConfig.fit
;
; CALLING SEQUENCE:
;   sdssproc, infile, outfile, [hdr=hdr, $
;     configfile=configfile, ccdfile=ccdfile, varfile=varfile]
;
; INPUTS:
;   infile     - Raw SDSS frame
;   outfile    - Calibrated 2d frame, after processing
;
; OPTIONAL KEYWORDS:
;   configfile - "opConfig.fit" 
;   ccdfile    - "opECalib.fit" 
;   hdr        - Header returned in memory
;   varfile    - Inverse Variance Frame after processing
;
;
pro sdssproc, infile, outfile, hdr=hdr, $
   configfile=configfile, ccdfile=ccdfile, varfile=varfile

   if (NOT keyword_set(configfile)) then configfile = 'opConfig.fit'
   if (NOT keyword_set(ccdfile)) then ccdfile = 'opECalib.fit'
   
   junk = findfile(configfile,count=ct)
   if (ct EQ 0) then $
    message, 'No configuration file ' + string(configfile)

   junk = findfile(ccdfile,count=ct)
   if (ct EQ 0) then $
    message, 'No ECalib file ' + string(ccdfile)

   rawdata = rdss_fits(infile, hdr)

   cards = sxpar(hdr,'NAXIS*')
;   if (cards[0] NE 2128 OR cards[1] NE 2069) then $
;      message, 'Expecting 2128x2069, found '+string(cards[0])+','$
;               +string(cards[1])

   names = str_sep(infile,'/')
   camrow = long( strmid( names[N_elements(names) -1], 4,1 ))
   camcol = long( strmid( names[N_elements(names) -1], 5,1 ))

   if (camrow NE 0) then message, 'not a spectroscopic image'
   if (camcol LT 1 OR camcol GT 4) then message, 'not a spectroscopic image'
 
;
;	Read in opConfig.fit
; 
   config = mrdfits(configfile, 1, confighdr)
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
   calib = mrdfits(ccdfile, 1, ccdhdr)
   calib = calib[ where(calib.camrow EQ camrow AND calib.camcol EQ camcol) ]

   gain = [calib.gain0,calib.gain1,calib.gain2,calib.gain3]
   readnoiseDN = [calib.readnoiseDN0,calib.readnoiseDN1,calib.readnoiseDN2, $
    calib.readnoiseDN3]

   ; Construct the final image
   igood = where(qexist)
   nr = max((srow+nrow)[igood])
   nc = max((scol+ncol)[igood])
   finalimg = fltarr(nc, nr)
   invvar = fltarr(nc, nr)

   for iamp=0, 3 do begin
      if (qexist[iamp] EQ 1) then begin
         if (nover[iamp] NE 0) then begin ; Use the "overscan" region
            biasval = median( $
             rawdata[sover[iamp]:sover[iamp]+nover[iamp]-1, $
             sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] )
         endif else if (nmapover[iamp] NE 0) then begin ; Use the "mapped overscan"
            biasval = median( $
             rawdata[smapover[iamp]:smapover[iamp]+nmapover[iamp]-1, $
             sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] )
         endif

         ; Copy the data for this amplifier into the final image
         finalimg[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
          rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
                  sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] - biasval

         invvar[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
           1.0/(abs(finalimg[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
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
    writefits, outfile, finalimg, hdr

   if (keyword_set(varfile)) then begin
      varhdr = hdr
      sxaddpar, varhdr, 'LOOKHERE', 'INVERSE VARIANCE of ' + outfile
      writefits, varfile, invvar, varhdr
   endif

   return
end
