;+
; NAME:
;   quickproc
;
; PURPOSE:
;   Read in Raw SDSS files, and keep as UINT, process with opConfig.par
;
; CALLING SEQUENCE:
;   image = quickproc(infile, hdr=hdr, configfile=configfile)
;    
;
; INPUTS:
;   infile     - Raw SDSS frame
;
; OPTIONAL KEYWORDS:
;   hdr        - Header returned in memory
;   configfile - Default to "opConfig.par"
;
; OUTPUTS:
;   image      - Processed 2d image (UINT) zero level 1000
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
;   10-Nov-1999  Written by Scott Burles, modified from sdssproc
;-
;------------------------------------------------------------------------------

function quickproc, infile, hdr=hdr, configfile=configfile

   if (N_params() LT 1) then begin
      print, 'Syntax - quickproc(infile, [hdr=hdr, configfile=configfile])'
      return, -1
   endif
   if (NOT keyword_set(configfile)) then configfile = 'opConfig.par'

   junk = findfile(configfile, count=ct)
   if (ct NE 1) then begin
     pp = getenv('EVIL_PAR') 
     tempname = findfile(filepath(configfile, root_dir=pp), count=ct)
   endif
   if (ct NE 1) then $
     message, 'No configuration file ' + string(configfile)

   realconfig = tempname[0]

   rawdata = rdss_fits(infile, hdr, /nofloat)

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
   camcol = fix( strmid(infile, i-2, 2) )

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
        sxaddpar, hdr, 'CAMERAS', 'b2       ', 'Guessed b2 by default'
;        message, 'Cannot determine CCD number from file name ' + infile
        end
   endcase
   camrow = 0

   sxaddpar, hdr, 'CAMROW', camrow
   sxaddpar, hdr, 'CAMCOL', camcol
 
   ; Read in opConfig.par file

   yanny_read, realconfig, pdata
   config = *pdata[0]
   ptr_free,pdata
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

   ; Construct the final image
   igood = where(qexist)
   nr = max((srow+nrow)[igood])
   nc = max((scol+ncol)[igood])
   image = uintarr(nc, nr)

   for iamp=0, 3 do begin
      if (qexist[iamp] EQ 1) then begin
         if (nover[iamp] NE 0) then begin
            ; Use the "overscan" region
            biasval = median( $
             rawdata[sover[iamp]:sover[iamp]+nover[iamp]-1, $
             sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] ) - 1000
         endif else if (nmapover[iamp] NE 0) then begin
            ; Use the "mapped overscan" region
            biasval = median( $
             rawdata[smapover[iamp]:smapover[iamp]+nmapover[iamp]-1, $
             sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] ) - 1000
         endif

         ; Copy the data for this amplifier into the final image
         image[scol[iamp]:scol[iamp]+ncol[iamp]-1, $
                  srow[iamp]:srow[iamp]+nrow[iamp]-1] = $
         rawdata[sdatacol[iamp]:sdatacol[iamp]+ncol[iamp]-1, $
                  sdatarow[iamp]:sdatarow[iamp]+nrow[iamp]-1] - biasval

         ; Add to the header
         sxaddpar, hdr, 'BIAS'+string(iamp,format='(i1)'), biasval
      endif
   endfor

   return, image
end
