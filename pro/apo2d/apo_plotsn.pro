;+
; NAME:
;   apo_plotsn
;
; PURPOSE:
;   Generate S/N plot for one plate from a FITS logfile written by APOREDUCE.
;
; CALLING SEQUENCE:
;   apo_plotsn, logfile, plate, [ plugdir=, plotfile= ]
;
; INPUTS:
;   logfile    - Logfile as written by APOREDUCE.  This is a FITS file
;                with an HDU of information for each reduced frame.
;   plate      - Plate number to plot.
;
; OPTIONAL KEYWORDS:
;   plugdir    - Input directory for PLUGFILE; default to '.'
;                The name of the plugmap file is taken from the first
;                structure in the LOGFILE.
;   plotfile   - Name of plot file; if not specified then send plot to
;                current device.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   apo_checklimits()
;   djs_lockfile()
;   djs_unlockfile
;   mrdfits
;   plotsn
;   sortplugmap()
;   splog
;   readplugmap()
;
; REVISION HISTORY:
;   02-May-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
pro apo_plotsn, logfile, plate, plugdir=plugdir, plotfile=plotfile

   if (NOT keyword_set(plate)) then return
   if (NOT keyword_set(plugdir)) then plugdir = './'

   platestr = string(plate, format='(i4.4)')
   splog, 'Generating S/N plot for plate '+platestr
   bands = [1,3]

   ;----------
   ; Read the science frames for this plate

   while(djs_lockfile(logfile) EQ 0) do wait, 5
   PPSCIENCE = mrdfits(logfile, 4)
   djs_unlockfile, logfile

   if (NOT keyword_set(PPSCIENCE)) then return
   ii = where(PPSCIENCE.plate EQ plate)
   if (ii[0] EQ -1) then return
   PPSCIENCE = PPSCIENCE[ii]
   mjd = PPSCIENCE[0].mjd
   plugfile = PPSCIENCE[0].plugfile

   ;----------
   ; Read the plug map file for all 640 fibers

   fullplugfile = filepath(plugfile, root_dir=plugdir)
   plugmap = readplugmap(fullplugfile, /deredden, /apotags)
   plugsort = sortplugmap(plugmap, fibermask=fibermask)

   ;----------
   ; Loop through reductions for all science frames, and add S/N
   ; in quadrature, e.g. sum (S/N)^2.
   ; Only add (S/N)^2 that is not flagged as anything bad in
   ; the opLimits file (currently anything < 2.0 is bad).

   sn2array = fltarr(2, 640)
   for ii=0, n_elements(PPSCIENCE)-1 do begin
      meansn2 = PPSCIENCE[ii].sn2vector
      if (apo_checklimits('science', 'SN2', PPSCIENCE[ii].camera, $
                          PPSCIENCE[ii].sn2) EQ '') then begin
         case PPSCIENCE[ii].camera of
            'b1': sn2array[0,0:319] =  sn2array[0,0:319] + meansn2
            'b2': sn2array[0,320:639] =  sn2array[0,320:639] + meansn2
            'r1': sn2array[1,0:319] =  sn2array[1,0:319] + meansn2
            'r2': sn2array[1,320:639] =  sn2array[1,320:639] + meansn2
         endcase
      endif
   endfor

   ;----------
   ; Make the plot

   ; Lock the file to do this.
   if (keyword_set(plotfile)) then $
    while(djs_lockfile(plotfile, lun=plot_lun) EQ 0) do wait, 5

   plottitle = 'APO Spectro MJD=' + strtrim(string(mjd),2) $
    + ' Plate=' + strtrim(string(plate),2)
   plotsn, sqrt(sn2array), plugsort, bands=bands, snmag=[0,20.33,0,20.06,0], $
    plottitle=plottitle, plotfile=plotfile

   if (keyword_set(plotfile)) then $
    djs_unlockfile, plotfile, lun=plot_lun
   splog, 'Plot finished'

   return
end
;------------------------------------------------------------------------------
