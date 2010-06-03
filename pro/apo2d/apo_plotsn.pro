;+
; NAME:
;   apo_plotsn
;
; PURPOSE:
;   Generate S/N plot for one plate from a FITS logfile written by APOREDUCE.
;
; CALLING SEQUENCE:
;   apo_plotsn, logfile, plate, [ expnum=, plugdir=, plotfile= ]
;
; INPUTS:
;   logfile    - Logfile as written by APOREDUCE.  This is a FITS file
;                with an HDU of information for each reduced frame.
;   plate      - Plate number to plot.
;
; OPTIONAL KEYWORDS:
;   expnum     - If set, then make plot with S/N from exposures matching
;                this value (which can be an array of exposure numbers)
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
;   splog
;   readplugmap()
;
; REVISION HISTORY:
;   02-May-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
pro apo_plotsn, logfile, plate, expnum=expnum, plugdir=plugdir, $
 plotfile=plotfile

   if (NOT keyword_set(plate)) then return
   if (NOT keyword_set(plugdir)) then plugdir = './'

   platestr = string(plate, format='(i4.4)')
   splog, 'Generating S/N plot for plate '+platestr

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
   ; Read the plug map file for all fibers (both spectrographs)

   fullplugfile = filepath(plugfile, root_dir=plugdir)
   plugmap = readplugmap(fullplugfile, /deredden, /apotags, fibermask=fibermask)

   ;----------
   ; Loop through reductions for all science frames, and add S/N
   ; in quadrature, e.g. sum (S/N)^2.
   ; Only add (S/N)^2 that is not flagged as anything bad in the opLimits file

   for ii=0, n_elements(PPSCIENCE)-1 do begin
      meansn2 = PPSCIENCE[ii].sn2vector

      if (ii EQ 0) then begin
         nfiber = n_elements(meansn2) ; per spectrograph
         sn2array = fltarr(2, 2*nfiber)
      endif

      ; Test that the exposure falls within valid S/N^2 limits
      qkeep = apo_checklimits('science', 'SN2', PPSCIENCE[ii].camera, $
       PPSCIENCE[ii].sn2) EQ ''

      ; If EXPNUM is specified, then only use data from those exposure(s)
      if (keyword_set(expnum)) then $
       qkeep = qkeep AND (total(expnum EQ PPSCIENCE[ii].expnum) GT 0)

      if (qkeep) then begin
         case PPSCIENCE[ii].camera of
            'b1': sn2array[0,0:nfiber-1] += meansn2
            'b2': sn2array[0,nfiber:2*nfiber-1] += meansn2
            'r1': sn2array[1,0:nfiber-1] += meansn2
            'r2': sn2array[1,nfiber:2*nfiber-1] += meansn2
         endcase
      endif
   endfor

   ;----------
   ; Make the plot

   ; Lock the file to do this.
   if (keyword_set(plotfile)) then $
    while(djs_lockfile(plotfile, lun=plot_lun) EQ 0) do wait, 5

   plottitle = 'BOSS Spectro MJD=' + strtrim(string(mjd),2) $
    + ' Plate=' + strtrim(string(plate),2)
   if (keyword_set(expnum)) then $
    plottitle += ' exp=' + strtrim(expnum[0],2)
   plotsn, sqrt(sn2array), plugmap, sncode='sos', filter=['g','i'], $
    plottitle=plottitle, plotfile=plotfile

   if (keyword_set(plotfile)) then $
    djs_unlockfile, plotfile, lun=plot_lun
   splog, 'Plot finished'

   return
end
;------------------------------------------------------------------------------
