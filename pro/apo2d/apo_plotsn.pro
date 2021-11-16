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
 plotfile=plotfile, fps=fps
  
   if (NOT keyword_set(plate)) then return
   if (NOT keyword_set(plugdir)) then plugdir = './'
   
   if (NOT keyword_set(fps)) then begin
      var_str='Plate'
   endif else begin
      var_str='Configuration'
   endelse

   platestr = plate_to_string(plate)
   splog, 'Generating S/N plot for plate '+platestr

   ;----------
   ; Read the science frames for this plate

   while(djs_lockfile(logfile) EQ 0) do wait, 5
   PPSCIENCE = mrdfits(logfile, 4)
   djs_unlockfile, logfile

   if (NOT keyword_set(PPSCIENCE)) then return
   ii = where(PPSCIENCE.config EQ plate)
   if (ii[0] EQ -1) then return
   PPSCIENCE = PPSCIENCE[ii]
   mjd = PPSCIENCE[0].mjd
   plugfile = PPSCIENCE[0].plugfile

   ;----------
   ; Read the plug map file for all fibers (both spectrographs)
   spd1=1
   fullplugfile = filepath(plugfile, root_dir=plugdir)
   
   if (Not keyword_set(fps)) then begin
      plugmap = readplugmap(fullplugfile,spd1,/deredden,/apotags, fibermask=fibermask,hdr=plhdr, /plates); included /deredden to match the SN2 in the html and plot-vivek
   endif else begin
      plugmap = readplugmap(fullplugfile, spd1, /deredden, /apotags, fibermask=fibermask, hdr=plhdr); included /deredden to match the SN2 in the html and plot-vivek
   endelse
   ;----------
   ; Loop through reductions for all science frames, and add S/N
   ; in quadrature, e.g. sum (S/N)^2.
   ; Only add (S/N)^2 that is not flagged as anything bad in the opLimits file

   for ii=0, n_elements(PPSCIENCE)-1 do begin
      meansn2 = PPSCIENCE[ii].sn2vector

      if (ii EQ 0) then begin
         nfiber = n_elements(meansn2) ; per spectrograph
         ;sn2array = fltarr(2, 2*nfiber)
         sn2array = fltarr(2, nfiber)
      endif

      ; Test that the exposure falls within valid S/N^2 limits
      qkeep = apo_checklimits('science', 'SN2', PPSCIENCE[ii].camera, $
       PPSCIENCE[ii].sn2) NE 'red'

      ; If EXPNUM is specified, then only use data from those exposure(s)
      if (keyword_set(expnum)) then $
       qkeep = qkeep AND (total(expnum EQ PPSCIENCE[ii].expnum) GT 0)

      if (qkeep) then begin
         case PPSCIENCE[ii].camera of
            'b1': sn2array[0,0:nfiber-1] += meansn2
            ;'b2': sn2array[0,nfiber:2*nfiber-1] += meansn2
            'r1': sn2array[1,0:nfiber-1] += meansn2
            ;'r2': sn2array[1,nfiber:2*nfiber-1] += meansn2
         endcase
      endif
   endfor
   ;----------
   ; Make the plot
   ; Treat the single exposure plot and total exposure plot differently -vivek
   ; In the blue cameras the SN plot for the plate, SN2 = n_exp * 4.0, -vivek
   ; just for consistency with the  html page -vivek  
   if (NOT keyword_set(expnum)) and n_elements(PPSCIENCE.expnum) gt 1 then begin
	nexp= n_elements(uniq(PPSCIENCE.expnum))
   endif else begin
	nexp=1
  endelse
   ; Lock the file to do this.
   if (keyword_set(plotfile)) then $
    while(djs_lockfile(plotfile, lun=plot_lun) EQ 0) do wait, 5

   plottitle = 'BOSS Spectro MJD=' + strtrim(string(mjd),2) $
    + ' '+var_str+'=' + strtrim(string(plate),2)
   if (keyword_set(expnum)) then $
    plottitle += ' exp=' + strtrim(expnum[0],2)
   ;; For ELG plates, call plotsn_elg rather than plotsn -vivek
   ;; plotsn_elg can handle n_exp ; filter =['g','z'] -vivek
   ;;if (plcnt gt 0) and (programname eq 'ELG_SGC' or programname eq 'ELG_NGC') then begin 
   ;;plotsn_elg, sqrt(sn2array), plugmap, sncode='sos', filter=['g','z'], $
   ;; plottitle=plottitle, plotfile=plotfile,snmin=0.2,nexp=nexp
   ;;endif else begin
       if (plate eq 7338 or plate eq 7339 or plate eq 7340) then begin
   ; Modifications for RM plates, plotsn_rm will have a scaled SN2 values
   ; RM plates to have higher depth than eBOSS plates -- vivek
   plotsn_rm, sqrt(sn2array), plugmap, sncode='sos', filter=['g','i'], $
    plottitle=plottitle, plotfile=plotfile,snmin=0.2
       endif else begin
    ;print, plugmap.mag
    ;print, sqrt(sn2array)
    plotsn, sqrt(sn2array), plugmap, sncode='sos', filter=['g','i'], $
    plottitle=plottitle, plotfile=plotfile,snmin=0.2
   endelse
   ;endelse

   if (keyword_set(plotfile)) then $
    djs_unlockfile, plotfile, lun=plot_lun
   splog, 'Plot finished'

   return
end
;------------------------------------------------------------------------------
