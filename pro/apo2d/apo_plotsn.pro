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
;   apo_readlog
;   djs_lockfile()
;   djs_unlockfile
;   plotsn
;   sortplugmap()
;   yanny_read
;   yanny_free
;
; REVISION HISTORY:
;   02-May-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
pro apo_plotsn, logfile, plate, plugdir=plugdir, plotfile=plotfile

   if (NOT keyword_set(plate)) then return
   if (NOT keyword_set(plugdir)) then plugdir = './'

   platestr = string(plate, format='(i4.4)')
   bands = [1,3]

   ;----------
   ; Read the science frames for this plate

   pp = apo_readlog(logfile, plate=plate, flavor='science')
   if (NOT keyword_set(pp)) then return
   mjd = (*pp[0]).mjd
   plugfile = (*pp[0]).plugfile

   ;----------
   ; Read the plug map file for all 640 fibers

   fullplugfile = filepath(plugfile, root_dir=plugdir)
   yanny_read, fullplugfile, pdata
   plugmap = *pdata[0]
   yanny_free, pdata
   plugsort = sortplugmap(plugmap, fibermask=fibermask)

   ;----------
   ; Loop through reductions for all science frames, and add S/N
   ; in quadrature, e.g. sum (S/N)^2

   sn2array = fltarr(2, 640)
   for ii=0, n_elements(pp)-1 do begin
      meansn2 = (*pp[ii]).sn2vector
      case (*pp[ii]).camera of
         'b1': sn2array[0,0:319] =  sn2array[0,0:319] + meansn2
         'b2': sn2array[0,320:639] =  sn2array[0,320:639] + meansn2
         'r1': sn2array[1,0:319] =  sn2array[1,0:319] + meansn2
         'r2': sn2array[1,320:639] =  sn2array[1,320:639] + meansn2
      endcase
   endfor

   ;----------
   ; Make the plot

   ; Lock the file to do this.
   if (keyword_set(plotfile)) then $
    while(djs_lockfile(plotfile, lun=plot_lun) EQ 0) do wait, 5

   plottitle = 'APO SPECTRO MJD=' + strtrim(string(mjd),2) $
    + ' PLATE=' + strtrim(string(plate),2)
   plotsn, sqrt(sn2array), plugsort, bands=bands, plottitle=plottitle, $
    plotfile=plotfile

   if (keyword_set(plotfile)) then $
    djs_unlockfile, plotfile, lun=plot_lun

   return
end
