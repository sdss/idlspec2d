;+
; NAME:
;   spallreduce
;
; PURPOSE:
;   This is a Fermi-only routine.
;   Calling script for SPREDUCE and COMBINE2DOUT that reduces a night
;   of data according to a plan file.
;
; CALLING SEQUENCE:
;   spallreduce, [ planfile, docams=, /ncombine, /xdisplay ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name(s) of output plan file; default to reducing all
;                plan files matching 'spPlan2d*.par'
;   docams     - Cameras to reduce; default to ['b1', 'b2', 'r1', 'r2'];
;                set to 0 or '' to disable running SPREDUCE and to only
;                combine files.
;   ncombine   - How many exposures to combine, the best N
;   xdisplay   - Send plots to X display rather than to plot file
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Try to avoid using SPAWN.
;   How do I just combine files? and not re-run spreduce?
;
;   What happens when the gain has changed?
;   We need to list opECalib.par in spPlan file.
;   Use default opECalib if not found...
;
; PROCEDURES CALLED:
;   spreduce2d
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;   26-Sep-2000  Reinvoked after the incredible screw up.
;-
;------------------------------------------------------------------------------

pro spallreduce, planfile, docams=docams, ncombine=ncombine, $
 xdisplay=xdisplay

   if (NOT keyword_set(planfile)) then planfile = findfile('spPlan2d*.par')
   if n_elements(ncombine) EQ 0 then ncombine=7

   if (N_elements(planfile) NE 1) then $
      message, 'Please just give me one plan file'

   ;-------------------------------------------------------------
   ;  First just do standard 2d reduction

   ; Read path names from plan file
   yanny_read, planfile[0], hdr=hdr
   inputdir = yanny_par(hdr, 'inputdir')
   plugdir = yanny_par(hdr, 'plugdir')
   flatdir = yanny_par(hdr, 'flatdir')
   mjd     = yanny_par(hdr, 'MJD')

   ; Remove the MJD from INPUTDIR, PLUGDIR
   junk = fileandpath(inputdir, path=inputdir)
   junk = fileandpath(plugdir, path=plugdir)

   ; Set environment variables for call to SPREDUCE2D
   setenv, 'RAWDATA_DIR=' + inputdir
   setenv, 'SPECLOG_DIR=' + plugdir
   setenv, 'SPECFLAT_DIR=' + flatdir

   spreduce2d, planfile[0], docams=docams, xdisplay=xdisplay

   ;-------------------------------------------------------------
   ;  Now coadd, and also produce stopgap 2dmerge directory
   ;

   if NOT keyword_set(ncombine) then return
   spallcombine, mjd, topindir='.', ncombine=ncombine

   return
end
;------------------------------------------------------------------------------
