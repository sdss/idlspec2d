;+
; NAME:
;   spallreduce
;
; PURPOSE:
;   Calling script for SPREDUCE and COMBINE2DOUT that reduces a night
;   of data according to a plan file.
;
; CALLING SEQUENCE:
;   spallreduce, [ planfile, docams=, /nocombine, /xdisplay ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name(s) of output plan file; default to reducing all
;                plan files matching 'spPlan2d*.par'
;   docams     - Cameras to reduce; default to ['b1', 'b2', 'r1', 'r2'];
;                set to 0 or '' to disable running SPREDUCE and to only
;                combine files.
;   nocombine  - Only run SPREDUCE, not COMBINE2DOUT.
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

pro spallreduce, planfile, docams=docams, nocombine=nocombine, $
 xdisplay=xdisplay

   if (NOT keyword_set(planfile)) then planfile = findfile('spPlan2d*.par')
   if (N_elements(planfile) NE 1) then $
      message, 'Please just give me one plan file'

   ;-------------------------------------------------------------
   ;  First just do standard 2d reduction
   ;
   spreduce2d, planfile[0], docams=docams, xdisplay=xdisplay

        
   ;-------------------------------------------------------------
   ;  Now coadd, and also produce stopgap 2dmerge directory
   ;

   if NOT keyword_set(nocombine) then begin
     spplancomb, topindir='.'
     planfile = findfile('spPlancomb*.par')
     spcombine, planfile, adderr=0
   endif

   return
end
;------------------------------------------------------------------------------
