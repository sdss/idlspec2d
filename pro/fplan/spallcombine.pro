;+
; NAME:
;   spallcombine
;
; PURPOSE:
;   This is a Fermi-only routine.
;
; CALLING SEQUENCE:
;   spallcombine, mjd, [ topindir=, ncombine= ] 
;
; INPUTS:
;  mjd
;
; OPTIONAL INPUTS:
;   topindir   - Where should I start?  (default: '.')
;   ncombine   - How many exposures should be combined, the best N
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   spcombine
;   spplancomb
;   make2dmerge
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   24-Jan-2001  Written by Scott Burles.
;-
;------------------------------------------------------------------------------

pro spallcombine, mjd, topindir=topindir, ncombine=ncombine

   if NOT keyword_set(topindir) then topindir='.'
   if NOT keyword_set(ncombine) then ncombine=7

   cd, topindir, current=olddir

   spplancomb, topindir=topindir, /clobber, ncombine=ncombine
   planfile = findfile('spPlancomb*.par')
   if planfile[0] EQ '' then return
 
   nfile = n_elements(planfile)
   for i = 0, nfile - 1 do begin
     yanny_read, planfile[i], pdata, hdr=hdr
     hmm = where((*pdata[0]).mjd EQ mjd)
     yanny_free, pdata

     if hmm[0] NE -1 then begin
       spcombine, planfile[i]
       make2dmerge, planfile[i]  
     endif
   endfor

   cd, olddir
   return
end
;------------------------------------------------------------------------------
