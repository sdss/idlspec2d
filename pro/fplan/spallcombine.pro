;+
; NAME:
;   spallcombine
;
; PURPOSE:
;   This is a Fermi-only routine.
;
; CALLING SEQUENCE:
;   spallcombine, mjd, [ topindir ]
;
; INPUTS:
;  mjd
;
; OPTIONAL INPUTS:
;   topindir   - Where should I start?  (default: '.')
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

   cd, topindir, current=olddir

   spplancomb, topindir=topindir, mjd=mjd, /clobber, ncombine=ncombine
   planfile = findfile('spPlancomb*.par')
   if planfile[0] EQ '' then return
 
   nfile = n_elements(planfile)
   mjdlist = lonarr(nfile)
   
   for i = 0, nfile - 1 do begin
     yanny_read, planfile[i], pdata, hdr=hdr
     if (size(pdata, /tname) NE 'INT') then begin
       hmm = where((*pdata[0]).mjd EQ mjd)
       if hmm[0] NE -1 then mjdlist[i] = ((*pdata[0]).mjd)[hmm[0]]
       yanny_free, pdata
     endif
   endfor

   thelatest = max(mjdlist, place)

   spcombine, planfile[place]
   make2dmerge, planfile[place]  

   cd, olddir
   return
end
;------------------------------------------------------------------------------
