;+
; NAME:
;   spallflat
;
; PURPOSE:
;   Calling script for SPFLATTEN that generates pixel flats as specified in
;   a plan file.
;
; CALLING SEQUENCE:
;   spallflat, planfile=planfile
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name of output plan file; default to 'spPlanFlat.par'
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
;   spflatten
;   yanny_par()
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

pro spallflat, planfile=planfile

   if (NOT keyword_set(planfile)) then planfile = 'spPlanFlat.par'

   yanny_read, planfile, pdata, hdr=hdr

   ; Find the ONEEXP structure
   for i=0, N_elements(pdata)-1 do $
    if (tag_names(*pdata[i], /structure_name) EQ 'ONEEXP') then $
     allseq = *pdata[i]
   if (N_elements(allseq) EQ 0) then $
    message, 'No ONEEXP structures in plan file ' + planfile

   ; Find keywords from the header
   inputDir = yanny_par(hdr, 'inputDir')
   flatDir = yanny_par(hdr, 'flatDir')
   mjd = long(yanny_par(hdr, 'MJD'))

   camnames = ['b1', 'r2', 'b2', 'r1']
   camnums = ['01', '02', '03', '04']
   ncam = N_elements(camnames)

   for icam=0, ncam-1 do begin
      j = where(allseq[*].flavor EQ 'flat' $
       AND allseq[*].name[icam] NE 'UNKNOWN', ct)
      if (ct NE 0) then begin
         flatname = allseq[j].name[icam]
         outfile = 'pixflat-' + string(mjd,format='(i5.5)') $
          + '-' + camnums[icam] + '.fits'

      print, ''
      print, 'Generating pixel flat ' + flatDir+outfile
      spflatten, flatname, outfile=outfile, $
       indir=inputDir, outdir=flatDir, tmpdir=flatDir

      endif
   endfor

   return
end
;------------------------------------------------------------------------------
