;+
; NAME:
;   spallreduce
;
; PURPOSE:
;   Calling script for SPREDUCE and COMBINE2DOUT that reduces a night
;   of data according to a plan file.
;
; CALLING SEQUENCE:
;   spallreduce, planfile=planfile
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name of output plan file; default to 'spPlan2d.par'
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Rather than use FINDFILE, we should know what all the idlout* files are.
;   Need to implement lampfile.
;   Pass QA flags to SPREDUCE
;
; PROCEDURES CALLED:
;   combine2dout
;   spreduce
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

pro spallreduce, planfile=planfile

   if (NOT keyword_set(planfile)) then planfile = 'spPlan2d.par'

   yanny_read, planfile, pdata, hdr=hdr

   ; Find the ONEEXP structure
   for i=0, N_elements(pdata)-1 do begin
      if (tag_names(*pdata[i], /structure_name) EQ 'PIXFLATS') then $
       pixflats = *pdata[i]
      if (tag_names(*pdata[i], /structure_name) EQ 'ONEPLUG') then $
       allplug = *pdata[i]
      if (tag_names(*pdata[i], /structure_name) EQ 'ONEEXP') then $
       allseq = *pdata[i]
   endfor
   if (N_elements(allseq) EQ 0) then $
    message, 'No ONEEXP structures in plan file ' + planfile

   ; Find keywords from the header
   plugDir = yanny_par(hdr, 'plugDir')
   extractDir = yanny_par(hdr, 'extractDir')
   inputDir = yanny_par(hdr, 'inputDir')
   flatDir = yanny_par(hdr, 'flatDir')
   mjd = yanny_par(hdr, 'MJD')

   camnames = ['b1', 'r2', 'b2', 'r1']
   camnums = ['01', '02', '03', '04']
   ncam = N_elements(camnames)

   ; Find all the sequence IDs
   seqid = allseq[ sort(allseq.seqid) ].seqid
   seqid = allseq[ uniq(allseq.seqid) ].seqid

   for iseq=0, N_elements(seqid)-1 do begin

      ; Get the plate ID number from any (e.g., the first) exposure with
      ; this sequence ID number
      j = where(allseq.seqid EQ seqid[iseq])
      plateid = allseq[j[0]].plateid

      ; Find the corresponding plug map file
      j = where(allplug.seqid EQ seqid[iseq] $
            AND allplug.plateid EQ plateid )
      if (j[0] NE -1) then plugfile = allplug[j[0]].name $
       else message, 'No plug map file for SEQID= ' $
        + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
        + strtrim(string(plateid),2)

      for icam=0, ncam-1 do begin

         ; Find the corresponding pixel flat
         pixflatname = pixflats.name[icam]
         if (pixflatname EQ 'UNKNOWN') then $
          message, 'No pixel flat for CAMERA= ' + camnums[icam]

         j = where(allseq.seqid EQ seqid[iseq] $
               AND allseq.flavor EQ 'science' $
               AND allseq.name[icam] NE 'UNKNOWN' )

         if (j[0] NE -1) then begin

            ; String array with all science exposures at this sequence + camera
            objname = allseq[j].name[icam]

            ; Select the first flat exposure at this sequence + camera
            j = where(allseq.seqid EQ seqid[iseq] $
                  AND allseq.flavor EQ 'flat' $
                  AND allseq.name[icam] NE 'UNKNOWN' )
            if (j[0] NE -1) then flatname = allseq[j[0]].name[icam] $
             else message, 'No flat for SEQID= ' $
              + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
              + strtrim(string(plateid),2) + ', CAMERA= ' + camnums[icam]

            ; Select the first arc exposure at this sequence + camera
            j = where(allseq.seqid EQ seqid[iseq] $
                  AND allseq.flavor EQ 'arc' $
                  AND allseq.name[icam] NE 'UNKNOWN' )
            if (j[0] NE -1) then arcname = allseq[j[0]].name[icam] $
             else message, 'No arc for SEQID= ' $
              + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
              + strtrim(string(plateid),2) + ', CAMERA= ' + camnums[icam]

            spreduce, flatname, arcname, objname, $
             pixflatname=flatDir+pixflatname, $
             plugfile=plugfile, lampfile=lampfile, $
             indir=inputDir, plugdir=plugDir, outdir=extractDir, $
             qadir=extractDir

         endif

      endfor ; End loop for camera number

      ; Combine all red+blue exposures for a given sequence

      for side = 1, 2 do begin
         for i=1, 320 do begin

            outputfile = 'idlout-'+string(format='(i1,a,i4.4,a,i3.3,a)',side, $
             '-',plateid,'-',i,'.fit')

            expres = string(format='(a,i1,a,i3.3,a)', 's-', side, '*', i,'.fit')
            files = findfile(extractDir+expres)

            combine2dout, files, extractDir+outputfile, wavemin = alog10(3750.0)
         endfor
      endfor

   endfor ; End loop for sequence number

   return
end
;------------------------------------------------------------------------------
