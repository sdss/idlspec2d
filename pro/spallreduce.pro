;+
; NAME:
;   spallreduce
;
; PURPOSE:
;   Calling script for SPREDUCE and COMBINE2DOUT that reduces a night
;   of data according to a plan file.
;
; CALLING SEQUENCE:
;   spallreduce, planfile=planfile, combineonly=combineonly, docams=docams
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

pro spallreduce, planfile=planfile, combineonly=combineonly, docams=docams

   if (NOT keyword_set(planfile)) then planfile = 'spPlan2d.par'

   docomb = 0
   if (NOT keyword_set(docams)) then begin
       docams = ['b1', 'r2', 'b2', 'r1'] ; do all cameras
       docomb = 1
   endif
   ndo = N_elements(docams)

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

   ptr_free, pdata

   if (N_elements(allseq) EQ 0) then $
    message, 'No ONEEXP structures in plan file ' + planfile

   ; Find keywords from the header
   plugDir = yanny_par(hdr, 'plugDir')
   extractDir = yanny_par(hdr, 'extractDir')
   combDir = yanny_par(hdr, 'combDir')
   if (combDir EQ '') then combDir=extractDir

   inputDir = yanny_par(hdr, 'inputDir')
   flatDir = yanny_par(hdr, 'flatDir')
   mjd = yanny_par(hdr, 'MJD')
   run = strtrim(yanny_par(hdr, 'run'),2)
   if (run EQ '') then run = '0'

  
   camnames = ['b1', 'r2', 'b2', 'r1']
   camnums = ['01', '02', '03', '04']
   ncam = N_elements(camnames)

   ; Find all the sequence IDs
   seqid = allseq[ sort(allseq.seqid) ].seqid
   seqid = allseq[ uniq(allseq.seqid) ].seqid

   set_plot,'ps'
   for iseq=0, N_elements(seqid)-1 do begin

      ; Get the plate ID number from any (e.g., the first) exposure with
      ; this sequence ID number
      j = where(allseq.seqid EQ seqid[iseq])
      plateid = allseq[j[0]].plateid

      plateDir=filepath(strtrim(string(plateid),2)+'/2d_'+ $
           run,root_dir=extractDir)
      combineDir=filepath(strtrim(string(plateid),2)+'/comb_'+ $
           run,root_dir=combDir)

    if (NOT keyword_set(combineonly)) then begin

      ; Find the corresponding plug map file
      j = where(allplug.seqid EQ seqid[iseq] $
            AND allplug.plateid EQ plateid )
      if (j[0] NE -1) then plugfile = allplug[j[0]].name $
       else message, 'No plug map file for SEQID= ' $
        + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
        + strtrim(string(plateid),2)

      for ido=0, ndo-1 do begin

         icam = where(camnames EQ docams[ido], camct)
         if (camct NE 1) then message, 'Non-unique camera id'

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

;
;	Changed flats and arcs to take the list of names in extract_image
;
            if (j[0] NE -1) then tempflatspot = j $
             else message, 'No flat for SEQID= ' $
              + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
              + strtrim(string(plateid),2) + ', CAMERA= ' + camnums[icam]

            ; Select the first arc exposure at this sequence + camera
            j = where(allseq.seqid EQ seqid[iseq] $
                  AND allseq.flavor EQ 'arc' $
                  AND allseq.name[icam] NE 'UNKNOWN' , numarcs)
            if (j[0] NE -1) then temparcspot = j $
             else message, 'No arc for SEQID= ' $
              + strtrim(string(seqid[iseq]),2) + ', PLATEID= ' $
              + strtrim(string(plateid),2) + ', CAMERA= ' + camnums[icam]


	    spawn, 'mkdir -p '+plateDir	

;
;	Need to make sure we have a one-to-some match between
;	arcs and flats 
;	  step through arcs and pick flat with closest seqid number

          newflatspot = temparcspot
	  for i=0,numarcs-1 do begin
	     mindist = min(abs(tempflatspot-temparcspot[i]),closestspot)
             newflatspot[i] = tempflatspot[closestspot]
	  endfor

          ; Use reverse to check last arc first
	  arcname = allseq[reverse(temparcspot)].name[icam]
	  flatname = allseq[reverse(newflatspot)].name[icam]

            spreduce, flatname, arcname, objname, $
             pixflatname=filepath(pixflatname,root_dir=flatDir), $
             plugfile=plugfile, lampfile=lampfile, $
             indir=inputDir, plugdir=plugDir, $
             outdir=plateDir, qadir=plateDir

          heap_gc   ; garbage collection
         endif

      endfor ; End loop for camera number

    endif
      ; Combine all red+blue exposures for a given sequence

      if (docomb) then begin

	spawn, 'mkdir -p '+combineDir	
	startcombtime = systime(1)

        for side = 1, 2 do begin
          for i=1, 320 do begin

            outputfile = 'idlout-'+string(format='(i1,a,i4.4,a,i3.3,a)',side, $
             '-',plateid,'-',i,'.fit')

            expres = string(format='(a,i1,a,i3.3,a)', 's-', side, '*', i,'.fit')
            files = findfile(filepath(expres,root_dir=plateDir))

	    if (files[0] EQ '') then $
               print, 'No files found for ', i, ' side ', side $
            else $
              combine2dout, files, filepath(outputfile,root_dir=combineDir), $
               wavemin = alog10(3750.0)
         endfor
       endfor
       print, 'Finished combining sequence', seqid[iseq], ' in', $
             systime(1)-startcombtime, ' seconds'
     endif

   endfor ; End loop for sequence number

   device, /close
   set_plot, 'x'
   return
end
;------------------------------------------------------------------------------
