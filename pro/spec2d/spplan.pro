;begin+
; NAME:
;   spplan
;
; PURPOSE:
;   Create a plan file for running the Spectro-2D pipeline.
;
; CALLING SEQUENCE:
;   spplan, [ indir=indir, plugdir=plugdir, flatdir=flatdir, $
;    mjd=mjd, planfile=planfile, run=run, /flats, /checkstats ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   indir      - Input directory to look for files; default to '.'
;   plugdir    - Directory for plug map files; default to same as INDIR
;   flatdir    - Directory for pixel flat files; default to same as INDIR
;                unless there are wildcards -- in that case, default to '.'
;   mjd        - Modified Julian date; default to the MJD listed in the
;                header of the first FITS file
;   planfile   - Name of output plan file; default to 'spPlan2d.par'
;                unless FLATS is set, then default to 'spPlanFlat.par'
;   flats      - Set this keyword to plan the spPlanFlat file
;   checkstats - Set this keyword to read in full frame and do statistics
;	         on the frames to verify flavor
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;   Create the plan file "spPlanFlat.par" for making pixel flats by looking
;   for all flats in all subdirectories beginning with "514*" :
;     spplan, indir='514*/', /flats
;
;   Create the plan file "spPlan2d.par" for reducing one night of data
;   in the directory "51441/" :
;     spplan, indir='51441/'
;
; BUGS:
;   Need to implement lampfile.
;
; PROCEDURES CALLED:
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;   spplan_create_exp
;   spplan_create_plug
;   spplan_create_pixflats
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;   10-Nov-1999    SMB added checkstats to verify flavors
;-
;------------------------------------------------------------------------------

function spplan_create_exp, seqid, plateid, flavor
   badname = 'UNKNOWN            ' ; Make same length as other file names
   oneexp = {oneexp, seqid: 0L, plateid: 0L, flavor: '', $
    name: [badname, badname, badname, badname] }
   oneexp.seqid = seqid
   oneexp.plateid = plateid
   oneexp.flavor = flavor
   return, oneexp
end

;------------------------------------------------------------------------------

function spplan_create_plug, seqid, plateid
   badname = 'UNKNOWN'
   oneplug = {oneplug, seqid: 0L, plateid: 0L, name: badname}
   oneplug.seqid = seqid
   oneplug.plateid = plateid
   return, oneplug
end

;------------------------------------------------------------------------------

function spplan_create_pixflats
   badname = 'UNKNOWN'
   pixflats = {pixflats, name: [badname, badname, badname, badname] }
   return, pixflats
end

;------------------------------------------------------------------------------

pro spplan, indir=indir, plugdir=plugdir, flatdir=flatdir, $
 mjd=mjd, planfile=planfile, flats=flats, checkstats=checkstats, $
 run=run, root2d=root2d, combroot = combroot

   if (NOT keyword_set(run)) then run = 0
   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(root2d)) then begin 
      if (strpos(indir, '*') EQ -1) then root2d = indir $
       else root2d = '.'
   endif 

   if (NOT keyword_set(combroot)) then begin
      if (strpos(indir, '*') EQ -1) then combroot = indir $
       else combroot = '.'
   endif 


   if (NOT keyword_set(flatdir)) then begin
      if (strpos(indir, '*') EQ -1) then flatdir = indir $
       else flatdir = '.'
   endif 

   if (NOT keyword_set(plugdir)) then begin
      if (strpos(indir, '*') EQ -1) then plugdir = indir $
       else plugdir = '.'
   endif 

   if (NOT keyword_set(planfile)) then begin
      if (keyword_set(flats)) then planfile = 'spPlanFlat.par' $
       else planfile = 'spPlan2d.par'
   endif

   camnames = ['b1', 'r2', 'b2', 'r1']
   camnums = ['01', '02', '03', '04']
   ncam = N_elements(camnames)

   fullname = findfile(filepath('*.fit',root_dir=indir), count=nfile)

   if (nfile EQ 0) then $
    message, 'No files found.'

   ; Remove the path from the file names
   shortname = strarr(nfile)
   for i=0, nfile-1 do begin
      res = str_sep(fullname[i],'/')
      shortname[i] = res[N_elements(res)-1]
   endfor

   ; Sort the files based upon exposure number + camera number
   isort = sort( strmid(shortname,7,8) + strmid(shortname,4,2) )
   fullname = fullname[isort]
   shortname = shortname[isort]

   ; Find all useful header keywords
   PLATEID = lonarr(nfile)
   EXPTIME = fltarr(nfile)
   EXPOSURE = fltarr(nfile)
   FLAVOR = strarr(nfile)
   CAMERAS = strarr(nfile)
   for i=0, nfile-1 do begin
      if (NOT keyword_set(checkstats)) then hdr = headfits(fullname[i]) $
      else image = quickproc(fullname[i],hdr=hdr)


      PLATEID[i] = long( sxpar(hdr, 'PLATEID') )
      EXPTIME[i] = sxpar(hdr, 'EXPTIME')
      EXPOSURE[i] = long( sxpar(hdr, 'EXPOSURE') )
      FLAVOR[i] = strtrim(sxpar(hdr, 'FLAVOR'),2)
      CAMERAS[i] = strmid(shortname[i],4,2) ; Camera number from file name


      goodcamera = where(CAMERAS[i] EQ camnames,ct)
      if (ct EQ 1) then CAMERAS[i] = camnums[goodcamera] $
      else begin
        goodcamera = where(CAMERAS[i] EQ camnums OR CAMERAS[i] EQ camnames,ct)
        if (ct NE 1) then message, 'Camera number is not in file name'
      endelse

      if (keyword_set(checkstats)) then $
        newflavor = checkflavor(image, flavor[i], camnames[goodcamera])
	
      ; Rename 'target' as 'science'
      if (FLAVOR[i] EQ 'target') then FLAVOR[i] = 'science'

      ; Read the MJD from the header of the first file
      if (NOT keyword_set(mjd) AND i EQ 0) then mjd = sxpar(hdr, 'MJD')
   endfor

   ; Find the start of each new sequence not by looking at the recorded
   ; sequence number, which may be wrong, but by looking at where the
   ; PLATEID has changed.
   begseq = 0L
   lastplate = PLATEID[0]
   for i=1, nfile-1 do begin
      if (PLATEID[i] NE lastplate) then begin
         begseq = [begseq, i]
         lastplate = PLATEID[i]
      endif
   endfor

   ; Loop through each plate sequence
   allseq = 0
   nseq = N_elements(begseq)
   for iseq=0, nseq-1 do begin

      oneseq = 0

      ; Determine all files that are in this sequence
      if (iseq LT nseq-1) then $
       ifile = lindgen(begseq[iseq+1]-begseq[iseq]) + begseq[iseq] $
      else $
       ifile = lindgen(nfile-begseq[iseq]) + begseq[iseq]

      ; Set the sequence ID equal to the first exposure number
      ; Also, get the plate ID from that first frame
      seqid = EXPOSURE[ifile[0]]
      pltid = PLATEID[ifile[0]]

      qdone = bytarr(N_elements(ifile)) ; Set equal to 1 as each seq frame is
                                        ; written to a structure

      ; Only look at those frames labelled as 'flat', 'arc', or 'science'
      ignore = where(FLAVOR[ifile] NE 'flat' $
                 AND FLAVOR[ifile] NE 'arc' $
                 AND FLAVOR[ifile] NE 'science')
      if (ignore[0] NE -1) then qdone[ignore] = 1

      while (min(qdone) EQ 0) do begin
         inotdone = where(qdone EQ 0)
         indx = ifile[inotdone]
         indx = indx[ where( EXPOSURE[indx] EQ EXPOSURE[indx[0]] $
                         AND FLAVOR[indx] EQ FLAVOR[indx[0]] ) ]
         oneexp = spplan_create_exp(seqid, pltid, FLAVOR[indx[0]])

         for icam=0, ncam-1 do begin
            j = where(CAMERAS[indx] EQ camnums[icam], ct)
            if (ct EQ 1) then begin
               oneexp.name[icam] = shortname[indx[j[0]]]
               qdone[inotdone[j[0]]] = 1
            endif else if (ct GT 1) then begin
               message, 'Several frames with EXPOSURE=' $
                + string(EXPOSURE[indx[0]]) $
                + ' and CAMERA=' + camnums[icam]
            endif
         endfor

         ; Add all of these exposures to the ONESEQ structure
         if (NOT keyword_set(flats) OR $
             (keyword_set(flats) AND oneexp.flavor EQ 'flat') ) then begin
            if (keyword_set(oneseq)) then oneseq = [oneseq, oneexp] $
             else oneseq = oneexp
         endif

      endwhile

      if (NOT keyword_set(flats) AND keyword_set(oneseq)) then begin
         ; Test that a flat and an arc exist for any camera with a
         ; science exposure
         for icam=0, ncam-1 do begin
            junk = where(oneseq.flavor EQ 'science' $
                     AND oneseq.name[icam] NE 'UNKNOWN', nscience)
            junk = where(oneseq.flavor EQ 'flat' $
                     AND oneseq.name[icam] NE 'UNKNOWN', nflat)
            junk = where(oneseq.flavor EQ 'arc' $
                     AND oneseq.name[icam] NE 'UNKNOWN', narc)
            if (nscience GT 0 AND nflat EQ 0) then $
             print, 'No flat for SEQID=' + string(oneseq[0].seqid) $
             + ' and PLATEID=' + string(oneseq[0].plateid)
            if (nscience GT 0 AND narc EQ 0) then $
             print, 'No arc for SEQID= ' + strtrim(string(oneseq[0].seqid),2) $
             + ', PLATEID= ' + strtrim(string(oneseq[0].plateid),2) $
             + ', CAMERA= ' + strtrim(string(camnums[icam]),2)
         endfor

         ; Find a plug map file
         oneplug = spplan_create_plug(seqid, pltid)
         if (pltid GT 0 AND pltid LT 9999) then $
          platestr = string(pltid,format='(i04.4)') $
          else platestr = '0000'
         files = 'plPlugMapM-' + platestr + '*.par'
         files = findfile(filepath(files,root_dir=plugdir), count=ct)
         if (ct GT 1) then begin
            print, 'Several plug map files found for plate number ' $
             + string(platenum)
            print, 'Using file ', files[0]
            files = files[0]
            ct = 1
         endif
         if (ct EQ 1) then begin
            ; Take the only plugmap file
            j = rstrpos(files[0], '/')
            if (j EQ -1) then oneplug.name = files[0] $
             else oneplug.name = strmid(files[0],j+1,strlen(files[0])-j)

         endif
         if (ct EQ 0) then begin
            print, 'No plug map files found for plate number ' + string(pltid)
         endif

         ; Add this plug map file to the ALLPLUG structure
         if (N_elements(allplug) EQ 0) then allplug = oneplug $
          else allplug = [allplug, oneplug]

      endif

      ; Add all exposures in ONESEQ to the ALLSEQ structure
      if (keyword_set(allseq) AND keyword_set(oneseq)) then $
       allseq = [allseq, oneseq] $
      else allseq = oneseq

   endfor

   if (NOT keyword_set(flats)) then begin
      ; Look for the pixel flats
      pixflats = spplan_create_pixflats()
      for icam=0, ncam-1 do begin
         files = 'pixflat-*-' + camnums[icam] + '.fits'
         files = findfile(filepath(files,root_dir=flatdir), count=ct)
         if (ct EQ 0) then begin
           files = 'pixflat-*-' + camnames[icam] + '.fits'
           files = findfile(filepath(files,root_dir=flatdir), count=ct)
         endif
         if (ct GT 1) then begin
            print, 'Several pixel flats found for CAMERA= ' + camnums[icam]
            print, 'Using file ', files[0]
            files = files[0]
            ct = 1
         endif
         if (ct EQ 1) then begin
            ; Take the only pixel flat
            j = rstrpos(files[0], '/')
            if (j EQ -1) then pixflats.name[icam] = files[0] $
             else pixflats.name[icam] = strmid(files[0],j+1,strlen(files[0])-j)
         endif
         if (ct EQ 0) then begin
            print, 'No pixel flats found for CAMERA= ' + camnums[icam]
         endif
      endfor

   endif

   ; Create keyword pairs for plan file
   hdr = ''


   if (NOT keyword_set(flats)) then begin
   hdr = [hdr, "configDir   'XXX'     # Directory for CCD calibration files"]
   hdr = [hdr, "ccdConfig   'XXX'     # File name for amplifier configuration"]
   hdr = [hdr, "ccdECalib   'XXX'     # File name for electronic calibrations"]
   hdr = [hdr, "ccdBC       'XXX'     # File name for bad pixels"]
   hdr = [hdr, "plugDir     '" + plugdir + "'  # Directory for plugmap files"]
   hdr = [hdr, "extractDir  '" + root2d + "'  # Root dir for 2d spectra"]
   hdr = [hdr, "combDir     '" + combroot + $
                                          "'  # Root dir for combined spectra"]
   endif

   hdr = [hdr, "inputDir    '" + indir + "'  # Directory for raw images"]
   hdr = [hdr, "flatDir     '" + flatdir + "'     # Directory for pixel flats"]
   hdr = [hdr, "MJD     " + string(mjd) + "  # Modified Julian Date"]
   hdr = [hdr, "run         " + string(run) + "  # Modified Julian Date"]

   if (NOT keyword_set(flats)) then $
    yanny_write, planfile, [ptr_new(pixflats), ptr_new(allplug), ptr_new(allseq)], hdr=hdr $
    else $
    yanny_write, planfile, [ptr_new(allseq)], hdr=hdr

   return
end
;------------------------------------------------------------------------------
