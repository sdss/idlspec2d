;+
; NAME:
;   spplan
;
; PURPOSE:
;   Create plan file(s) for running the Spectro-2D pipeline.
;
; CALLING SEQUENCE:
;   spplan, [ rawdir, astrolog=, mjd=, flatdir=, minexp= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   rawdir     - Search for raw data files in RAWDIR/MJD/*.
;                This should be an absolute file path, and we default to
;                '/usr/sdss/data05/spectro/rawdata'.
;   astrolog   - Search for plug-map files in PLUGDIR/MJD/*.
;                This should be an absolute file path, and we default to
;                '/usr/sdss/data05/spectro/astrolog'.
;   mjd        - Look for raw data files in RAWDIR/MJD; default to '*' to
;                search all subdirectories.  Note that this need not be
;                integer-valued, but could be for example '51441_test'.
;   flatdir    - Directory for pixel flat files.  For now, default
;                to 'pixflat'.
;   minexp     - Minimum exposure time for science frames; default to 300 sec.
;
; OUTPUT:
;
; COMMENTS:
;   Look for the input files in:
;     RAWDIR/MJD/sdR-cs-eeeeeeee.fit          - Raw frames, c=color, s=spec #
;     RAWDIR/MJD/plPlugMapM-pppp-mmmmm-aa.par - Plug map files, aa=mapper
;
;   The top-level of the output directory structure, TOPDIR, is 2d_VERSION,
;   where we read the version with IDLSPEC2d_VERSION().  If this is not a CVS-
;   declared version of the code, then the output directory is '2d_test'.
;   The files created are:
;     TOPDIR/PLATE/spPlan2d-mmmmm-pppp.par    - 2D plan file (could be several)
;
;   If an output plan file already exists, this procedure will crash rather
;   than overwrite that file.
;
;   Note we assume Unix directory names and we SPAWN the Unix 'ls' command.
;
; EXAMPLES:
;   Create the plan file(s) for reducing the data for MJD=51441, with that
;   raw data in the directory '/u/schlegel/rawdata/51441'
;     spplan, '/u/schlegel/rawdata', mjd=51441
;
; BUGS:
;
; PROCEDURES CALLED:
;   headfits()
;   idlspec2d_version()
;   splog
;   sxpar()
;   yanny_read
;   yanny_write
;
; INTERNAL SUPPORT ROUTINES:
;   spplan_create_exp
;   spplan_create_plug
;   spplan_create_pixflats
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;   10-Nov-1999  SMB added checkstats to verify flavors.
;   18-Jan-2000  Re-written to conform to new file names + directory structure.
;                Not using CHECKSTATS any more at the moment.
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

pro spplan, rawdir, astrolog=astrolog, mjd=mjd, flatdir=flatdir, minexp=minexp

   if (NOT keyword_set(rawdir)) then begin
      if ((findfile('/usr/sdss/data05/spectro/rawdata'))[0] NE '') then $
       rawdir = '/usr/sdss/data05/spectro/rawdata' $
      else if ((findfile('/home/schlegel/data/rawdata'))[0] NE '') then $
       rawdir = '/home/schlegel/data/rawdata' $
      else $
       message, 'Must specify RAWDIR'
   endif
   if (NOT keyword_set(astrolog)) then $
    astrolog = strmid(rawdir, 0, rstrpos(rawdir,'/')+1) + 'astrolog'
   if (NOT keyword_set(flatdir)) then flatdir = 'pixflat'
   if (NOT keyword_set(minexp)) then minexp = 300

   ;----------
   ; Determine the top-level of the output directory tree, and quit if
   ; it already exists and contains some files.

   vers = idlspec2d_version()
   if (strpos(vers, 'NOCVS') EQ -1) then topdir = '2d_' + vers $
    else topdir = '2d_test'
   splog, 'Setting top-level of output directory to ' + topdir

   ; Create a list of the MJD directories (as strings)
   if (NOT keyword_set(mjd)) then begin
      spawn, '\ls -d '+rawdir+'/*', mjdlist
   endif else begin
      tmpstring = ''
      for i=0, N_elements(mjd)-1 do $
       tmpstring = tmpstring + ' ' + rawdir+'/'+strtrim(string(mjd[i]),2)
      spawn, '\ls -d '+tmpstring, mjdlist
   endelse

   ; Strip leading directory names from MJDLIST
   for imjd=0, N_elements(mjdlist)-1 do begin
      i = rstrpos(mjdlist[imjd], '/') > 0
      mjdlist[imjd] = strmid(mjdlist[imjd], i)
   endfor

   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   ;---------------------------------------------------------------------------
   ; Loop through each input directory

   for imjd=0, N_elements(mjdlist)-1 do begin

      mjddir = mjdlist[imjd]
      inputdir = rawdir+'/'+mjddir
      plugdir = astrolog+'/'+mjddir
      splog, ''
      splog, 'Data directory ', inputdir
      splog, 'Astrolog directory ', plugdir

      ; Find all raw FITS files in this directory
      fullname = findfile(inputdir+'/sdR*.fit', count=nfile)
      splog, 'Number of FITS files found: ', nfile

      if (nfile GT 0) then begin

         ;----------
         ; Remove the path from the file names

         shortname = strarr(nfile)
         for i=0, nfile-1 do begin
            res = str_sep(fullname[i],'/')
            shortname[i] = res[N_elements(res)-1]
         endfor

         ;----------
         ; Sort the files based upon exposure number + camera number

;         isort = sort( strmid(shortname,7,8) + strmid(shortname,4,2) )
;         fullname = fullname[isort]
;         shortname = shortname[isort]

         ;----------
         ; Find all useful header keywords

         PLATEID = lonarr(nfile)
         EXPTIME = fltarr(nfile)
         EXPOSURE = fltarr(nfile)
         FLAVOR = strarr(nfile)
         CAMERAS = strarr(nfile)
         for i=0, nfile-1 do begin

            hdr = headfits(fullname[i])

            PLATEID[i] = long( sxpar(hdr, 'PLATEID') )
            EXPTIME[i] = sxpar(hdr, 'EXPTIME')
            EXPOSURE[i] = long( sxpar(hdr, 'EXPOSURE') )
            FLAVOR[i] = strtrim(sxpar(hdr, 'FLAVOR'),2)
            CAMERAS[i] = strmid(shortname[i],4,2) ; Camera number from file name

            ; Rename 'target' -> 'science', and 'calibration' -> 'arc'
            if (FLAVOR[i] EQ 'target') then FLAVOR[i] = 'science'
            if (FLAVOR[i] EQ 'calibration') then FLAVOR[i] = 'arc'

            ; Read the MJD from the header of the first file
            ; This should usually be the same as MJDDIR, though an integer
            ; rather than a string variable.
            if (i EQ 0) then thismjd = sxpar(hdr, 'MJD')
         endfor

         ;----------
         ; Determine all the plate numbers

         platenums = PLATEID[ uniq(PLATEID, sort(PLATEID)) ]
         nseq = N_elements(platenums)

         qdone = bytarr(nfile) ; Set equal to 1 as each seq
                               ; frame is written to a structure

         ;----------
         ; Loop through each plate

         for iseq=0, nseq-1 do begin

            ; Set plate number for this sequence (one sequence per plate)
            pltid = platenums[iseq]
            if (pltid GT 0 AND pltid LT 9999) then $
             platestr = string(pltid,format='(i04.4)') $
             else platestr = '0000'
            splog, camname=0, ''
            splog, camname='Plate '+platestr

            ; Zero-out data structure
            oneseq = 0

            ; Determine all files that are in this sequence
            ifile = where(pltid EQ PLATEID)

            ; Set the sequence ID equal to the first exposure number
            seqid = EXPOSURE[ifile[0]]

            ; Only look at those frames labelled as 'flat', 'arc', or 'science'
            ignore = where(FLAVOR[ifile] NE 'flat' $
                       AND FLAVOR[ifile] NE 'arc' $
                       AND FLAVOR[ifile] NE 'science', ct)
            if (ct GT 0) then qdone[ifile[ignore]] = 1
            if (ct GT 0) then $
             splog, 'Ignore ' + strtrim(string(ct),2) $
              + ' frames with unusable flavor'
;            for i=0, ct-1 do $
;             splog, 'Ignore file ', shortname[ifile[ignore[i]]], $
;              ' FLAVOR=', FLAVOR[ifile[ignore[i]]]

            ; Ignore short science exposures
            ignore = where(FLAVOR[ifile] EQ 'science' $
             AND EXPTIME[ifile] LT minexp, ct)
            if (ct GT 0) then qdone[ifile[ignore]] = 1
;            for i=0, ct-1 do $
;             splog, 'Ignore file ', shortname[ifile[ignore[i]]], $
;              ' EXPTIME=', EXPTIME[ifile[ignore[i]]]
            if (ct GT 0) then $
             splog, 'Ignore ' + strtrim(string(ct),2) $
              + ' science frames with EXPTIME < ', minexp

            while (min(qdone[ifile]) EQ 0) do begin
               inotdone = where(qdone[ifile] EQ 0)
               indx = ifile[inotdone]
               indx = indx[ where( EXPOSURE[indx] EQ EXPOSURE[indx[0]] $
                               AND FLAVOR[indx] EQ FLAVOR[indx[0]] ) ]
               oneexp = spplan_create_exp(seqid, pltid, FLAVOR[indx[0]])

               for icam=0, ncam-1 do begin
                  j = where(CAMERAS[indx] EQ camnames[icam], ct)
                  if (ct EQ 1) then begin
                     oneexp.name[icam] = shortname[indx[j[0]]]
                     qdone[indx[j[0]]] = 1
                  endif else if (ct GT 1) then begin
                     message, 'Several frames with EXPOSURE=' $
                      + string(EXPOSURE[indx[0]]) $
                      + ' and CAMERA=' + camnames[icam]
                  endif
               endfor

               ; Add all of these exposures to the ONESEQ structure
               if (NOT keyword_set(flats) OR $
                   (keyword_set(flats) AND oneexp.flavor EQ 'flat') ) then begin
                  if (keyword_set(oneseq)) then oneseq = [oneseq, oneexp] $
                   else oneseq = oneexp
               endif
            endwhile

            if (keyword_set(oneseq)) then begin
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
                   splog, 'No flat for SEQID=' + string(oneseq[0].seqid) $
                   + ' and PLATEID=' + string(oneseq[0].plateid)
                  if (nscience GT 0 AND narc EQ 0) then $
                   splog, 'No arc for SEQID= ' $
                    + strtrim(string(oneseq[0].seqid),2) $
                    + ', PLATEID= ' + strtrim(string(oneseq[0].plateid),2) $
                    + ', CAMERA= ' + strtrim(string(camnames[icam]),2)
               endfor

               ; Find a plug map file
               oneplug = spplan_create_plug(seqid, pltid)
               files = 'plPlugMapM-' + platestr + '*.par'
               files = findfile(filepath(files, root_dir=plugdir), count=ct)
               if (ct GT 1) then begin
                  splog, 'WARNING: Several plug map files found for plate number ' $
                   + platestr + ' in ' + plugdir
                  ; Use last one, because they should be ordered by MJD-rerun
                  files = files[ (sort(files))[ct-1]]
                  ct = 1
               endif
               if (ct EQ 1) then begin
                  ; Take the only plugmap file -- remove leading directory name
                  j = rstrpos(files[0], '/')
                  if (j EQ -1) then oneplug.name = files[0] $
                   else oneplug.name = strmid(files[0],j+1,strlen(files[0])-j)
                  splog, 'Setting plug map file = ', oneplug.name
               endif
               if (ct EQ 0) then begin
                  splog, 'WARNING: No plug map files found for plate number ' $
                   + platestr + ' in ' + plugdir
               endif

            endif

            ;----------
            ; Look for the pixel flats

            pixflats = spplan_create_pixflats()
            for icam=0, ncam-1 do begin
               files = 'pixflat-*-' + camnames[icam] + '.fits'
               files = findfile(filepath(files, root_dir=flatdir), count=ct)
               if (ct EQ 0) then begin
                 files = 'pixflat-*-' + camnames[icam] + '.fits'
                 files = findfile(filepath(files, root_dir=flatdir), count=ct)
               endif
               if (ct GT 1) then begin
                  splog, 'Several pixel flats found for CAMERA= ' + camnames[icam]
                  splog, 'Using pixel flat ', files[0]
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
                  splog, 'No pixel flats found for CAMERA= ' + camnames[icam]
               endif

            endfor

            ;----------
            ; Determine names of output files

            outdir = topdir + '/' + platestr
            if (strmid(flatdir,0,1) EQ '/') then fullflatdir = flatdir $
             else fullflatdir = '../../' + flatdir
            planfile = string( 'spPlan2d-', thismjd, '-', pltid, '.par', $
             format='(a,i5.5,a,i4.4,a)' )
            logfile = string( 'spDiag2d-', thismjd, '-', pltid, '.log', $
             format='(a,i5.5,a,i4.4,a)' )
            plotfile = string( 'spDiag2d-', thismjd, '-', pltid, '.ps', $
             format='(a,i5.5,a,i4.4,a)' )

            ;----------
            ; Create keyword pairs for plan file

            hdr = ''
            hdr = [hdr, "MJD     " + string(thismjd) + "  # Modified Julian Date"]
            hdr = [hdr, "inputDir    '" + inputdir + "'  # Directory for raw images"]
            hdr = [hdr, "plugDir     '" + plugdir + "'  # Directory for plugmap files"]
            hdr = [hdr, "flatDir     '" + fullflatdir + "'  # Directory for pixel flats"]
            hdr = [hdr, "extractDir  '2d'       # Directory for 2d spectra"]
            hdr = [hdr, "combineDir  '2dmerge'  # Directory for combined spectra"]
            hdr = [hdr, "logfile     '" + logfile + "'  # Text log file"]
            hdr = [hdr, "plotfile    '" + plotfile + "'  # PostScript log file"]

            ; Only output plan file if some raw FITS data files exist
            if (keyword_set(oneseq)) then begin
               spawn, 'mkdir -p '+outdir
               fullplanfile = filepath(planfile, root_dir=outdir)
               if (keyword_set(findfile(fullplanfile))) then $
                message, 'Output plan file already exists: '+planfile
               splog, 'Writing plan file ', fullplanfile
               yanny_write, fullplanfile, [ptr_new(pixflats), $
                ptr_new(oneplug), ptr_new(oneseq)], hdr=hdr
            endif

         endfor ; End loop through sequence number (one plate)

         splog, camname=0

      endif

   endfor ; End loop through input directory names (one MJD)

   return
end
;------------------------------------------------------------------------------
