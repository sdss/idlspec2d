;+
; NAME:
;   spplan2d
;
; PURPOSE:
;   Create plan file(s) for running the Spectro-2D pipeline.
;
; CALLING SEQUENCE:
;   spplan2d, [ topoutdir=, mjd=, mjstart=, mjend=, minexp=, $
;    /clobber ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   topoutdir  - Top directory name for output files; default to ''.
;   mjd        - Look for raw data files in RAWDATA_DIR/MJD; default to
;                search all subdirectories.  Note that this need not be
;                integer-valued, but could be for example '51441_test'.
;   mjstart    - Starting MJD.
;   mjend      - Ending MJD.
;   minexp     - Minimum exposure time for science frames; default to 300 sec.
;   clobber    - If set, then over-write conflicting plan files; default to
;                not over-write files.
;
; OUTPUT:
;
; COMMENTS:
;   The environment variable RAWDATA_DIR must be set.
;
;   Look for the raw FITS data files in:
;     RAWDATA_DIR/MJD/sdR-cs-eeeeeeee.fit
;   where c=color, s=spectrograph number, eeeeeeee=exposure number.
;
;   The output 2D plan files created are:
;     TOPOUTDIR/PLATE/spPlan2d-pppp-mmmmm.par
;   where pppp=plate number, mmmmm=MJD.
;
;   For the earliest data (before MJD=51455), then NAME keyword in the FITS
;   files did not properly describe the plug-map name.  In those cases,
;   look for the actual plug-map files in SPECLOG_DIR/MJD.
;
; EXAMPLES:
;   Create the plan file(s) for reducing the data for MJD=51441, with that
;   top level directory set to '/u/schlegel/2d_test'
;   > spplan2d, '/u/schlegel/2d_test', mjd=51441
;
; BUGS:
;   This routine spawns the Unix command 'mkdir'.
;   The use of CONCAT_DIR may not be valid for non-Unix OS's.
;
; PROCEDURES CALLED:
;   concat_dir()
;   fileandpath()
;   get_mjd_dir()
;   idlspec2d_version()
;   splog
;   sdsshead()
;   spplan_findrawdata
;   spplan_create_spexp
;   sxpar()
;   yanny_write
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro spplan2d, topoutdir=topoutdir, mjd=mjd, $
 mjstart=mjstart, mjend=mjend, minexp=minexp, clobber=clobber

   if (NOT keyword_set(minexp)) then minexp = 300

   ;----------
   ; Determine the top-level of the output directory tree

   if (NOT keyword_set(topoutdir)) then topoutdir = ''
   splog, 'Setting top-level of output directory to ' + topoutdir

   ;----------
   ; Read environment variable for RAWDATA_DIR for finding raw data files.

   rawdata_dir = getenv('RAWDATA_DIR')
   if (NOT keyword_set(rawdata_dir)) then $
    message, 'Must set environment variable RAWDATA_DIR'

   speclog_dir = getenv('SPECLOG_DIR')
   if (NOT keyword_set(speclog_dir)) then $
    message, 'Must set environment variable SPECLOG_DIR'
    
   ;----------
   ; Create a list of the MJD directories (as strings)

   mjdlist = get_mjd_dir(rawdata_dir, mjd=mjd, mjstart=mjstart, mjend=mjend)

   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   ;---------------------------------------------------------------------------
   ; Loop through each input MJD directory

   for imjd=0, N_elements(mjdlist)-1 do begin

      mjddir = mjdlist[imjd]
      inputdir = concat_dir(rawdata_dir, mjddir)
      plugdir = concat_dir(speclog_dir, mjddir)

      splog, ''
      splog, 'Data directory ', inputdir

      
      ; Find all raw FITS files in this directory
      fullname = spplan_findrawdata(inputdir, nfile)
      splog, 'Number of FITS files found: ', nfile

      if (nfile GT 0) then begin

         ;----------
         ; Remove the path from the file names

         shortname = strarr(nfile)
         for i=0, nfile-1 do shortname[i] = fileandpath(fullname[i])

         ;----------
         ; Find all useful header keywords

         PLATEID = lonarr(nfile)
         EXPTIME = fltarr(nfile)
         EXPOSURE = lonarr(nfile)
         FLAVOR = strarr(nfile)
         CAMERAS = strarr(nfile)
         MAPNAME = strarr(nfile)

         for i=0, nfile-1 do begin
            hdr = sdsshead(fullname[i])

            if (size(hdr,/tname) EQ 'STRING') then begin
               PLATEID[i] = long( sxpar(hdr, 'PLATEID') )
               EXPTIME[i] = sxpar(hdr, 'EXPTIME')
               EXPOSURE[i] = long( sxpar(hdr, 'EXPOSURE') )
               FLAVOR[i] = strtrim(sxpar(hdr, 'FLAVOR'),2)
               CAMERAS[i] = strtrim(sxpar(hdr, 'CAMERAS'),2)
               MAPNAME[i] = strtrim(sxpar(hdr, 'NAME'),2)

               ; Read the MJD from the header of the first file
               ; This should usually be the same as MJDDIR, though an integer
               ; rather than a string variable.
               if (i EQ 0) then thismjd = sxpar(hdr, 'MJD')

               ; MAPNAME should be of the form '0306-51683-01'.
               ; If it only contains the PLATEID (for MJD <= 51454),
               ; then find the actual plug-map file.
               if (strlen(MAPNAME[i]) LE 4) then begin
                  plugfile = 'plPlugMapM-' $
                   + string(long(MAPNAME[i]), format='(i4.4)') + '-*.par'
                  plugfile = (findfile(filepath(plugfile, root_dir=plugdir), $
                   count=ct))[0]
                  if (ct EQ 1) then $
                   MAPNAME[i] = strmid(fileandpath(plugfile), 11, 13)
               endif
            endif
         endfor

         ;----------
         ; Determine all the plate plugging names

         allmaps = MAPNAME[ uniq(MAPNAME, sort(MAPNAME)) ]

         ;----------
         ; Loop through all plate plugging names

         for imap=0, n_elements(allmaps)-1 do begin

            spexp = 0 ; Zero-out this output structure

            ;----------
            ; Loop through all exposure numbers for this plate-plugging

            theseexp = EXPOSURE[ where(MAPNAME EQ allmaps[imap]) ]
            allexpnum = theseexp[ uniq(theseexp, sort(theseexp)) ]

            for iexp=0, n_elements(allexpnum)-1 do begin
               indx = where(MAPNAME EQ allmaps[imap] $
                AND EXPOSURE EQ allexpnum[iexp])

               spexp1 = spplan_create_spexp(allexpnum[iexp], $
                PLATEID[indx[0]], thismjd, $
                MAPNAME[indx[0]], FLAVOR[indx[0]], EXPTIME[indx[0]], $
                shortname[indx], CAMERAS[indx], minexp=minexp)

               if (keyword_set(spexp1)) then begin
                  if (keyword_set(spexp)) then spexp = [spexp, spexp1] $
                   else spexp = spexp1
               endif
            endfor

            ;----------
            ; Discard these observations if the plate number is not
            ; in the range 1 to 9990.

            pltid = PLATEID[indx[0]]
            if (pltid GT 0 AND pltid LT 9990) then begin
               platestr = string(pltid, format='(i04.4)')
            endif else begin
               platestr = '0000'
               spexp = 0
            endelse
            mjdstr = string(thismjd, format='(i05.5)')

            ;----------
            ; Discard these observations if there is not at least one flat,
            ; one arc, and one science exposure

            if (keyword_set(spexp)) then begin
               junk = where(spexp.flavor EQ 'flat', ct)
               if (ct EQ 0) then begin
                  splog, 'WARNING: No flats for MAPNAME=' + allmaps[imap]
                  spexp = 0
               endif
            endif

            if (keyword_set(spexp)) then begin
               junk = where(spexp.flavor EQ 'arc', ct)
               if (ct EQ 0) then begin
                  splog, 'WARNING: No arcs for MAPNAME=' + allmaps[imap]
                  spexp = 0
               endif
            endif

            if (keyword_set(spexp)) then begin
               junk = where(spexp.flavor EQ 'science', ct)
               if (ct EQ 0) then begin
                  splog, 'WARNING: No science frames for MAPNAME=' + allmaps[imap]
                  spexp = 0
               endif
            endif

            if (keyword_set(spexp)) then begin

               ;----------
               ; Determine names of output files

               outdir = concat_dir(topoutdir, platestr)

               planfile = 'spPlan2d-' + platestr + '-' + mjdstr + '.par'
               logfile = 'spDiag2d-' + platestr + '-' + mjdstr + '.log'
               plotfile = 'spDiag2d-' + platestr + '-' + mjdstr + '.ps'

               ;----------
               ; Create keyword pairs for plan file

               hdr = ''
               hdr = [hdr, "plateid  " + platestr + "  # Plate number"]
               hdr = [hdr, "MJD     " + mjdstr $
                + "  # Modified Julian Date"]
               hdr = [hdr, "planfile2d  '" + planfile $
                + "'  # Plan file for 2D spectral reductions"]
               hdr = [hdr, "extractdir ''" $
                + "  # Directory for spFrame files"]
               hdr = [hdr, "logfile    '" $
                + logfile + "'  # Text log file"]
               hdr = [hdr, "plotfile   '" + plotfile $
                + "'  # PostScript log file"]

               ;----------
               ; Write output file

               spawn, 'mkdir -p ' + outdir
               fullplanfile = filepath(planfile, root_dir=outdir)
               qexist = keyword_set(findfile(fullplanfile))
               if (qexist) then begin
                  if (keyword_set(clobber)) then $
                   splog, 'WARNING: Over-writing plan file: ' + planfile $
                  else $
                   splog, 'WARNING: Will not over-write plan file: ' + planfile
               endif
               if ((NOT qexist) OR keyword_set(clobber)) then begin
                  splog, 'Writing plan file ', fullplanfile
                  yanny_write, fullplanfile, ptr_new(spexp), hdr=hdr
               endif
            endif

         endfor ; End loop through plate plugging names
      endif
   endfor

   return
end
;------------------------------------------------------------------------------
