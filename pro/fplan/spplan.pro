;+
; NAME:
;   spplan
;
; PURPOSE:
;   Create plan file(s) for running the Spectro-2D pipeline.
;
; CALLING SEQUENCE:
;   spplan, [ topindir=, topoutdir=, mjd=, mjstart=, mjend=, minexp=, $
;    flatdir=, /clobber ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   topindir   - Top directory name for input files;
;                default to '/data/spectro'.
;   topoutdir  - Top directory name for output files; default to the
;                subdirectory '2d_' + VERSION under the current directory.
;   flatdir    - override default pixflats directory  $TOPINDIR/pixflats
;   mjd        - Look for raw data files in RAWDIR/MJD; default to
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
;   Look for the input files in:
;     TOPINDIR/rawdata/MJD/sdR-cs-eeeeeeee.fit           - Raw frames
;     TOPINDIR/astrolog/MJD/plPlugMapM-pppp-mmmmm-aa.par - Plug map files
;     TOPINDIR/pixflats/pixflat-mmmmm-cs.fits               - Pixel flats
;   where c=color, s=spectrograph number, pppp=plate number, aa=mapper ID,
;   mmmmm=MJD.
;
;   The top-level of the output directory structure, TOPOUTDIR, is 2d_VERSION,
;   where we read the version with IDLSPEC2D_VERSION().  If this is not a CVS-
;   declared version of the code, then the output directory is '2d_test'.
;   The files created are:
;     TOPOUTDIR/PLATE/spPlan2d-pppp-mmmmm.par  - 2D plan file (could be several)
;
; EXAMPLES:
;   Create the plan file(s) for reducing the data for MJD=51441, with that
;   top level directory set to '/u/schlegel/spectro'
;   > spplan, '/usr/sdss/data05/spectro', mjd=51441
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
;   sxpar()
;   yanny_write
;
; INTERNAL SUPPORT ROUTINES:
;   spplan_create_exp
;   spplan_create_plug
;   spplan_create_pixflats
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

function spplan_create_spexp, expnum, plateid, mjd, mapname, flavor, exptime, $
 filename, cameras, minexp=minexp

   if (flavor NE 'flat' AND flavor NE 'arc' $
    AND flavor NE 'science' AND flavor NE 'smear') then $
    return, 0
   if (keyword_set(minexp)) then begin
      if (flavor EQ 'science' AND exptime LT minexp) then $
       return, 0 
   endif

   badname = 'UNKNOWN'
   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   spexp = {spexp, $
    plateid : long(plateid), $
    mjd     : long(mjd), $
    mapname : string(mapname), $
    flavor  : string(flavor), $
    exptime : float(exptime), $
    name    : strarr(ncam) }

   for icam=0, ncam-1 do begin
      ii = where(cameras EQ camnames[icam], ct)
      if (ct GT 1) then $
       message, 'Multiple files with EXPOSURE=' + string(expnum) $
        + ' CAMERAS=' + camnames[icam]
      if (ct EQ 1) then $
       spexp.name[icam] = filename[ii[0]] $
      else $
       spexp.name[icam] = badname
   endfor

   return, spexp
end

;------------------------------------------------------------------------------

pro spplan, topindir=topindir, topoutdir=topoutdir, mjd=mjd, $
 mjstart=mjstart, mjend=mjend, minexp=minexp, clobber=clobber, $
 flatdir=flatdir

;   if (NOT keyword_set(topindir)) then topindir = '/data/spectro'
   if (NOT keyword_set(topindir)) then topindir = '/home/data'
   if (NOT keyword_set(minexp)) then minexp = 300

   ;----------
   ; Determine the top-level of the output directory tree

   if (NOT keyword_set(topoutdir)) then begin
      vers = idlspec2d_version()
      if (strpos(vers, 'NOCVS') NE -1) then vers = 'test'
      topoutdir = '2d_' + vers
   endif
   splog, 'Setting top-level of output directory to ' + topoutdir

   ;----------
   ; Set directory names RAWDIR, ASTROLOG, FLATDIR

   rawdir = concat_dir(topindir, 'rawdata')
   astrolog = concat_dir(topindir, 'astrolog')

   if NOT keyword_set(flatdir) then  $
   flatdir = concat_dir(topindir, 'pixflats')

   ;----------
   ; Create a list of the MJD directories (as strings)

   mjdlist = get_mjd_dir(rawdir, mjd=mjd, mjstart=mjstart, mjend=mjend)

   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   ;---------------------------------------------------------------------------
   ; Loop through each input MJD directory

   for imjd=0, N_elements(mjdlist)-1 do begin

      mjddir = mjdlist[imjd]
      inputdir = concat_dir(rawdir, mjddir)
      plugdir = concat_dir(astrolog, mjddir)

      splog, ''
      splog, 'Data directory ', inputdir
      splog, 'Astrolog directory ', plugdir

      ; Find all raw FITS files in this directory
      fullname = findfile(filepath('sdR*.fit', root_dir=inputdir), count=nfile)

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
               hdr = [hdr, "MJD     " + mjdstr $
                + "  # Modified Julian Date"]
               hdr = [hdr, "inputdir   '" + inputdir $
                + "'  # Directory for raw images"]
               hdr = [hdr, "plugdir    '" + plugdir $
                + "'  # Directory for plugmap files"]
               hdr = [hdr, "flatdir    '" + flatdir $
                + "'  # Directory for pixel flats"]
               hdr = [hdr, "extractdir ''" $
                + "  # Directory for extracted spectra"]
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
