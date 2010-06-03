;+
; NAME:
;   spplan1d
;
; PURPOSE:
;   Create plan file(s) for combining Spectro-2D outputs into one plate.
;
; CALLING SEQUENCE:
;   spplan1d, [ topdir=, run2d=, mjd=, mjstart=, mjend=, $
;    platenum=, platestart=, plateend=, /clobber ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   topdir     - Optional override value for the environment
;                variable $BOSS_SPECTRO_REDUX.
;   run2d      - Optional override value for the environment variable $RUN2D
;   mjd        - Use data from these MJD's.
;   mjstart    - Starting MJD.
;   mjend      - Ending MJD.
;   platenum   - Look for input data files in TOPINDIR/PLATENUM; default to
;                search all subdirectories.  Note that this need not be
;                integer-valued, but could be for example '0306_test'.
;   platestart - Starting plate number.
;   plateend   - Ending plate number.
;   clobber    - If set, then over-write conflicting plan files; default to
;                not over-write files.
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
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
;   idlutils_version()
;   mjd_match()
;   splog
;   sdsshead()
;   sxpar()
;   yanny_free
;   yanny_par()
;   yanny_read
;   yanny_write
;
; REVISION HISTORY:
;   04-Jul-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro spplan1d, topdir=topdir1, run2d=run2d1, $
 mjd=mjd, mjstart=mjstart, mjend=mjend, $
 platenum=platenum, platestart=platestart, plateend=plateend, $
 clobber=clobber

   ;----------
   ; Determine the top-level of the directory tree

   if (keyword_set(topdir1)) then topdir = topdir1 $
    else topdir = getenv('BOSS_SPECTRO_REDUX')
   splog, 'Setting TOPDIR=', topdir
   if (keyword_set(run2d1)) then run2d = strtrim(run2d1,2) $
    else run2d = getenv('RUN2D')
   splog, 'Setting RUN2D=', run2d
   if (keyword_set(run2d)) then $
    topdir = djs_filepath('', root_dir=topdir, subdir=run2d)
   finfo = file_info(topdir)
   if (finfo.directory EQ 0) then begin
      splog, 'Directory does not exist: '+topdir
      return
   endif

   ;----------
   ; Create a list of the plate directories (as strings)

   platelist = get_mjd_dir(topdir, mjd=platenum, mjstart=platestart, $
    mjend=plateend)
   splog, 'Number of plate directories = ', n_elements(platelist)

   camnames = ['b1', 'b2', 'r1', 'r2']
   ncam = N_elements(camnames)

   ;---------------------------------------------------------------------------
   ; Loop through each input plate directory

   for iplate=0, N_elements(platelist)-1 do begin

      platedir = platelist[iplate]

      splog, ''
      splog, 'Plate directory ', platedir

      ;----------
      ; Find all 2D plan files

      allplan = findfile(djs_filepath('spPlan2d*.par', root_dir=topdir, $
       subdirectory=platedir), count=nplan)

      ;----------
      ; Read all the 2D plan files
      ; The string array PLANLIST keeps a list of the plan file that each
      ; element of the ALLEXP structure came from, and MJDLIST keeps the
      ; list of each MJD

      allexp = 0
      planlist = 0

      for iplan=0, nplan-1 do begin
         yanny_read, allplan[iplan], pp, hdr=hdr
         thismjd = long(yanny_par(hdr, 'MJD'))
         for ii=0, n_elements(pp)-1 do begin
            sname = tag_names(*pp[ii], /structure_name)
            if (sname EQ 'SPEXP') then begin
               nadd = n_elements(*pp[ii])
               if (NOT keyword_set(allexp)) then begin
                  allexp = *pp[ii]
                  planlist = replicate(fileandpath(allplan[iplan]), nadd)
                  mjdlist = replicate(thismjd, nadd)
               endif else begin
                  allexp = [allexp, *pp[ii]]
                  planlist = [planlist, $
                   replicate(fileandpath(allplan[iplan]), nadd)]
                  mjdlist = [mjdlist, replicate(thismjd, nadd)]
               endelse
            endif
         endfor
         yanny_free, pp
      endfor

      if (keyword_set(allexp)) then begin

         ;----------
         ; Determine all the plate plugging names

         allmaps = allexp.mapname
         allmaps = allmaps[ uniq(allmaps, sort(allmaps)) ]

         ;----------
         ; Loop through all plate plugging names

         for imap=0, n_elements(allmaps)-1 do begin

            indx = where(allexp.mapname EQ allmaps[imap] $
             AND (allexp.flavor EQ 'science' OR allexp.flavor EQ 'smear'))
            if (indx[0] NE -1) then spexp = allexp[indx] $
             else spexp = 0

            ;----------
            ; Decide if any of these MJD's are within the bounds
            ; specified by MJD,MJSTART,MJEND.  If so, then set QMJD=1

            qmjd = 1B
            if (keyword_set(spexp)) then begin
               mjdlist1 = mjdlist[indx]
               qmjd = mjd_match(mjdlist1, mjd=mjd, mjstart=mjstart, mjend=mjend)

               if (qmjd EQ 0) then $
                splog, 'Skip MAP=', allmaps[imap], ' with MJD=', $
                 mjdlist1[ uniq(mjdlist1, sort(mjdlist1)) ]
            endif

            if (keyword_set(spexp) AND qmjd) then begin

               ;----------
               ; Determine the 2D plan files that are relevant

               planlist1 = planlist[indx]
               planlist1 = planlist1[ uniq(planlist1, sort(planlist1)) ]

               ;----------
               ; Replace the prefix 'sdR' with 'spFrame' in the science frames
               ; and the suffix '.fit' with '.fits'

               newnames = spexp.name
               for i=0, n_elements(newnames)-1 do begin
                  jj = strpos(newnames[i], '-')
                  kk = strpos(newnames[i], '.', /reverse_search)
                  if (jj NE -1 AND kk NE -1) then $
                   newnames[i] = 'spFrame' + strmid(newnames[i], jj, kk-jj) $
                    + '.fits'
               endfor
               spexp.name = newnames

               ;----------
               ; Determine names of output files

               pltid = spexp[0].plateid
               platestr = string(pltid, format='(i04.4)')
               thismjd = max(spexp.mjd)
               mjdstr = string(thismjd, format='(i05.5)')
               outdir = concat_dir(topdir, platedir)

               planfile = 'spPlancomb-' + platestr + '-' + mjdstr + '.par'

               ;----------
               ; Create keyword pairs for plan file

               hdr = ''
               hdr = [hdr, "plateid  " + platestr + "  # Plate number"]
               hdr = [hdr, "MJD      " + mjdstr $
                + "  # Modified Julian Date for most recent observation"]
               hdr = [hdr, "RUN2D  " + run2d + "  # 2D reduction name"]
               sq = "'"
               hdr = [hdr, "planfile2d  " $
                + string(sq+planlist1+sq+' ', format='(99a)') $
                + " # Plan file(s) for 2D spectral reductions"]
               hdr = [hdr, "planfilecomb '" + planfile $
                + "'  # Plan file for combine (this file)"]
               hdr = [hdr, "idlspec2dVersion '" + idlspec2d_version() $
                + "'  # Version of idlspec2d when building plan file"]
               hdr = [hdr, "idlutilsVersion '" + idlutils_version() $
                + "'  # Version of idlutils when building plan file"]
               spawn, 'speclog_version', logvers, /noshell
               hdr = [hdr, "speclogVersion '" + logvers $
                + "'  # Version of speclog when building plan file"]

               ;----------
               ; Write output file

               spawn, 'mkdir -p ' + outdir
               fullplanfile = djs_filepath(planfile, root_dir=outdir)
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
