;+
; NAME:
;   spplan1d
;
; PURPOSE:
;   Create plan file(s) for combining Spectro-2D outputs into one plate.
;
; CALLING SEQUENCE:
;   spplan1d, [ topindir=, topoutdir=, mjd=, mjstart=, mjend=, $
;    platenum=, platestart=, plateend=, /clobber ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   topindir   - Top directory name for 2D outputs; default to ''
;   topoutdir  - Top directory name for 2D outputs; default to the
;                same as TOPINDIR.
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

pro spplan1d, topindir=topindir, topoutdir=topoutdir, $
 mjd=mjd, mjstart=mjstart, mjend=mjend, $
 platenum=platenum, platestart=platestart, plateend=plateend, $
 clobber=clobber

   ;----------
   ; Determine the top-level of the input and output directory tree

   if (NOT keyword_set(topindir) OR NOT keyword_set(topoutdir)) then begin
      defaultdir = ''
   endif

   if (NOT keyword_set(topindir)) then topindir = defaultdir
   if (NOT keyword_set(topoutdir)) then topoutdir = topindir

   splog, 'Setting top-level of input directory to ' + topindir
   splog, 'Setting top-level of output directory to ' + topoutdir

   ;----------
   ; Create a list of the plate directories (as strings)

   platelist = get_mjd_dir(topindir, mjd=platenum, mjstart=platestart, $
    mjend=plateend)

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

      allplan = findfile(filepath('spPlan2d*.par', root_dir=topindir, $
       subdirectory=platedir), count=nplan)

      ;----------
      ; Read all the 2D plan files
      ; The string array PLANLIST keeps a list of the plan file that each
      ; element of the ALLEXP structure came from.

      allexp = 0
      planlist = 0
      extractdir = ''

      for iplan=0, nplan-1 do begin
         yanny_read, allplan[iplan], pp, hdr=hdr
         thismjd = yanny_par(hdr, 'MJD')
         for ii=0, n_elements(pp)-1 do begin
            sname = tag_names(*pp[ii], /structure_name)
            if (sname EQ 'SPEXP') then begin
               nadd = n_elements(*pp[ii])
               if (NOT keyword_set(allexp)) then begin
                  allexp = *pp[ii]
                  planlist = replicate(fileandpath(allplan[iplan]), nadd)
               endif else begin
                  allexp = [allexp, *pp[ii]]
                  planlist = [planlist, $
                   replicate(fileandpath(allplan[iplan]), nadd)]
               endelse
               if (NOT keyword_set(extractdir)) then $
                extractdir = yanny_par(hdr, 'extractdir')
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

            if (keyword_set(spexp)) then begin

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
                  kk = rstrpos(newnames[i], '.')
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
               outdir = concat_dir(topoutdir, platedir)

               planfile = 'spPlancomb-' + platestr + '-' + mjdstr + '.par'
               logfile = 'spDiagcomb-' + platestr + '-' + mjdstr + '.log'
               plotfile = 'spDiagcomb-' + platestr + '-' + mjdstr + '.ps'
               combinefile = 'spPlate-' + platestr + '-' + mjdstr + '.fits'
               fcalibprefix = 'spFluxcalib-' + platestr + '-' + mjdstr
               plotsnfile = 'spSN2d-' + platestr + '-' + mjdstr + '.ps'
               snfits = 'spSN2d-' + platestr + '-' + mjdstr + '.fits'

               ;----------
               ; Create keyword pairs for plan file

               hdr = ''
               hdr = [hdr, "planfile2d  '" + planlist1 $
                + "'  # Plan file for 2D spectral reductions"]
               hdr = [hdr, "planfilecomb '" + planfile $
                + "'  # Plan file for combining spectra"]
               hdr = [hdr, "extractdir '" + extractdir $
                + "'  # Directory for extracted spectra"]
               hdr = [hdr, "combinedir ''" $
                + "  # Directory for combined spectra"]
               hdr = [hdr, "logfile    '" + logfile $
                + "'  # Text log file"]
               hdr = [hdr, "plotfile   '" + plotfile $
                + "'  # PostScript log file"]
               hdr = [hdr, "fcalibprefix '" + fcalibprefix $
                + "'  # Prefix for flux-calibration files"]
               hdr = [hdr, "combinefile   '" + combinefile $
                + "'  # Output combined spectra file"]
               hdr = [hdr, "plotsnfile        '" + plotsnfile $
                + "'  # Two page S/N and magnitude plot"]
               hdr = [hdr, "snfile        '" + snfits $
                + "'  # Small fits file with S/N numbers "]

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
