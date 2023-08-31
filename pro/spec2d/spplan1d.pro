;+
; NAME:
;   spplan1d
;
; PURPOSE:
;   Create plan file(s) for combining Spectro-2D outputs.
;
; CALLING SEQUENCE:
;   spplan1d, [ topdir=, run2d=, mjd=, mjstart=, mjend=, $
;    platenum=, platestart=, plateend=, /clobber ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   topdir     - Optional override value for the environment
;                variable $BHM_SPECTRO_REDUX.
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
;   15-Nov-2018  Modified by Hector Ibarra for the BHM
;-
;------------------------------------------------------------------------------


pro spplan1d, topdir=topdir1, run2d=run2d1, $
 mjd=mjd, mjstart=mjstart, mjend=mjend, $
 confinum=confinum, confistart=confistart, confiend=confiend, $
 fieldnum=fieldnum, fieldstart=fieldstart, fieldend=fieldend, $
 clobber=clobber, lco=lco, legacy=legacy, plates=plates

 RESOLVE_ALL, /QUIET, /SKIP_EXISTING, /CONTINUE_ON_ERROR

   ;----------
   ; Determine the top-level of the directory tree
   ;----------
   if (keyword_set(topdir1)) then topdir = topdir1 else topdir = getenv('BOSS_SPECTRO_REDUX')
   splog, 'Setting TOPDIR=', topdir

   if (keyword_set(run2d1)) then run2d = strtrim(run2d1,2) else run2d = getenv('RUN2D')
   splog, 'Setting RUN2D=', run2d
   if (keyword_set(run2d)) then $
        topdir = djs_filepath('', root_dir=topdir, subdir=run2d)

   finfo = file_info(topdir)
   if (finfo.directory EQ 0) then begin
      splog, 'Directory does not exist: '+topdir
      return
   endif
    
    

   fieldlist = get_mjd_dir(topdir, mjd=fieldnum, mjstart=fieldstart, $
      mjend=fieldend)
   splog, 'Number of field directories = ', n_elements(fieldlist)


   ; Loop through each input configuration directory
   for ifield=0, N_elements(fieldlist)-1 do begin
      fielddir = fieldlist[ifield]
      fieldid = long(fileandpath(fielddir))
      get_field_type, fieldid=fieldid, legacy=legacy, plates=plates, fps=fps
      splog, ''
      splog, 'Field directory ', fielddir
      ;----------
      ; Find all 2D plan files
      allplan = findfile(djs_filepath('spPlan2d*.par', root_dir=topdir, $
           subdirectory=fielddir), count=nplan)
      ;----------
      ; Read all the 2D plan files
      ; The string array PLANLIST keeps a list of the plan file that each
      ; element of the ALLEXP structure came from, and MJDLIST keeps the
      ; list of each MJD
      allexp = 0
      planlist = 0
      ditherlist = 0
      for iplan=0, nplan-1 do begin
             yanny_read, allplan[iplan], pp, hdr=hdr
             thismjd = long(yanny_par(hdr, 'MJD'))
             thisdither = yanny_par(hdr, 'DITHER')
             for ii=0, n_elements(pp)-1 do begin
                sname = tag_names(*pp[ii], /structure_name)
                if (sname EQ 'SPEXP') then begin
                   nadd = n_elements(*pp[ii])
                   if (NOT keyword_set(allexp)) then begin
                      allexp = *pp[ii]
                      planlist = replicate(fileandpath(allplan[iplan]), nadd)
                      mjdlist = replicate(thismjd, nadd)
                      ditherlist = replicate(thisdither, nadd)
                   endif else begin
                      allexp = [allexp, *pp[ii]]
                      planlist = [planlist, $
                       replicate(fileandpath(allplan[iplan]), nadd)]
                      mjdlist = [mjdlist, replicate(thismjd, nadd)]
                      ditherlist = [ditherlist, replicate(thisdither, nadd)]
                   endelse
                endif
             endfor
             yanny_free, pp
      endfor
    
      if keyword_set(legacy) or keyword_set(plates) then begin
         if (keyword_set(allexp)) then begin
           ;----------
           ; Determine all the plate plugging names
           allmaps = allexp.mapname
           allmaps = allmaps[ uniq(allmaps, sort(allmaps)) ]
           ;----------
           ; Loop through all plate plugging names
           ;print,allmaps
           ;stop
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
               ditherlist1 = ditherlist[indx]
               ;----------
               ; Replace the prefix 'sdR' with 'spFrame' in the science frames
               ; and the suffix '.fit' with '.fits'
               ;HJIM epoch definition for the plate program
               RM_plates=[15000,15001,15002,15038,15070,15071,15171,15172,15173,15252,15253]
               plt=spexp[0].fieldid
               junk = where(RM_plates eq plt ,ct)
               epc=1000;this number force to coadd during all the pluggin mapname
               if ct gt 0 then epc=3
               spexp_org=spexp
               stop_loop=0
               while stop_loop eq 0 do begin
                   indx_ep = where(spexp.mjd lt min(spexp.mjd)+epc)
                   spexp=spexp[indx_ep]
                   epoch_tag = replicate( $
                         {epoch_combine: min(spexp.mjd)}, n_elements(spexp))
                   spexp = struct_addtags(spexp, epoch_tag)
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
                   pltid = spexp[0].fieldid
                   platestr = field_to_string(pltid)
                   thismjd = max(spexp.mjd)
                   mjdstr = string(thismjd, format='(i05.5)')
                   outdir = concat_dir(topdir, fielddir)
                   planfile = 'spPlancomb-' + platestr + '-' + mjdstr + '.par'
                   ;----------
                   ; Create keyword pairs for plan file
                   hdr = ''
                   hdr = [hdr, "fieldid  " + platestr + "  # Field number"]
                   hdr = [hdr, "MJD      " + mjdstr $
                     + "  # Modified Julian Date for most recent observation"]
                   hdr = [hdr, "OBS     " + "APO  # Observatory"]
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
                   hdr = [hdr, "DITHER '"+ ditherlist1[0] $
                    + " # Is the Field Dithered (T: True, F: False)"]
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
                     yanny_write, fullplanfile, ptr_new(spexp), hdr=hdr, stnames='SPEXP'
                   endif
                   if max(spexp_org.mjd) eq max(spexp.mjd) then begin
                     stop_loop=1
                   endif else begin
                     stop_loop=0
                   endelse
                   indx_ep = where(spexp_org.mjd gt max(spexp.mjd))
                   if (indx_ep[0] NE -1) then spexp = spexp_org[indx_ep]
               endwhile
             endif
           endfor ; End loop through plate plugging names
         endif
      endif else begin
          if (keyword_set(allexp)) then begin
                indx = where((allexp.flavor EQ 'science' OR allexp.flavor EQ 'smear'))
                if (indx[0] NE -1) then begin
                  spexp = allexp[indx]
                  planlist = planlist[indx]
                  ditherlist1 = ditherlist[indx]
                endif else begin
                  spexp = 0
                endelse
                ;----------
                ; Decide if any of these MJD's are within the bounds
                ; specified by MJD,MJSTART,MJEND.  If so, then set QMJD=1
                qmjd = 1B
                if (keyword_set(spexp)) then begin
                   mjdlist1 = mjdlist[indx]
                   qmjd = mjd_match(mjdlist1, mjd=mjd, mjstart=mjstart, mjend=mjend)
                   if (qmjd EQ 0) then $
                    splog, 'Skip FIELD=', spexp[0].fieldid, ' with MJD=', $
                     mjdlist1[ uniq(mjdlist1, sort(mjdlist1)) ]
                endif
                if (keyword_set(spexp) AND qmjd) then begin
                  RM_fields=[]
                  fld=spexp[0].fieldid
                  junk = where(RM_fields eq fld ,ct)
                  epc=1 ; only coadd single night
                  if ct gt 0 then epc=3
                  spexp_org=spexp
                  planlist_org=planlist
                  stop_loop=0
                  while stop_loop eq 0 do begin
                   indx_ep = where(spexp.mjd lt min(spexp.mjd)+epc)
                   spexp=spexp[indx_ep]
                   epoch_tag = replicate( $
                         {epoch_combine: min(spexp.mjd)}, n_elements(spexp))
                   spexp = struct_addtags(spexp, epoch_tag)
  
                   ;----------
                   ; Determine the 2D plan files that are relevant
                   planlist1 = planlist[indx_ep]
                   ;----------
                   ; Replace the prefix 'sdR' with 'spFrame' in the science frames
                   ; and the suffix '.fit' with '.fits'
                   if keyword_set(mjd) then begin
                      runspexp=[]
                      planlist1f=[]
                      foreach runmjd, mjd do begin
                          nind=where(spexp.mjd EQ mjd,ct)
                          if ct ne 0 then begin
                              runspexp=[runspexp,spexp[nind]]
                              planlist1f=[planlist1f,planlist1[nind]]
                          endif
                      endforeach
                      if n_elements(planlist1f) eq 0 then goto, NoMatch
                      planlist1=planlist1f[0]
                   endif
                   if (keyword_set(mjstart) AND keyword_set(mjend) AND qmjd) then begin
                     nind=where((spexp.mjd LE mjend) and (spexp.mjd GE mjstart))
                     spexp=spexp[nind]
                     planlist1f=planlist1[nind]
                     planlist1=planlist1f[0]
                   endif
                   planlist1 = planlist1[ uniq(planlist1, sort(planlist1)) ]
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
                   fieldstr = spexp[0].fieldid
                   thismjd = max(spexp.mjd)
                   mjdstr = string(thismjd, format='(i05.5)')
                   outdir = concat_dir(topdir, fielddir)
                   planfile = 'spPlancomb-' + fieldstr + '-' + mjdstr + '.par'
                   ;----------
                   ; Create keyword pairs for plan file
                   hdr = ''
                   hdr = [hdr, "fieldid  " + fieldstr + "  # Field number"]
                   hdr = [hdr, "MJD      " + mjdstr $
                    + "  # Modified Julian Date for most recent observation"]
                   if keyword_set(lco) then $
                      hdr = [hdr, "OBS     " + "LCO  # Observatory"] $
                   else hdr = [hdr, "OBS     " + "APO  # Observatory"]
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
                   hdr = [hdr, "DITHER '"+ ditherlist1[0] $
                    + " # Is the Field Dithered (T: True, F: False)"]
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
                      yanny_write, fullplanfile, ptr_new(spexp), hdr=hdr, stnames='SPEXP'
                   endif
                   NoMatch: 
                   if max(spexp_org.mjd) eq max(spexp.mjd) then begin
                     stop_loop=1
                   endif else begin
                     stop_loop=0
                   endelse
                   indx_ep = where(spexp_org.mjd gt max(spexp.mjd))
                   if (indx_ep[0] NE -1) then begin
                      spexp = spexp_org[indx_ep]
                      planlist=planlist_org[indx_ep]
                   endif
                  endwhile
                endif
          endif
      endelse
    endfor
   return
end
;------------------------------------------------------------------------------
