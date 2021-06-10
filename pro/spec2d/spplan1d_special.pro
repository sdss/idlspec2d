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
pro spplan1d_special, topdir=topdir1, run2d=run2d1, $
 mjd=mjd, mjstart=mjstart, mjend=mjend, $
 confinum=confinum, confistart=confistart, confiend=confiend, $
 fieldnum=fieldnum, fieldstart=fieldstart, fieldend=fieldend, $
 clobber=clobber, lco=lco, legacy=legacy, plates=plates, $
 platnum=platnum, platstart=platstart, platend=platend

   ;----------
   ; Determine the top-level of the directory tree
   ; 
   if keyword_set(lco) then begin
     obsdir='LCO'
     BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_S'
   endif else begin
     obsdir='APO'
     BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_N'
   endelse
   obsdir='';coment this line for the final version HJIM
   ;----------
   if (keyword_set(topdir1)) then topdir = topdir1 $
    else topdir = getenv('BOSS_SPECTRO_REDUX')
   topdir=concat_dir(topdir, obsdir)
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
   
   ;reduxdir = getenv('BOSS_SPECTRO_REDUX')
   topdir2=topdir+'/eFEDS'
   path=topdir+'/conflist.fits'
   data=mrdfits(path,1)
   ind=uniq(data.plate)
   plates_f=data[ind].plate
   program=data[ind].PROGRAMNAME
   ind=where(strmatch(program, '*eFED*', /fold_case) eq 1,nplates)
   plates_f=plates_f[ind]
   platelist=plate_to_string(plates_f)
   ;if keyword_set(platstart) or keyword_set(platend) then begin
   ; if (not keyword_set(legacy)) and (not keyword_set(plates)) then begin
   ;  splog, 'To reduce the plate format you must set the /legacy or /plates key word, process halt'
   ;  exit
   ; endif
   ;endif
    
   if keyword_set(legacy) or keyword_set(plates) then begin
     ;----------
     ; Create a list of the plate directories (as strings) and select only the plate directories
     ;platelist = get_mjd_dir(topdir, mjd=platnum, mjstart=platstart, $
     ;  mjend=platend,/alldirs)
       for ili=0, n_elements(platelist)-1 do begin
        platelist[ili]=platelist[ili]+'p'
        if strmid(strtrim(platelist[ili],2),4,1) ne 'p' then begin
          if strmid(strtrim(platelist[ili],2),5,1) ne 'p' then begin
            platelist[ili]=''
          endif
        endif
       endfor
     ii = where(platelist NE '', ct)
     if (ct EQ 0) then begin
         splog, 'There is no plate directories'
         return
     endif else begin
        platelist = platelist[ii]
     endelse
     splog, 'Number of plate directories = ', n_elements(platelist)
   endif else begin
     ;----------
     ; Create a list of the configuration directories (as strings)
     ; HJIM-- Change fiber by confi
     ;fieldlist = get_mjd_dir(topdir, mjd=fieldnum, mjstart=fieldstart, $
     ; mjend=fieldend)
     ;splog, 'Number of bhm field directories = ', n_elements(fieldlist)
   endelse
  
   ;HJIM -- reduce the number of spectrographs to one
   camnames = ['b1', 'r1']
   ncam = N_elements(camnames)
   if keyword_set(legacy) then begin
     camnames = ['b1', 'r1', 'b2', 'r2']
   endif


   if keyword_set(legacy) or keyword_set(plates) then begin
     ;---------------------------------------------------------------------------
     ; Loop through each input plate directory
     spawn, 'mkdir -p ' + topdir2
     for iplate=0, N_elements(platelist)-1 do begin
       platedir = '00000p';
       platedirt = platelist[iplate]
       splog, ''
       splog, 'Plate directory ', platedirt;,platelist[iplate]
       outdir = concat_dir(topdir2, platedir)
       spawn, 'mkdir -p ' + outdir
       ;----------
       ; Find all 2D plan files
       allplan = findfile(djs_filepath('spPlan2d*.par', root_dir=topdir, $
         subdirectory=platedirt), count=nplan)
       for ifile=0, nplan-1 do begin
         ;print,allplan[ifile]
         spawn, 'cp ' + allplan[ifile] + '  ' + outdir + '/'
       endfor
       tempf = findfile(djs_filepath('spArc-*', root_dir=topdir, $
         subdirectory=platedirt), count=ntemp)
       for ifile=0, ntemp-1 do begin
         spawn, 'cp ' + tempf[ifile] + '  ' + outdir + '/'
       endfor
       tempf = findfile(djs_filepath('spDiag2d*', root_dir=topdir, $
         subdirectory=platedirt), count=ntemp)
       for ifile=0, ntemp-1 do begin
         spawn, 'cp ' + tempf[ifile] + '  ' + outdir + '/'
       endfor
       tempf = findfile(djs_filepath('spFlat-*', root_dir=topdir, $
         subdirectory=platedirt), count=ntemp)
       for ifile=0, ntemp-1 do begin
         spawn, 'cp ' + tempf[ifile] + '  ' + outdir + '/'
       endfor
       tempf = findfile(djs_filepath('spFrame-*', root_dir=topdir, $
         subdirectory=platedirt), count=ntemp)
       for ifile=0, ntemp-1 do begin
         spawn, 'cp ' + tempf[ifile] + '  ' + outdir + '/'
       endfor
     endfor   
     ;for iplate=0, N_elements(platelist)-1 do begin   
       allplan = findfile(djs_filepath('spPlan2d*.par', root_dir=topdir2, $
         subdirectory=platedir), count=nplan)
       ;----------
       ; Read all the 2D plan files
       ; The string array PLANLIST keeps a list of the plan file that each
       ; element of the ALLEXP structure came from, and MJDLIST keeps the
       ; list of each MJD
       allexp = 0
       planlist = 0
       for iplan=0, nplan-1 do begin
         temp_plan=strsplit(allplan[iplan],'/',/extract)
         temp_mjd=strsplit(temp_plan[n_elements(temp_plan)-1],'-',/extract)
         temp_mjd=long(repstr(temp_mjd[n_elements(temp_mjd)-1],'.par'))
         if keyword_set(plates) then begin
           min_mjd=59030;LIMIT the use of single spectrograph after mjd 59005
           max_mjd=70000; Needs to change
         endif else begin
           min_mjd=0
           max_mjd=59030
         endelse
         if temp_mjd ge min_mjd and temp_mjd lt max_mjd then begin
         ;print,temp_mjd,nplan,temp_plan,'HOLA'
         yanny_read, allplan[iplan], pp, hdr=hdr
         thismjd = long(yanny_par(hdr, 'MJD'))
         ;print,thismjd,temp_mjd
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
         endif
         yanny_free, pp
       endfor
       ;print,mjdlist
       ;stop
       if (keyword_set(allexp)) then begin
         ;----------
         ; Determine all the plate plugging names
         allmaps = allexp.mapname
         allmaps = allmaps[ uniq(allmaps, sort(allmaps)) ]
         ;----------
         ; Loop through all plate plugging names
         ;print,allmaps
         ;stop
         ;for imap=0, n_elements(allmaps)-1 do begin
           ;indx = where(allexp.mapname EQ allmaps[imap] $
           ;  AND (allexp.flavor EQ 'science' OR allexp.flavor EQ 'smear'))
           indx = where(allexp.flavor EQ 'science' OR allexp.flavor EQ 'smear')
           if (indx[0] NE -1) then spexp = allexp[indx] $
           else spexp = 0
           ;----------
           ; Decide if any of these MJD's are within the bounds
           ; specified by MJD,MJSTART,MJEND.  If so, then set QMJD=1
           qmjd = 1B
           if (keyword_set(spexp)) then begin
             mjdlist1 = mjdlist[indx]
             qmjd = mjd_match(mjdlist1, mjd=mjd, mjstart=mjstart, mjend=mjend)
             ;if (qmjd EQ 0) then $
               ;splog, 'Skip MAP=', allmaps[imap], ' with MJD=', $
               ;mjdlist1[ uniq(mjdlist1, sort(mjdlist1)) ]
           endif
           ;print, allmaps[imap], qmjd, keyword_set(spexp)
           if (keyword_set(spexp) AND qmjd) then begin
             ;----------
             ; Determine the 2D plan files that are relevant
             planlist1 = planlist[indx]
             planlist1 = planlist1[ uniq(planlist1, sort(planlist1)) ]
             ;----------
             ; Replace the prefix 'sdR' with 'spFrame' in the science frames
             ; and the suffix '.fit' with '.fits'
             ;print,planlist1
             ;print,spexp.mjd
             ;print,max(spexp.mjd)-min(spexp.mjd)
             ;HJIM epoch definition for the plate program
             ;planlist1_org=planlist1
             RM_plates=plates_f;[15000,15001,15002,15038,15070,15071,15171,15172,15173,15252,15253]
             plt=spexp[0].plateid
             junk = where(RM_plates eq plt ,ct)
             epc=1000;this number force to coadd during all the pluggin mapname
             ;if ct gt 0 then epc=3
             spexp_org=spexp
             stop_loop=0
             ;print,spexp.mjd
             ;if max(spexp_org.mjd)-min(spexp_org.mjd) gt 3 then begin
               ;while max(spexp.mjd)-min(spexp.mjd) gt 3  do begin
               while stop_loop eq 0 do begin
                 ;minmjdepoch=min(spexp.mjd)
                 ;if max(spexp_org.mjd) eq max(spexp.mjd) then begin
                 ;  stop_loop=1
                 ;endif else begin
                 ;  stop_loop=0
                 ;endelse
                 indx_ep = where(spexp.mjd lt min(spexp.mjd)+epc)
                 ;print,max(spexp.mjd)-min(spexp.mjd)
                 spexp=spexp[indx_ep]
                 epoch_tag = replicate( $
                       {epoch_combine: min(spexp.mjd)}, n_elements(spexp))
                 spexp = struct_addtags(spexp, epoch_tag)
                 ;print,ptr_new(spexp)
                 ;print,temp.mjd
                 ;print,spexp.mjd
               ;endwhile
                 ;endif
                 ;print,max(spexp[0].mjd)
                 ;stop
                 newnames = spexp.name
                 for i=0, n_elements(newnames)-1 do begin
                   jj = strpos(newnames[i], '-')
                   kk = strpos(newnames[i], '.', /reverse_search)
                   if (jj NE -1 AND kk NE -1) then $
                     newnames[i] = 'spFrame' + strmid(newnames[i], jj, kk-jj) $
                     + '.fits'
                 endfor
                 spexp.name = newnames
                 spexp.plateid='00000';platedir
                 ;----------
                 ; Determine names of output files
                 pltid = spexp[0].plateid
                 platestr = plate_to_string(pltid)
                 thismjd = max(spexp.mjd)
                 ;print,min(spexp.mjd),max(spex.mjd)
                 mjdstr = string(thismjd, format='(i05.5)')
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
                 ;iprint,spexp_org.mjd,'HOLA1'
                 ;print,spexp.mjd,'HOLA1'
                 ;print,indx_ep
                 if (indx_ep[0] NE -1) then spexp = spexp_org[indx_ep]
                 ;print,spexp.mjd,'HOLA2'
               endwhile
             ;endif;End of the epoch definition for the plate program
           endif
         ;endfor ; End loop through plate plugging names
       endif
     ;endfor
   endif else begin
     ;---------------------------------------------------------------------------
     ; Loop through each input configuration directory
     for ifield=0, N_elements(fieldlist)-1 do begin
        fielddir = fieldlist[ifield]
        splog, ''
        splog, 'bhm field directory ', fielddir
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
           ; Determine all the fps configuration names
           ;allmaps = allexp.mapname
           ;allmaps = allmaps[ uniq(allmaps, sort(allmaps)) ]
           ;----------
           ; Loop through all fps configuration names
           ;for imap=0, n_elements(allmaps)-1 do begin
              ;indx = where(allexp.mapname EQ allmaps[imap] $
              ; AND (allexp.flavor EQ 'science' OR allexp.flavor EQ 'smear'))
              indx = where((allexp.flavor EQ 'science' OR allexp.flavor EQ 'smear'))
              if (indx[0] NE -1) then begin
                spexp = allexp[indx]
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
              ;print, allmaps[imap], qmjd, keyword_set(spexp)
              if (keyword_set(spexp) AND qmjd) then begin
                 ;spexp=spexp[qmjd]
                 ;----------
                 ; Determine the 2D plan files that are relevant
                 planlist1 = planlist[indx]
                 ;planlist1 = planlist1[qmjd]
                 ;print, planlist1u
                 ;----------
                 ; Replace the prefix 'sdR' with 'spFrame' in the science frames
                 ; and the suffix '.fit' with '.fits'
                 if keyword_set(mjd) then begin
                    nind=where(spexp.mjd EQ mjd)
                    spexp=spexp[nind]
                    planlist1f=planlist1[nind]
                    planlist1=planlist1f[0]
                 endif
                 if (keyword_set(mjstart) AND keyword_set(mjend) AND qmjd) then begin
                   nind=where((spexp.mjd LE mjend) and (spexp.mjd GE mjstart))
                   spexp=spexp[nind]
                   planlist1f=planlist1[nind]
                   planlist1=planlist1f[0]
                 endif
                 planlist1 = planlist1[ uniq(planlist1, sort(planlist1)) ]
                 ;print, planlist1
                 ;exit
                 newnames = spexp.name
                 ;print, newnames
                 for i=0, n_elements(newnames)-1 do begin
                    jj = strpos(newnames[i], '-')
                    kk = strpos(newnames[i], '.', /reverse_search)
                    if (jj NE -1 AND kk NE -1) then $
                     newnames[i] = 'spFrame' + strmid(newnames[i], jj, kk-jj) $
                      + '.fits'
                 endfor
                 spexp.name = newnames
                 ;exit
                 ;----------
                 ; Determine names of output files
                 ; HJIM -- change plate by confi
                 ;conid = spexp[0].confiid
                 ;confistr = config_to_string(conid)
                 ;confistr = spexp[0].confiid
                 fieldstr = spexp[0].fieldid
                 thismjd = max(spexp.mjd)
                 ;thismjd =spexp.mjd
                 mjdstr = string(thismjd, format='(i05.5)')
                 outdir = concat_dir(topdir, fielddir)
                 planfile = 'spPlancomb-' + fieldstr + '-' + mjdstr + '.par'
                 print, planfile
                 ;----------
                 ; Create keyword pairs for plan file
                 hdr = ''
                 ;hdr = [hdr, "confiid  " + confistr + "  # FPS Configuration number"]
                 hdr = [hdr, "fieldid  " + fieldstr + "  # BHM Field number"]
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
           ;endfor ; End loop through fps configuration names
        endif
     endfor
   endelse
   
   return
end
;------------------------------------------------------------------------------
