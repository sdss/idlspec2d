;+
; NAME:
;   spplan2d
;
; PURPOSE:
;   Create plan file(s) for running the Spectro-2D pipeline.
;
; CALLING SEQUENCE:
;   spplan2d, [ topdir=, run2d=, mjd=, mjstart=, mjend=, minexp=, $
;    /clobber, /dr13 ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   topdir     - Optional override value for the environment
;                variable $BOSS_SPECTRO_REDUX.
;   run2d      - Optional override value for the environment variable $RUN2D
;   mjd        - Look for raw data files in BOSS_SPECTRO_DATA/MJD; default to
;                search all subdirectories.  Note that this need not be
;                integer-valued, but could be for example '51441_test'.
;   mjstart    - Starting MJD.
;   mjend      - Ending MJD.
;   minexp     - Minimum exposure time for science frames; default to 1 sec.
;   clobber    - If set, then over-write conflicting plan files; default to
;                not over-write files.
;
; OUTPUT:
;
; COMMENTS:
;   The environment variables SDSSCORE and BOSS_SPECTRO_DATA must be set.
;
;   Look for the raw FITS data files in:
;     BOSS_SPECTRO_DATA/MJD/sdR-cs-eeeeeeee.fit
;   where c=color, s=spectrograph number, eeeeeeee=exposure number.
;
;   The output 2D plan files created are:
;     TOPOUTDIR/PLATE/spPlan2d-cccc-mmmmm.par
;   where cccc=configuration number, mmmmm=MJD.
;
;   For the earliest data (before MJD=51455), then NAME keyword in the FITS
;   files did not properly describe the confsummary name.  In those cases,
;   look for the actual confsummary files in SDSSCORE.
;
;   Exclude all files where the QUALITY header keyword is not 'excellent'.
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
;   idlutils_version()
;   splog
;   sdsshead()
;   spplan_findrawdata
;   spplan_create_spexp
;   sxpar()
;   yanny_write
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;   05-Sep-2015 Plate number armaggedon correction by JEB.
;   15-Nov-2018  Modified by Hector Ibarra for the BHM
;-
;------------------------------------------------------------------------------
pro spplan2d_new, topdir=topdir1, run2d=run2d1, mjd=mjd, lco=lco, $
 mjstart=mjstart, mjend=mjend, minexp=minexp, clobber=clobber, dr13=dr13, $
 _extra=foo, legacy=legacy, plates=plates, nocomm=nocomm, $
 twilight_flats=twilight_flats, test_twilight=test_twilight

 RESOLVE_ALL, /QUIET, /SKIP_EXISTING, /CONTINUE_ON_ERROR

   if (NOT keyword_set(minexp)) then minexp = 1
   if keyword_set(lco) then begin
     obsdir='LCO'
     sdsscore_obs='lco'
     BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_S'
   endif else begin
     obsdir='APO'
     sdsscore_obs='apo'
     BOSS_SPECTRO_DATA='BOSS_SPECTRO_DATA_N'
   endelse
   obsdir='';coment this line for the final version HJIM
   ;----------
   ; Determine the top-level of the output directory tree
   if (keyword_set(topdir1)) then topdir = topdir1 $
    else topdir = getenv('BOSS_SPECTRO_REDUX')
   ;topdir=concat_dir(topdir, obsdir)
   splog, 'Setting TOPDIR=', topdir
   if (keyword_set(run2d1)) then run2d = strtrim(run2d1,2) $
    else run2d = getenv('RUN2D')
   splog, 'Setting RUN2D=', run2d

   ;----------
   ; Read environment variable for BOSS_SPECTRO_DATA for finding raw data files.

   rawdata_dir = getenv(BOSS_SPECTRO_DATA)
   if (NOT keyword_set(rawdata_dir)) then $
    message, 'Must set environment variable BOSS_SPECTRO_DATA'
   ;rawdata_dir = concat_dir(rawdata_dir, obsdir)
   splog, 'Setting BOSS_SPECTRO_DATA=', rawdata_dir
   
   if keyword_set(legacy) or keyword_set(plates) then begin
      speclog_dir = getenv('SPECLOG_DIR')
      if (NOT keyword_set(speclog_dir)) then $
        message, 'Must set environment variable SPECLOG_DIR'
      splog, 'Setting SPECLOG_DIR=', speclog_dir
   endif else begin
      sdsscore_dir = getenv('SDSSCORE_DIR')
      if (NOT keyword_set(sdsscore_dir)) then $
        message, 'Must set environment variable SDSSCORE_DIR'
      sdsscore_dir  = concat_dir(sdsscore_dir, sdsscore_obs)
      sdsscore_dir  = concat_dir(sdsscore_dir, 'summary_files')
      splog, 'Setting SDSSCORE_DIR=', sdsscore_dir
   endelse
   spawn, 'speclog_version', logvers, /noshell

   ;----------
   ; Create a list of the MJD directories (as strings)

   mjdlist = get_mjd_dir(rawdata_dir, mjd=mjd, mjstart=mjstart, mjend=mjend)
   nmjd = n_elements(mjdlist)
   splog, 'Number of MJDs = ', nmjd
   ;;HJIM -- reduce the number of spectrographs to one
   camnames = ['b1', 'r1']
      plateflavor0='BHM'
      plateflavor1='BHM&MWM'
   if keyword_set(plates) then begin
      plateflavor0='BHM';'BOSSHALF'
      plateflavor1='BHM&MWM';'APOGEE-BOSS'
   endif
   if keyword_set(legacy) then begin
      camnames = ['b1', 'r1', 'b2', 'r2']
      plateflavor0='EBOSS'
      plateflavor1='BOSS'
   endif
   ncam = N_elements(camnames)

   ;---------------------------------------------------------------------------
   ; Loop through each input MJD directory

   for imjd=0, nmjd-1 do begin

      mjddir = mjdlist[imjd]
      thismjd = long(mjdlist[imjd])
      inputdir = concat_dir(rawdata_dir, mjddir)
      if keyword_set(legacy) or keyword_set(plates) then begin 
         plugdir = concat_dir(speclog_dir, mjddir)
      endif else begin
         confdir = sdsscore_dir; concat_dir(sdsscore_dir, mjddir);HJIM Needs to check the final path for the confsummary file
      endelse
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
         ; HJIM-- Change FIBERID by CONFID
         
         EXPTIME = fltarr(nfile)
         EXPOSURE = lonarr(nfile)
         FLAVOR = strarr(nfile)
         CAMERAS = strarr(nfile)
         MAPNAME = strarr(nfile)
         TAI = fltarr(nfile)
         if keyword_set(legacy) or keyword_set(plates) then begin
           PLATEID = lonarr(nfile)
         endif else begin
           FIELDID = strarr(nfile) ; Added by HJIM
           CONFID = strarr(nfile) ; Added by HJIM
           CONFNAME = strarr(nfile) ; Added by HJIM          
         endelse


         for i=0, nfile-1 do begin
            if (i eq 0) then begin
		hdr = sdsshead(fullname[i])
	    endif else begin
		 hdr = sdsshead(fullname[i],/silentwarn)
            endelse
            if (size(hdr,/tname) EQ 'STRING') then begin

               EXPTIME[i] = sxpar(hdr, 'EXPTIME')
               EXPOSURE[i] = long( sxpar(hdr, 'EXPOSURE') )
               FLAVOR[i] = strtrim(sxpar(hdr, 'FLAVOR'),2)
               CAMERAS[i] = strtrim(sxpar(hdr, 'CAMERAS'),2)
               TAI[i] = sxpar(hdr, 'TAI-BEG')
               if keyword_set(legacy) or keyword_set(plates) then begin
                 MAPNAME[i] = strtrim(sxpar(hdr, 'NAME'),2)
                 PLATEID[i] = long( sxpar(hdr, 'PLATEID') )
                 platetype = sxpar(hdr, 'PLATETYP', count=nhdr)
               endif else begin
                 MAPNAME[i] = strtrim(sxpar(hdr, 'CONFID'),2)
                 map_name = strarr(1)
                 map_name[0]=MAPNAME[i] ; strsplit(MAPNAME[i],'-',/extract)
                 CONFNAME[i] = map_name[0]
                 CONFID[i] = strtrim( sxpar(hdr, 'CONFID') );change long plate  format to string format
                 platetype = 'BHM&MWM' ;sxpar(hdr, 'CONFTYP', count=nhdr)
                 nhdr = 1
               endelse
               ;; Check CONFTYP for BOSS or EBOSS (e.g. not MANGA)
               ;; If keyword is missing (older data), assume this is BOSS
               if (nhdr GT 0) then begin
                   platetype = strupcase(strtrim(platetype,2))
                   if (platetype NE plateflavor0) && (platetype NE plateflavor1) then begin
                    if keyword_set(legacy) or keyword_set(plates) then begin
                       splog, 'Skipping ' + platetype + $
                           ' plate ', PLATEID[i], $
                           ' exposure ', EXPOSURE[i]
                       FLAVOR[i] = 'unknown'
                    endif else begin
                       splog, 'Skipping ' + platetype + $
                           ' configuration '+ CONFID[i] + $
                           ' exposure ', EXPOSURE[i]
                       FLAVOR[i] = 'unknown'
                    endelse
                   endif else begin
                    if keyword_set(legacy) then begin                  
                      ;; Skip also eBOSS plates and some RM plates for DR13
                      if keyword_set(dr13) then begin
                        if (platetype NE 'BOSS') OR $
                        ( (PLATEID[i] EQ 7338 OR PLATEID[i] EQ 7339 OR PLATEID[i] EQ 7340) AND thismjd GT 57000) $
                        then begin
                           splog, 'Skipping ' + platetype + $
                           ' plate ', PLATEID[i], $
                           ' exposure ', EXPOSURE[i], ' for DR13', thismjd
                        FLAVOR[i] = 'unknown'
                        endif
                      endif 
                    endif 
                   endelse 
               endif
               if keyword_set(legacy) then begin
                  ;-- Removing exposure 258988 of plate 9438 mjd 58125 because of trail in data
                  if sxpar( hdr, 'EXPOSURE') EQ 258988L then FLAVOR[i] = 'unknown'
               endif
               ; Exclude all files where the QUALITY keyword is not 'excellent'.
               quality = strtrim(sxpar(hdr, 'QUALITY'),2)
               if (quality NE 'excellent') then begin
                  splog, 'Warning: Non-excellent quality '+FLAVOR[i]+' file ' $
                   + fileandpath(fullname[i]) + ' ('+quality+')'
                  FLAVOR[i] = 'unknown'
               endif
               ; Exclude files where the plate or configuration number does not match that
               ; in the map name
               if keyword_set(legacy) or keyword_set(plates) then begin
                    ; JEB -- plate number
                    if (plate_to_string(PLATEID[i]) NE strmid(MAPNAME[i],0,strpos(MAPNAME[i],'-')))  $
                        && (FLAVOR[i] NE 'bias') then begin
                        platestr = strtrim(string(PLATEID[i]), 2)
                        splog, 'Warning: Plate number ' + platestr $
				             + ' flavor '+ FLAVOR[i] $
				             + ' inconsistent with map name ' + fileandpath(fullname[i])
                        FLAVOR[i] = 'unknown'
                    endif
               endif else begin
                    ; HJIM -- configuration number
                    if (strtrim(CONFID[i],2) NE strtrim((map_name[0]),2)) $
                        && (FLAVOR[i] NE 'bias') then begin
                        platestr = strtrim(CONFID[i])

                        splog, 'Warning: Configuration number ' + platestr $
                            + ' flavor '+ FLAVOR[i] $
                            + ' inconsistent with map name ' + fileandpath(fullname[i])
                        FLAVOR[i] = 'unknown'
                    endif
               endelse
               if (sxpar(hdr, 'MJD') NE thismjd) then $
                splog, 'Warning: Wrong MJD in file '+fileandpath(fullname[i])
               if keyword_set(legacy) or keyword_set(plates) then begin
                 ; MAPNAME should be of the form '000000-51683-01'.
                 ; If it only contains the PLATEID ;; (for MJD <= 51454),
                 ; then find the actual plug-map file.
                 if (strlen(MAPNAME[i]) LE 4) then begin
                   plugfile = 'plPlugMapM-' $
                     + string(long(MAPNAME[i]), format='(i4.4)') + '-*.par'
                   plugfile = (findfile(filepath(plugfile, root_dir=plugdir), $
                     count=ct))[0]
                   if (ct EQ 1) then $
                     MAPNAME[i] = strmid(fileandpath(plugfile), 11, 13)
                 endif
               endif else begin
                 if (strlen(MAPNAME[i]) LE 15) then begin
                    confile = 'confSummaryF-' + map_name[0] + '-*.par'
                    confile = (findfile(filepath(confile, root_dir=confdir, subdir='*'), count=ct))[0]
                    if (ct ne 0) then begin
                       MAPNAME[i] = strmid(fileandpath(confile), 11, 15)
                       confile = 'confSummaryF-'+ map_name[0]+'.par'
                    endif else begin
                      confile = 'confSummary-' + map_name[0] + '-*.par'
                      confile = (findfile(filepath(confile, root_dir=confdir, subdir='*'), count=ct))[0]
                      if (ct EQ 1) then MAPNAME[i] = strmid(fileandpath(confile), 11, 15)
                      confile = 'confSummary-'+ map_name[0]+'.par'
                    endelse
                    thisplan=(findfile(filepath(confile, root_dir=confdir,subdir='*')))[0]
                    allseq = yanny_readone(thisplan, 'SPEXP', hdr=hdr1, /anon)
                    thisfield=field_to_string(yanny_par(hdr1,'field_id'))
                    if yanny_par(hdr1,'field_id')  eq -999 then thisfield=field_to_string(0)
                    if strlen(thisfield) eq 0 then thisfield=field_to_string(0)
                    FIELDID[i]=thisfield
                 endif
               endelse
            endif
         endfor

         ;----------
         ;
         if keyword_set(legacy) or keyword_set(plates) then begin
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
                 AND EXPOSURE EQ allexpnum[iexp] $
                 AND FLAVOR NE 'unknown', ct)

               if (ct GT 0) then begin
                 spexp1 = spplan_create_spexp_legacy(allexpnum[iexp], $
                   PLATEID[indx[0]], thismjd, $
                   MAPNAME[indx[0]], FLAVOR[indx[0]], EXPTIME[indx[0]], $
                   shortname[indx], CAMERAS[indx], minexp=minexp, plates=plates)
                 if (keyword_set(spexp1)) then begin
                   if (keyword_set(spexp)) then spexp = [spexp, spexp1] $
                   else spexp = spexp1
                 endif
               endif
             endfor
             ;----------
             ; Discard these observations if the plate number is not
             ; in the range 1 to 9990.
             if (keyword_set(spexp)) then begin
               pltid = long(spexp[0].plateid)
               ;if (pltid GT 0 AND pltid LT 9990) then begin
               ;   platestr = string(pltid, format='(i04.4)')
               if (pltid GT 0) then begin
                 platestr = plate_to_string(pltid)
               endif else begin
                 splog, 'WARNING: Plate number '+strtrim(string(pltid),2)+' invalid for MAPNAME=' + allmaps[imap]
                 platestr = '0000'
                 spexp = 0
               endelse
             endif else begin
               splog, "WARNING: no good exposures for MAPNAME="+allmaps[imap]
               pltid = 0
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
               junk = where(spexp.flavor EQ 'science' $
                 OR spexp.flavor EQ 'smear', ct)
               if (ct EQ 0) then begin
                 splog, 'WARNING: No science frames for MAPNAME=' + allmaps[imap]
                 spexp = 0
               endif
             endif
             if (keyword_set(spexp)) then begin
               ;----------
               ; Determine names of output files
               outdir = djs_filepath('', root_dir=topdir, $
                 subdir=[run2d,platestr]);HJIM add p to identify plate diretories
               planfile = 'spPlan2d-' + platestr + '-' + mjdstr + '.par'
               logfile = 'spDiag2d-' + platestr + '-' + mjdstr + '.log'
               plotfile = 'spDiag2d-' + platestr + '-' + mjdstr + '.ps'
               ;----------
               ; Create keyword pairs for plan file
               hdr = ''
               hdr = [hdr, "plateid  " + platestr + "  # Plate number"]
               hdr = [hdr, "MJD     " + mjdstr $
                 + "  # Modified Julian Date"]
               hdr = [hdr, "RUN2D  " + run2d + "  # 2D reduction name"]
               hdr = [hdr, "planfile2d  '" + planfile $
                 + "'  # Plan file for 2D spectral reductions (this file)"]
               hdr = [hdr, "idlspec2dVersion '" + idlspec2d_version() $
                 + "'  # Version of idlspec2d when building plan file"]
               hdr = [hdr, "idlutilsVersion '" + idlutils_version() $
                 + "'  # Version of idlutils when building plan file"]
               hdr = [hdr, "speclogVersion '" + logvers $
                 + "'  # Version of speclog when building plan file"]
               ;----------
               ; Write output file
               ; Create output directory if it does not yet exist
               if (file_test(outdir, /directory) EQ 0) then begin
                 ; Bad if this exists as a file
                 if (file_test(outdir)) then $
                   message, 'Expecting directory not file '+outdir
                 spawn, 'mkdir -p ' + outdir
               endif
               ; Bad if this exists as an unwriteable dir
               if (file_test(outdir, /directory, /write) EQ 0) then $
                 message, 'Cannot write to directory '+outdir
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
         endif else begin
           ; Determine all the conofiguration names
           allfield = FIELDID[ uniq(FIELDID, sort(FIELDID))]
           allconfs = CONFNAME[ uniq(FIELDID, sort(FIELDID))]
           ;----------
           ; Loop through all configuration pointing names
           ;for imap=0, n_elements(allconfs)-1 do begin
           for imap=0, n_elements(allfield)-1 do begin
              if keyword_set(nocomm) then begin
                if (long(allfield[imap]) gt 16000 and long(allfield[imap]) lt 100000) then continue
              endif
              spexp = 0 ; Zero-out this output structure
              ;----------
              ; Loop through all exposure numbers for this configuration pointing
              ;theseexp = EXPOSURE[ where(CONFNAME EQ allconfs[imap]) ]
              theseexp = EXPOSURE[ where(FIELDID EQ allfield[imap]) ]
              allexpnum = theseexp[ uniq(theseexp, sort(theseexp)) ]
              for iexp=0, n_elements(allexpnum)-1 do begin
                 ;indx = where(CONFNAME EQ allconfs[imap] $
                 indx = where(FIELDID EQ allfield[imap] $
                  AND EXPOSURE EQ allexpnum[iexp] $
                  AND FLAVOR NE 'unknown', ct)
                 if keyword_set(test_twilight) then $
                   indx = where(FIELDID EQ allfield[imap] $
                           AND EXPOSURE EQ allexpnum[iexp] $
                           AND FLAVOR NE 'unknown' $
                           AND FLAVOR NE 'flat', ct)
                 if (ct GT 0) then begin
                    spexp1 = spplan_create_spexp(allexpnum[iexp], $
                    CONFID[indx[0]], thismjd, FIELDID[indx[0]], $
                    MAPNAME[indx[0]], FLAVOR[indx[0]], EXPTIME[indx[0]], $
                    shortname[indx], CAMERAS[indx], minexp=minexp)
                    if (keyword_set(spexp1)) then begin
                       if (keyword_set(spexp)) then spexp = [spexp, spexp1] $
                       else spexp = spexp1
                    endif
                 endif
              endfor
              ;----------
              ; Discard these observations if the plate number is not
              ; in the range 1 to 9990.
              ; HJIM -- change plate by configuration
              if (keyword_set(spexp)) then begin
                conid =config_to_long(spexp[0].confid)
                fieid =config_to_long(spexp[0].fieldid)
                ;if (pltid GT 0 AND pltid LT 9990) then begin
                ;   platestr = string(pltid, format='(i04.4)')
                if (conid GE 0) then begin
                    if (fieid GE 0) then begin
                        fieldstr = field_to_string(fieid);Modified this to add the number of digits for the FieldID
                    endif else begin
                         splog, 'WARNING: Field number '+strtrim(string(fieid),2)+' invalid for COFNAME=' + allconfs[imap]
                         fieldstr = field_to_string(0)
                    endelse
                    ;confistr = spexp[0].confid;
                    confistr = config_to_string(conid) 
                    ;print, confistr
                endif else begin
                   if (fieid GE 0) then begin
                      fieldstr = field_to_string(fieid);Modified this to add the number of digits for the FieldID
                   endif else begin
                      splog, 'WARNING: Field number '+strtrim(string(fieid),2)+' invalid for COFNAME=' + allconfs[imap]
                       fieldstr = field_to_string(0);Modified this to add the number of digits for the FieldID
                   endelse
                   splog, 'WARNING: Configuration number '+strtrim(string(conid),2)+' invalid for COFNAME=' + allconfs[imap]
                   confistr = config_to_string(0)
                   spexp = 0
                endelse
              endif else begin
                splog, "WARNING: no good exposures for CONFNAME="+allconfs[imap]
                conid = 0
                fieid = 0
              endelse
  
              mjdstr = string(thismjd, format='(i05.5)')
              if keyword_set(twilight_flats) then begin
                if (keyword_set(spexp)) then begin
                    junk = where(spexp.flavor EQ 'science', cts)
                    if cts ne 0 then begin
                        junk = where(spexp.flavor EQ 'flat', ct)
                        if (ct EQ 0) then begin
                            f_indx = where(FLAVOR EQ 'flat')
                            f_tai = TAI[f_indx]
                            s_indx = where(shortname EQ (spexp.name[0])[0])
print, abs(f_tai - TAI[s_indx])
                            d_tai = min(abs(f_tai - TAI[s_indx]), match)
                            f_use=f_indx[match[0]]
                            f_use=where(EXPOSURE EQ EXPOSURE[f_use])
                            spexp1 = spplan_create_spexp(EXPOSURE[f_use[0]], CONFID[f_use[0]],
                                                         thismjd, (spexp.fieldid)[0], MAPNAME[f_use[0]],$
                                                         FLAVOR[f_use[0]], EXPTIME[f_use[0]], $
                                                        shortname[f_use], CAMERAS[f_use], minexp=minexp)
                            if (keyword_set(spexp1)) then begin
                                if (keyword_set(spexp)) then spexp = [spexp, spexp1] $
                            else spexp = spexp1
                        endif
                    endif
                endif
              endif
              ;----------
              ; Discard these observations if there is not at least one flat,
              ; one arc, and one science exposure
              if (keyword_set(spexp)) then begin
                junk = where(spexp.flavor EQ 'flat', ct)
                if (ct EQ 0) then begin
                ;   if (flatt EQ 0) then begin
                      splog, 'WARNING: No flats for CONFNAME=' + allconfs[imap]
                      spexp = 0
                ;   endif
                endif
                ;endif else begin
                ;   flatt =1
                ;endelse
              endif

              if (keyword_set(spexp)) then begin
                 junk = where(spexp.flavor EQ 'arc', ct)
                 if (ct EQ 0) then begin
                    splog, 'WARNING: No arcs for CONFNAME=' + allconfs[imap]
                    spexp = 0
                 endif
              endif

              if (keyword_set(spexp)) then begin
                 junk = where(spexp.flavor EQ 'science' $
                  OR spexp.flavor EQ 'smear', ct)
                 if (ct EQ 0) then begin
                    splog, 'WARNING: No science frames for CONFNAME=' + allconfs[imap]
                    spexp = 0
                 endif
              endif

              if (keyword_set(spexp)) then begin
                 ;----------
                 ; Determine names of output files
                 ; HJIM -- change fiberid by field_id
                 outdir = djs_filepath('', root_dir=topdir, $
                  subdir=[run2d,fieldstr])
                 planfile = 'spPlan2d-' + fieldstr + '-' + mjdstr + '.par'
                 logfile = 'spDiag2d-' + fieldstr + '-' + mjdstr + '.log'
                 plotfile = 'spDiag2d-' + fieldstr + '-' + mjdstr + '.ps'
                 ;----------
                 ; Create keyword pairs for plan file
                 hdr = ''
                 hdr = [hdr, "confname  " + confistr + "  # FPS configuration number"]
                 hdr = [hdr, "fieldname  " + fieldstr + "  # BHM field number"]
                 hdr = [hdr, "MJD     " + mjdstr $
                  + "  # Modified Julian Date"]
                 hdr = [hdr, "RUN2D  " + run2d + "  # 2D reduction name"]
                 hdr = [hdr, "planfile2d  '" + planfile $
                  + "'  # Plan file for 2D spectral reductions (this file)"]
                 hdr = [hdr, "idlspec2dVersion '" + idlspec2d_version() $
                  + "'  # Version of idlspec2d when building plan file"]
                 hdr = [hdr, "idlutilsVersion '" + idlutils_version() $
                  + "'  # Version of idlutils when building plan file"]
                 hdr = [hdr, "speclogVersion '" + logvers $
                  + "'  # Version of speclog when building plan file"]
                 ;----------
                 ; Write output file
                 ; Create output directory if it does not yet exist
                 if (file_test(outdir, /directory) EQ 0) then begin
                    ; Bad if this exists as a file
                    if (file_test(outdir)) then $
                      message, 'Expecting directory not file '+outdir
                    spawn, 'mkdir -p ' + outdir
                 endif
                 ; Bad if this exists as an unwriteable dir
                 if (file_test(outdir, /directory, /write) EQ 0) then $
                  message, 'Cannot write to directory '+outdir
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
           endfor ; End loop through configuration names
         endelse
      endif
   endfor
   return
end
;------------------------------------------------------------------------------
