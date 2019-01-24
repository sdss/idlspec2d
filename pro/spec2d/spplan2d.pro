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
;   files did not properly describe the obsSummary name.  In those cases,
;   look for the actual obsSummary files in SDSSCORE/MJD.
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
pro spplan2d, topdir=topdir1, run2d=run2d1, mjd=mjd, lco=lco $
 mjstart=mjstart, mjend=mjend, minexp=minexp, clobber=clobber, dr13=dr13, $
 _extra=foo

   if (NOT keyword_set(minexp)) then minexp = 1
   if keyword_set(lco) then begin
     obsdir='LCO'
   endif else begin
     obsdir='APO'
   endelse
   ;----------
   ; Determine the top-level of the output directory tree
   if (keyword_set(topdir1)) then topdir = topdir1 $
    else topdir = getenv('BOSS_SPECTRO_REDUX')
   topdir=concat_dir(topdir, obsdir)
   splog, 'Setting TOPDIR=', topdir
   if (keyword_set(run2d1)) then run2d = strtrim(run2d1,2) $
    else run2d = getenv('RUN2D')
   splog, 'Setting RUN2D=', run2d

   ;----------
   ; Read environment variable for BOSS_SPECTRO_DATA for finding raw data files.

   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   if (NOT keyword_set(rawdata_dir)) then $
    message, 'Must set environment variable BOSS_SPECTRO_DATA'
   rawdata_dir = concat_dir(rawdata_dir, obsdir)
   splog, 'Setting BOSS_SPECTRO_DATA=', rawdata_dir

   ;speclog_dir = getenv('SPECLOG_DIR')
   sdsscore_dir = getenv('SDSSCORE')
   if (NOT keyword_set(sdsscore_dir)) then $
    message, 'Must set environment variable SDSSCORE'
   sdsscore_dir  = concat_dir(sdsscore_dir, obsdir)
   splog, 'Setting SDSSCORE=', sdsscore_dir

   spawn, 'speclog_version', logvers, /noshell

   ;----------
   ; Create a list of the MJD directories (as strings)

   mjdlist = get_mjd_dir(rawdata_dir, mjd=mjd, mjstart=mjstart, mjend=mjend)
   nmjd = n_elements(mjdlist)
   splog, 'Number of MJDs = ', nmjd
   ;;HJIM -- reduce the number of spectrographs to one
   camnames = ['b1', 'r1']
   ncam = N_elements(camnames)

   ;---------------------------------------------------------------------------
   ; Loop through each input MJD directory

   for imjd=0, nmjd-1 do begin

      mjddir = mjdlist[imjd]
      thismjd = long(mjdlist[imjd])
      inputdir = concat_dir(rawdata_dir, mjddir)
      confdir = concat_dir(sdsscore_dir, mjddir);HJIM Needs to check the final path for the obsSummary file 
      plugdir='a'
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
         ; HJIM-- Change FIBERID by CONFIID
         CONFIID = strarr(nfile)
         EXPTIME = fltarr(nfile)
         EXPOSURE = lonarr(nfile)
         FLAVOR = strarr(nfile)
         CAMERAS = strarr(nfile)
         MAPNAME = strarr(nfile)
         FIELDID = strarr(nfile)
         CONFNAME = strarr(nfile) ; Added by HJIM 

         for i=0, nfile-1 do begin
            hdr = sdsshead(fullname[i])

            if (size(hdr,/tname) EQ 'STRING') then begin

               CONFIID[i] = strtrim( sxpar(hdr, 'CONFIID') );change long plate  format to string format 
               ;CONFIID[i] = long( sxpar(hdr, 'PLATEID') )
               EXPTIME[i] = sxpar(hdr, 'EXPTIME')
               EXPOSURE[i] = long( sxpar(hdr, 'EXPOSURE') )
               FLAVOR[i] = strtrim(sxpar(hdr, 'FLAVOR'),2)
               CAMERAS[i] = strtrim(sxpar(hdr, 'CAMERAS'),2)
               MAPNAME[i] = strtrim(sxpar(hdr, 'NAME'),2)
               map_name=strsplit(MAPNAME[i],'-',/extract)
               CONFNAME[i] = map_name[0]
               ;; Check CONFTYP for BOSS or EBOSS (e.g. not MANGA)
               ;; If keyword is missing (older data), assume this is BOSS
               platetype = sxpar(hdr, 'CONFTYP', count=nhdr)
               if (nhdr GT 0) then begin
                   platetype = strupcase(strtrim(platetype,2))
                   if (platetype NE 'BOSS') && (platetype NE 'EBOSS') then begin

                       splog, 'Skipping ' + platetype + $
                           ' configuration '+ CONFIID[i] + $
                           ' exposure ', EXPOSURE[i]
                       FLAVOR[i] = 'unknown'
                   endif ;else begin 
                    ;; Skip also eBOSS plates and some RM plates for DR13
                    ;if keyword_set(dr13) then begin
                       ;if (platetype NE 'BOSS') OR $
                       ;( (PLATEID[i] EQ 7338 OR PLATEID[i] EQ 7339 OR PLATEID[i] EQ 7340) AND thismjd GT 57000) $
                       ;then begin
                       ;    splog, 'Skipping ' + platetype + $
                       ;    ' plate ', PLATEID[i], $
                       ;    ' exposure ', EXPOSURE[i], ' for DR13', thismjd
                       ;FLAVOR[i] = 'unknown'
                       ;endif
                    ;endif 
                   ;endelse 
               endif

               ;-- Removing exposure 258988 of plate 9438 mjd 58125 because of trail in data
               ;if sxpar( hdr, 'EXPOSURE') EQ 258988L then FLAVOR[i] = 'unknown'

               ; Exclude all files where the QUALITY keyword is not 'excellent'.
               quality = strtrim(sxpar(hdr, 'QUALITY'),2)
               if (quality NE 'excellent') then begin
                  splog, 'Warning: Non-excellent quality '+FLAVOR[i]+' file ' $
                   + fileandpath(fullname[i]) + ' ('+quality+')'
                  FLAVOR[i] = 'unknown'
               endif

               ; Exclude files where the plate number does not match that
               ; in the map name
				; JEB -- plate number
				; HJIM -- configuration number
				
				       ;map_name=strsplit(MAPNAME[i],'-',/extract)
               ;if (CONFIID[i] NE (map_name[0] + '-' + map_name[1])) $
               if (CONFIID[i] NE (map_name[0])) $
                && (FLAVOR[i] NE 'bias') then begin
                  platestr = strtrim(CONFIID[i])
                  splog, 'Warning: Configuration number ' + platestr $
                   + ' flavor '+ FLAVOR[i] $
                   + ' inconsistent with map name ' + fileandpath(fullname[i])
                  FLAVOR[i] = 'unknown'
               endif

               if (sxpar(hdr, 'MJD') NE thismjd) then $
                splog, 'Warning: Wrong MJD in file '+fileandpath(fullname[i])

               ; MAPNAME should be of the form '000000-51683-01'.
               ; If it only contains the CONFIID ;; (for MJD <= 51454),
               ; then find the actual plug-map file.
               if (strlen(MAPNAME[i]) LE 15) then begin
                  confile = 'obsSummary-' $
                   ;+ string(long(MAPNAME[i]), format='(i4.4)') + '-*.par'
                   + map_name[0] + '-' + map_name[1] + '-*.par'
                  confile = (findfile(filepath(confile, root_dir=confdir), $
                   count=ct))[0]
                  if (ct EQ 1) then $
                   MAPNAME[i] = strmid(fileandpath(confile), 11, 15)
                   confile = 'obsSummary-'+MAPNAME[i]+'.par'
                   thisplan=filepath(confile, root_dir=confdir)
                   allseq = yanny_readone(thisplan, 'SPEXP', hdr=hdr1, /anon)
                   thisfield=strtrim(string(yanny_par(hdr1,'bhmfield_id')),2)
                   FIELDID[i]=thisfield;field_id
               endif
            endif
         endfor

         ;----------
         ;
         ; Determine all the conofiguration names
         ;allmaps = MAPNAME[ uniq(MAPNAME, sort(MAPNAME))]
         allfield = FIELDID[ uniq(FIELDID, sort(FIELDID))]
         allconfs = CONFNAME[ uniq(FIELDID, sort(FIELDID))]
         ;----------
         ; Loop through all configuration pointing names

         ;for imap=0, n_elements(allconfs)-1 do begin
         for imap=0, n_elements(allfield)-1 do begin 
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

               if (ct GT 0) then begin
                  spexp1 = spplan_create_spexp(allexpnum[iexp], $
                  CONFIID[indx[0]], thismjd, FIELDID[indx[0]], $
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
              conid =config_to_long(spexp[0].confiid)
              fieid =config_to_long(spexp[0].fieldid)
              ;if (pltid GT 0 AND pltid LT 9990) then begin
              ;   platestr = string(pltid, format='(i04.4)')
              if (conid GE 0) then begin
                  if (fieid GE 0) then begin
                      fieldstr = string(fieid, format='(i04.4)')
                  endif else begin
                       splog, 'WARNING: Field number '+strtrim(string(fieid),2)+' invalid for COFNAME=' + allconfs[imap]
                       fieldstr = '0000'
                  endelse
                  ;confistr = spexp[0].confiid; 
                  confistr = config_to_string(conid) 
                  ;print, confistr
              endif else begin
                 if (fieid GE 0) then begin
                    fieldstr = string(fieid, format='(i04.4)')
                 endif else begin
                    splog, 'WARNING: Field number '+strtrim(string(fieid),2)+' invalid for COFNAME=' + allconfs[imap]
                     fieldstr = '0000'
                 endelse
                 splog, 'WARNING: Configuration number '+strtrim(string(conid),2)+' invalid for COFNAME=' + allconfs[imap]
                 confistr = '000000'
                 spexp = 0
              endelse
            endif else begin
              splog, "WARNING: no good exposures for CONFNAME="+allconfs[imap]
              conid = 0
              fieid = 0
            endelse
             
            mjdstr = string(thismjd, format='(i05.5)')

            ;----------
            ; Discard these observations if there is not at least one flat,
            ; one arc, and one science exposure

            if (keyword_set(spexp)) then begin
               junk = where(spexp.flavor EQ 'flat', ct)
               if (ct EQ 0) then begin
               ;   if (flatt EQ 0) then begin
                     splog, 'WARNING: No flats for CONFNAME=' + allconfs[imap]
                     spexp = 0
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

         endfor ; End loop through plate plugging names
      endif
   endfor

   return
end
;------------------------------------------------------------------------------
