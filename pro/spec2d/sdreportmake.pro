;+
; NAME:
;   sdreportmake
;
; PURPOSE:
;   Generate sdReport files from the raw (idR) FITS files.
;
; CALLING SEQUENCE:
;   sdreportmake, [ speclog_dir=, mjd=, mjstart=, mjend=, /clobber ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   speclog_dir- Top directory name for output files; default to that
;                set by the environment variable SPECLOG_DIR.
;   mjd        - Look for raw data files in RAWDATA_DIR/MJD; default to
;                search all subdirectories.  Note that this need not be
;                integer-valued, but could be for example '51441_test'.
;   mjstart    - Starting MJD.
;   mjend      - Ending MJD.
;   clobber    - If set, then over-write conflicting sdReport files; default to
;                not over-write files.
;
; OUTPUT:
;
; COMMENTS:
;   The environment variables ASTROLOG_DIR and RAWDATA_DIR must be set.
;   Look for the original sdReport files and the raw FITS data files in:
;     ASTROLOG_DIR/MJD/sdReport-mmmmm.par
;     RAWDATA_DIR/MJD/sdR-cs-eeeeeeee.fit
;   where c=color, s=spectrograph number, eeeeeeee=exposure number.
;   At present, only the header keyword pairs are retained from the
;   original sdReport files.
;
;   The output sdReport files created are:
;     SPECLOG_DIR/MJD/sdReport-mmmmm.par
;   where SPECLOG_DIR is specified either by the environment variable
;   by that name, or by a keyword to this procedure.
;
; EXAMPLES:
;   Generate the sdReport files for the night of data from MJD=51441,
;   and output to the directory '/u/schlegel/test'
;   > sdreportmake, mjd=51441, speclog_dir='/u/schlegel/test'
;
; BUGS:
;   This routine spawns the Unix command 'mkdir'.
;   The use of CONCAT_DIR may not be valid for non-Unix OS's.
;
; PROCEDURES CALLED:
;   concat_dir()
;   fileandpath()
;   get_mjd_dir()
;   splog
;   sdsshead()
;   spplan_findrawdata
;   sxpar()
;   yanny_read
;   yanny_write
;
; REVISION HISTORY:
;   30-Apr-2001  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro sdreportmake, speclog_dir=speclog_dir, mjd=mjd, $
 mjstart=mjstart, mjend=mjend, minexp=minexp, clobber=clobber

   ;----------
   ; Set output structure

   outstruct1 = create_struct(name='SEXP', $
    'SEQID'     ,  0L, $
    'TILEID'    ,  0L, $
    'PLATEID'   ,  0L, $
    'MAPID'     ,  0L, $
    'FLAVOR'    ,  '                   ', $
    'QUALITY'   ,  '                   ', $
    'PROGRAM'   ,  '                   ', $
    'EXPOSURE'  ,  0L, $
    'EXPTIME'   , 0.d, $
    'MJD'       , 0.d, $
    'RA'        , 0.d, $
    'DEC'       , 0.d, $
    'CAMERAS'   ,  '                   ', $
    'TARGETNAME',  '                   ', $
    'TAIHMS'    ,  '                   ' )

   ;----------
   ; Read environment variables ASTROLOG_DIR and RAWDATA_DIR for
   ; finding raw data files.

   astrolog_dir = getenv('ASTROLOG_DIR')
   if (NOT keyword_set(astrolog_dir)) then $
    message, 'Must set environment variable ASTROLOG_DIR'

   rawdata_dir = getenv('RAWDATA_DIR')
   if (NOT keyword_set(rawdata_dir)) then $
    message, 'Must set environment variable RAWDATA_DIR'

   ;----------
   ; Determine the top-level of the output directory tree (SPECLOG_DIR)

   if (NOT keyword_set(speclog_dir)) then begin
      speclog_dir = getenv('SPECLOG_DIR')
      if (NOT keyword_set(speclog_dir)) then $
       message, 'Must set environment variable SPECLOG_DIR'
   endif
   splog, 'Setting output directory SPECLOG_DIR=' + speclog_dir

   ;----------
   ; Create a list of the MJD directories (as strings)

   mjdlist = get_mjd_dir(rawdata_dir, mjd=mjd, mjstart=mjstart, mjend=mjend)

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

         outstruct = replicate(outstruct1, nfile)

for i=0, 3 do begin ;???
;         for i=0, nfile-1 do begin
            hdr = sdsshead(fullname[i])

            if (size(hdr,/tname) EQ 'STRING') then begin
               ; Read the MJD from the header of the first file
               ; This should usually be the same as MJDDIR, though an integer
               ; rather than a string variable.
               if (i EQ 0) then mjdstr = strtrim(string(sxpar(hdr, 'MJD')),2)

               outstruct[i].seqid = 0 ; Always set to zero!!
               outstruct[i].tileid = sxpar(hdr, 'TILEID')
               outstruct[i].plateid = sxpar(hdr, 'PLATEID')
               outstruct[i].mapid = sxpar(hdr, 'MAPID')
               outstruct[i].flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)
               outstruct[i].quality = strtrim(sxpar(hdr, 'QUALITY'),2)
               outstruct[i].program = strtrim(sxpar(hdr, 'PROGRAM'),2)
               outstruct[i].exposure = sxpar(hdr, 'EXPOSURE')
               outstruct[i].exptime = sxpar(hdr, 'EXPTIME')
               tai = sxpar(hdr,'TAI-BEG')
               outstruct[i].mjd = tai / (24.D*3600.D) ; Fractional value!
               outstruct[i].ra = sxpar(hdr, 'RA')
               outstruct[i].dec = sxpar(hdr, 'DEC')
               outstruct[i].cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
               outstruct[i].targetname = strtrim(sxpar(hdr, 'NAME'),2)
               outstruct[i].taihms = strtrim(sxpar(hdr, 'TAIHMS'),2)

            endif
         endfor ; Loop over individual FITS files

         ; Read any existing sdReport file (at least the header!)
         sdreportname = string(mjdlist[imjd], $
          format='("sdReport-", i5.5, ".par")')
         oldfile = djs_filepath(sdreportname, $
          root_dir=concat_dir(astrolog_dir, mjdstr))
         yanny_read, oldfile, hdr=sdhdr, /anonymous

         ;----------
         ; Write output file

         outdir = concat_dir(speclog_dir, mjdstr)
         outfile = djs_filepath(sdreportname, root_dir=outdir)

         spawn, 'mkdir -p ' + outdir
         qexist = keyword_set(findfile(outfile))
         if (qexist) then begin
            if (keyword_set(clobber)) then $
             splog, 'WARNING: Over-writing sdReport file: ' + outfile $
            else $
             splog, 'WARNING: Will not over-write sdReport file: ' + outfile
         endif
         if ((NOT qexist) OR keyword_set(clobber)) then begin
            splog, 'Writing sdReport file ', outfile
            yanny_write, outfile, ptr_new(outstruct), hdr=sdhdr
         endif

      endif
   endfor ; End loop over MJD directories

   return
end
;------------------------------------------------------------------------------
