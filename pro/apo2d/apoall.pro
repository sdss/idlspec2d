;+
; NAME:
;   apoall
;
; PURPOSE:
;   Run APOREDUCE on one or many nights of data.
;
; CALLING SEQUENCE:
;   apoall, [ mjd=, mjstart=, mjend=, minexp=, copydir= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - Look for raw data files in RAWDIR/MJD; default to '*' to
;                search all subdirectories.  Note that this need not be
;                integer-valued, but could be for example '51441_test'.
;   mjstart    - Starting MJD.
;   mjend      - Ending MJD.
;   minexp     - Minimum exposure time for science frames; default to 0 sec.
;   copydir    - Copy the output log files to this directory; default to none.
;
; OUTPUT:
;
; COMMENTS:
;   The files are sorted before being sent to APOREDUCE.  For each plate,
;   reduce all the biases/darks, then all the flats, then all the arcs,
;   and finally all of the science/smear frames.
;
;   Look for the raw sdR files in $BOSS_SPECTRO_DATA/$MJD, the plPlugMapM files
;   in $SPECLOG_DIR/$MJD, and put the outputs in $SPECTROLOG_DIR/$MJD.
;   If that last environment variable is not set, then put the outputs
;   in the same directory as the sdR files.
;
; EXAMPLES:
;   Rerun the SOS code on the spectro data from MJD 53682, putting
;   the results in the current directory:
;     IDL> setenv, 'BOSS_SPECTRO_DATA=.'
;     IDL> apoall, mjd=53682
;
; BUGS:
;
; PROCEDURES CALLED:
;   aporeduce
;   djs_filepath()
;   get_mjd_dir()
;   sdsshead()
;   sxpar()
;
; REVISION HISTORY:
;   27-May-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro apoall, mjd=mjd, mjstart=mjstart, mjend=mjend, $
 minexp=minexp, copydir=copydir

   ;----------
   ; Set directory names, where $SPECTROLOG_DIR is used for the output dir

   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   if (NOT keyword_set(rawdata_dir)) then $
    message, 'BOSS_SPECTRO_DATA not set!'

   astrolog_dir = getenv('SPECLOG_DIR')
   if (NOT keyword_set(astrolog_dir)) then $
    message, 'SPECLOG_DIR not set!'

   spectrolog_dir = getenv('SPECTROLOG_DIR')
   if (NOT keyword_set(spectrolog_dir)) then begin
      splog, 'SPECTROLOG_DIR not set; putting outputs in input directories'
      spectrolog_dir = rawdata_dir
   endif

   ;----------
   ; Create a list of the MJD directories (as strings)

   mjdlist = get_mjd_dir(rawdata_dir, mjd=mjd, mjstart=mjstart, mjend=mjend)
   if (NOT keyword_set(mjdlist)) then begin
      splog, 'No matching MJD directories found'
      return
   endif
   splog, 'Number of MJDs = ', n_elements(mjdlist)

   ;---------------------------------------------------------------------------
   ; Loop through each input directory

   for imjd=0, n_elements(mjdlist)-1 do begin

      mjddir = mjdlist[imjd]
      inputdir = filepath('', root_dir=rawdata_dir, subdirectory=mjddir)
      plugdir = filepath('', root_dir=astrolog_dir, subdirectory=mjddir)
      outdir = filepath('', root_dir=spectrolog_dir, subdirectory=mjddir)

      splog, 'Data directory ', inputdir
      splog, 'Astrolog directory ', plugdir
      splog, 'Output directory ', outdir

      ;----------
      ; Make certain that the output directory exists.

      spawn, 'mkdir -p ' + outdir

      ;----------
      ; Find all raw FITS files in this directory

      fullname = findfile( djs_filepath('sdR*.fit*', root_dir=inputdir), $
       count=nfile)
      splog, 'Number of FITS files found: ', nfile

      if (nfile GT 0) then begin

         ;----------
         ; Find all useful header keywords

         PLATEID = lonarr(nfile)
         FLAVOR = strarr(nfile)
         CAMERAS = strarr(nfile)
         for ifile=0, nfile-1 do begin
            ; Print something since this might take a while to read all the
            ; FITS headers...
            print, format='(".",$)'

            hdr = sdsshead(fullname[ifile])

            if (size(hdr,/tname) EQ 'STRING') then begin
               PLATEID[ifile] = long( sxpar(hdr, 'PLATEID') )
               FLAVOR[ifile] = strtrim(sxpar(hdr, 'FLAVOR'),2)
               CAMERAS[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
            endif
         endfor

         ;----------
         ; Determine all the plate numbers

         platenums = PLATEID[ uniq(PLATEID, sort(PLATEID)) ]

         ;----------
         ; Loop through each plate, flavor, and camera.
         ; Must reduce arcs after flats, science after arcs.
         ; Reduce r2 camera last, so that HTML file is created at the end.

         flavlist = ['bias', 'dark', 'flat', 'arc', 'science', 'smear']
         camlist = ['b1', 'b2', 'r1', 'r2']

         for iplate=0, n_elements(platenums)-1 do begin
         for iflav=0, n_elements(flavlist)-1 do begin
         for icam=0, n_elements(camlist)-1 do begin
            ii = where(PLATEID EQ platenums[iplate] $
                   AND FLAVOR EQ flavlist[iflav] $
                   AND CAMERAS EQ camlist[icam])
            if (ii[0] NE -1) then $
             aporeduce, fileandpath(fullname[ii]), $
              indir=inputdir, outdir=outdir, $
              plugdir=plugdir, minexp=minexp, copydir=copydir, $
              /no_diskcheck, /no_lock
         endfor
         endfor
         endfor

      endif

   endfor ; End loop through input directory names (one MJD)

   return
end
;------------------------------------------------------------------------------
