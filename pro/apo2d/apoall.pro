;+
; NAME:
;   apoall
;
; PURPOSE:
;   Run APOREDUCE on one or many nights of data.
;
; CALLING SEQUENCE:
;   apoall, [ rawdir, astrolog=, flatdir=, mjd=, mjstart=, mjend=, $
;    minexp=, copydir= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   rawdir     - Search for raw data files in RAWDIR/MJD/*.
;                This should be an absolute file path, and we default to
;                '/usr/sdss/data05/spectro/rawdata'.
;   astrolog   - Search for plug-map files in ASTROLOG/MJD/*.
;                This should be an absolute file path, and we default to
;                '/usr/sdss/data05/spectro/astrolog'.
;   flatdir    - Directory for pixel flat files.  For now, default
;                to 'pixflat'.
;   mjd        - Look for raw data files in RAWDIR/MJD; default to '*' to
;                search all subdirectories.  Note that this need not be
;                integer-valued, but could be for example '51441_test'.
;   mjstart    - Starting MJD.
;   mjend      - Ending MJD.
;   minexp     - Minimum exposure time for science frames; default to 0 sec.
;   copydir    - Copy the output log files to this directory; default to
;                the current directory.
;
; OUTPUT:
;
; COMMENTS:
;   The files are sorted before being sent to APOREDUCE.  For each plate,
;   reduce the flats, arcs, and science frames in that order.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   aporeduce
;   get_mjd_dir()
;   sdsshead()
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   27-May-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

pro apoall, rawdir, astrolog=astrolog, flatdir=flatdir, mjd=mjd, $
 mjstart=mjstart, mjend=mjend, minexp=minexp, copydir=copydir

   ;----------
   ; Set directory names RAWDIR, ASTROLOG, FLATDIR

   if (NOT keyword_set(copydir)) then cd, current=copydir

   if (NOT keyword_set(rawdir)) then begin
      if ((findfile('/usr/sdss/data05/spectro/rawdata'))[0] NE '') then $
       rawdir = '/usr/sdss/data05/spectro/rawdata' $
      else if ((findfile('/home/schlegel/data/rawdata'))[0] NE '') then $
       rawdir = '/scr0/data/rawdata' $
      else if ((findfile('rawdata'))[0] NE '') then $
       rawdir = './rawdata' $
      else begin
        print, 'Must specify RAWDIR'
        return
      endelse
   endif

   if (NOT keyword_set(astrolog)) then $
    astrolog = strmid(rawdir, 0, rstrpos(rawdir,'/')+1) + 'astrolog'
   if (NOT keyword_set(flatdir)) then flatdir = 'pixflat'
   if (n_elements(minexp) EQ 0) then minexp = 0

   ;  This trick expands directories
   cd, rawdir, current=olddir
   cd, olddir, current=rawdir
   cd, astrolog, current=olddir
   cd, olddir, current=astrolog

   ;----------
   ; Create a list of the MJD directories (as strings)

   mjdlist = get_mjd_dir(rawdir, mjd=mjd, mjstart=mjstart, mjend=mjend)

   ;---------------------------------------------------------------------------
   ; Loop through each input directory

   for imjd=0, N_elements(mjdlist)-1 do begin

      mjddir = mjdlist[imjd]
      inputdir = filepath('', root_dir=rawdir, subdirectory=mjddir)
      plugdir = filepath('', root_dir=astrolog, subdirectory=mjddir)

      print, 'Data directory ', inputdir
      print, 'Astrolog directory ', plugdir

      ; Find all raw FITS files in this directory
      cd, inputdir, current=olddir
      fullname = findfile('sdR*.fit*', count=nfile)
      cd, olddir

      print, 'Number of FITS files found: ', nfile

      if (nfile GT 0) then begin

         ;----------
         ; Find all useful header keywords

         PLATEID = lonarr(nfile)
         FLAVOR = strarr(nfile)
         CAMERAS = strarr(nfile)
         for i=0, nfile-1 do begin
            ; Print something since this might take a while to read all the
            ; FITS headers...
            print, format='(".",$)'

            hdr = sdsshead(filepath(fullname[i], root_dir=inputdir))

            if (size(hdr,/tname) EQ 'STRING') then begin
               PLATEID[i] = long( sxpar(hdr, 'PLATEID') )
               FLAVOR[i] = strtrim(sxpar(hdr, 'FLAVOR'),2)
               CAMERAS[i] = strtrim(sxpar(hdr, 'CAMERAS'),2)
            endif
         endfor

         ;----------
         ; Determine all the plate numbers

         platenums = PLATEID[ uniq(PLATEID, sort(PLATEID)) ]

         ;----------
         ; Loop through each plate, flavor, and camera.
         ; Must reduce arcs after flats, science after arcs.
         ; Reduce r2 camera last, so that HTML file is created at the end.

         flavlist = ['flat', 'arc', 'science']
         camlist = ['b1', 'b2', 'r1', 'r2']

         for iseq=0, n_elements(platenums)-1 do begin
         for iflav=0, n_elements(flavlist)-1 do begin
         for icam=0, n_elements(camlist)-1 do begin

            ii = where(PLATEID EQ platenums[iseq] $
                   AND FLAVOR EQ flavlist[iflav] $
                   AND CAMERAS EQ camlist[icam])
            if (ii[0] NE -1) then $
             aporeduce, fullname[ii], indir=inputdir, outdir=inputdir, $
              plugdir=plugdir, minexp=minexp, copydir=copydir

         endfor
         endfor
         endfor

      endif

   endfor ; End loop through input directory names (one MJD)

   return
end
;------------------------------------------------------------------------------
