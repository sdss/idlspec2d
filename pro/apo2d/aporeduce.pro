;+
; NAME:
;   aporeduce
;
; PURPOSE:
;   Quick on-the-mountain reduction pipeline for 1 file at a time.
;
; CALLING SEQUENCE:
;   aporeduce, filename, [ indir=, outdir=, $
;    plugfile=, plugdir=, minexp=, $
;    copydir=copydir ]
;
; INPUTS:
;   filename   - Raw spectroscopic image file name(s) of any flavor; this
;                can be an array of file names, but cannot include wildcards
;
; OPTIONAL INPUTS:
;   indir      - Input directory for FILENAME; default to './'
;   outdir     - Output directory for reduced data and log files;
;                default to INDIR
;   plugfile   - Name of plugmap file (Yanny parameter file); default to
;                'plPlugMapM-'+NAME+'.par', where NAME is taken from that
;                keyword in the file FILENAME
;   plugdir    - Input directory for PLUGFILE; default to INDIR
;   minexp     - Minimum exposure time for science frames; default to 0 sec
;                so that any frame with a non-negative exposure time is
;                reduced.
;   copydir    - If set, then copy the output log files to this directory using
;                "scp1" copy.  Make an additional copy of the HTML file
;                called 'logsheet-current.html'.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   After reducing any 'r2' frame, we re-generate the HTML file and optionally
;   copy it to the file specified by COPYDIR.
;
; EXAMPLES:
;
; BUGS:
;   scp1 does not exist on sos.apo.nmsu.edu, reverted to scp
;
; INTERNAL SUPPORT ROUTINES:
;   apo_diskcheck
;
; PROCEDURES CALLED:
;   apo_appendlog
;   apo_log2html
;   apo_plotsn
;   djs_filepath()
;   fits_wait()
;   idlspec2d_version()
;   idlutils_version()
;   quickextract()
;   quicktrace()
;   quickwave()
;   splog
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel & S. Burles, APO
;-
;------------------------------------------------------------------------------
; Check disk space on the input or output disk.
pro apo_diskcheck, dirname

   if (NOT keyword_set(dirname)) then return

   spawn, 'df -k '+dirname, dfout
   if (size(dfout,/tname) EQ 'STRING') then begin
      dfout_entry = dfout[n_elements(dfout)-1]
      if (dfout_entry NE '') then begin
         perc  = str_sep(dfout_entry,'%')
         percentfull = long(strmid(perc[0],rstrpos(perc[0],' ')))
         if (percentfull GT 95) then $
          splog, 'WARNING: SOS disk '+dirname+' is ' $
           +strtrim(string(percentfull),2)+'% full'
      endif else splog, 'Warning: Could not check disk space on '+dirname
   endif else splog, 'Warning: Could not check disk space on '+dirname

   return
end
;------------------------------------------------------------------------------
pro aporeduce, filename, indir=indir, outdir=outdir, $
 plugfile=plugfile, plugdir=plugdir, minexp=minexp, $
 copydir=copydir

   if (n_params() LT 1) then begin
      print, 'Syntax - aporeduce, filename, [ indir=, outdir=, $'
      print, ' plugfile=, plugdir=, minexp= ]'
      return
   endif

   if (size(filename, /tname) NE 'STRING') then begin
      splog, 'FILENAME is not a string'
      return
   endif

   ;----------
   ; If multiple file names are passed, then call this script recursively
   ; for each file.

   if (n_elements(filename) GT 1) then begin
      for ifile=0, n_elements(filename)-1 do $
       aporeduce, filename[ifile], indir=indir, outdir=outdir, $
       plugfile=plugfile, plugdir=plugdir, minexp=minexp, $
       copydir=copydir
      return
   endif else begin
      filename = filename[0] ; Convert from an array to a scalar.
   endelse

   if (NOT keyword_set(indir)) then indir = './'
   if (NOT keyword_set(plugdir)) then plugdir = indir
   if (NOT keyword_set(outdir)) then outdir = indir
   if (n_elements(minexp) EQ 0) then minexp = 0

   filer = strmid(filename,0,3)  ; root 'sdR'
   filec = strmid(filename,4,2)  ; camera name
   filee = strmid(filename,7,8)  ; exposure number

   camnames = ['b1','b2','r1','r2']

   icam = (where(filec EQ camnames))[0]

   if (filer NE 'sdR' OR icam EQ -1) then begin
      splog, 'Cannot parse FILENAME '+filename
      return
   endif

; No need to do this test, since the FITS_WAIT() function below will do it.
;   fullname = (findfile(filepath(filename, root_dir=indir), count=ct))[0]
;   if (ct NE 1) then begin
;      splog, 'Found '+string(ct)+' instead of 1'
;      print, fullname
;      return
;   endif

   ;----------
   ; Open the log file to catch WARNINGs and ABORTs.

   splgfile = filepath('splog-'+filec+'-'+filee+'.log', root_dir=outdir)
   splog, filename=splgfile, prelog=filename
   splog, 'Log file ' + splgfile + ' opened' + systime()
   t0 = systime(1)

   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   splog, 'DISPLAY=' + getenv('DISPLAY')

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   ;----------
   ; Check disk space on both the input and the output disk.

   apo_diskcheck, indir
   apo_diskcheck, outdir

   ;----------
   ; Wait for an input FITS file to be fully written to disk, and exit
   ; if that doesn't happen within 3 minutes.

   fullname = djs_filepath(filename, root_dir=indir)
   spawn, 'ls -l '+fullname, lsstring
   splog, 'DIRLIST '+lsstring

   if (fits_wait(fullname, deltat=10, tmax=180) EQ 0) then begin
      splog, 'WARNING: File never fully written to disk: '+ fullname
      splog, /close
      return
   endif

   spawn, 'ls -l '+fullname, lsstring
   splog, 'DIRLIST '+lsstring

   ;----------
   ; Find flavor, plate and MJD

   splog, 'Using SDSSHEAD to read FITS header'
   hdr = sdsshead(fullname)

   flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)
   plate = sxpar(hdr, 'PLATEID')
   platestr = string(plate, format='(i4.4)')
   mjd = sxpar(hdr, 'MJD')
   mjdstr = strtrim(string(mjd),2)

   splog, 'FLAVOR=', flavor, ' PLATEID=', plate, ' MJD=', mjd

   ;----------
   ; Determine names for the FITS and HTML output log files

   logfile = filepath('logfile-' + mjdstr + '.fits', root_dir=outdir)
   htmlfile = filepath('logfile-' + mjdstr + '.html', root_dir=outdir)
   currentfile = filepath('logfile-current.html', root_dir=outdir)
    
   ;----------
   ; Find the full name of the plugmap file

   if (NOT keyword_set(plugfile)) then begin
       name = strtrim(sxpar(hdr,'NAME'),2)
       ; This string should contain PLATE-MJD-PLUGID, but it may not
       ; in some of the early data, in which case we're search using wildcards
       if (strlen(name) LT 13) then name = '*' + name + '*'
       plugfile = 'plPlugMapM-'+name+'.par'
   endif
   fullplugfile = findfile( filepath(plugfile, root_dir=plugdir) )
   ; If we found several plugmap files (using wildcards), take the most
   ; recent as determined by simply doing an ASCII sort of the file names.
   if (n_elements(fullplugfile) EQ 1) then fullplugfile = fullplugfile[0] $
    else fullplugfile = fullplugfile[ (reverse(sort(fullplugfile)))[0] ]

   ;----------
   ; Construct the names of the flat and arc output files if we generate
   ; them from this exposure.

   tsetfile1 = filepath( $
    'tset-'+mjdstr+'-'+platestr+'-'+filee+'-'+filec+'.fits', $
    root_dir=outdir)
   wsetfile1 = filepath( $
    'wset-'+mjdstr+'-'+platestr+'-'+filee+'-'+filec+'.fits', $
    root_dir=outdir)
   fflatfile1 = filepath( $
    'fflat-'+mjdstr+'-'+platestr+'-'+filee+'-'+filec+'.fits', $
    root_dir=outdir)

   ;----------
   ; Determine if a flat or arc for this plate has already been reduced,
   ; and test if the plugmap file exists.
   ; Use the last flat and arc files on disk.

   tsetfiles = findfile(filepath( $
    'tset-'+mjdstr+'-'+platestr+'-*-'+filec+'.fits', $
    root_dir=outdir))
   wsetfiles = findfile(filepath( $
    'wset-'+mjdstr+'-'+platestr+'-*-'+filec+'.fits', $
    root_dir=outdir))
   fflatfiles = findfile(filepath( $
    'fflat-'+mjdstr+'-'+platestr+'-*-'+filec+'.fits', $
    root_dir=outdir))

   tsetfile_last = (reverse(tsetfiles))[0]
   wsetfile_last = (reverse(wsetfiles))[0]
   fflatfile_last = (reverse(fflatfiles))[0]

   plugexist = keyword_set(fullplugfile)
   flatexist = keyword_set(tsetfile_last) AND $
    keyword_set( findfile(tsetfile_last) )
   arcexist = keyword_set(wsetfile_last) AND $
    keyword_set( findfile(wsetfile_last) )
   splog, 'PLUGEXIST = ', plugexist
   splog, 'FLATEXIST = ', flatexist
   splog, 'ARCEXIST = ', arcexist

   ;----------
   ; Reduce file depending on its flavor: bias/dark, flat, arc, or science/smear

   rstruct = 0
   myflavor = flavor
   if (myflavor EQ 'smear') then myflavor = 'science'
   if (myflavor EQ 'dark') then myflavor = 'bias'
   case myflavor of
      'bias' : begin
         rstruct = quickbias(fullname)
      end

      'flat' : begin
         if (plugexist) then begin
            rstruct = quicktrace(fullname, tsetfile1, fullplugfile)
         endif else begin
            splog, 'ABORT: Unable to reduce this flat exposure (need plug-map)'
         endelse
      end

      'arc' : begin
         if (flatexist) then begin
            rstruct = quickwave(fullname, tsetfile_last, wsetfile1, fflatfile1)
         endif else begin
             splog, 'ABORT: Unable to reduce this arc exposure (need flat)'
         endelse
      end

      'science': begin
          exptime = sxpar(hdr, 'EXPTIME')
          outsci = filepath('sci-'+platestr+'-'+filec+'-'+filee+'.fits',$
                 root_dir=outdir)

          if (flatexist AND arcexist AND exptime GE minexp) then begin
             rstruct = quickextract(tsetfile_last, wsetfile_last, $
              fflatfile_last, fullname, outsci)
          endif else begin
             if (NOT keyword_set(flatexist)) then $
              splog, 'ABORT: Unable to reduce this science exposure (need flat)'
             if (NOT keyword_set(arcexist)) then $
              splog, 'ABORT: Unable to reduce this science exposure (need arc)'
             if (exptime LT minexp) then $
              splog, 'ABORT: Exposure time = ' + string(exptime) $
               + ' < ' + string(minexp)
          endelse
       end

       else : begin
          splog, 'Unknown flavor: ', flavor
       end
   endcase

   ;----------
   ; Append to binary FITS log file a structure with info for this frame
   ; Lock the file to do this.

   i = rstrpos(fullplugfile,'/')
   if (i[0] EQ -1) then shortplugfile = fullplugfile $
    else shortplugfile = strmid(fullplugfile,i+1)

   ;----------
   ; Find WARNINGs and ABORTs from splog file.  Recast them as string
   ; arrays that are not empty (e.g., ''), or MWRFITS will fail.

   spawn, 'grep -e WARNING -e ABORT '+splgfile, tstring
   if (keyword_set(tstring)) then begin
      tstruct = create_struct('FILENAME', filename, $
                              'MJD', mjd, $
                              'PLATE', plate, $
                              'EXPNUM', filee, $
                              'TEXT', '' )
      tstruct = replicate(tstruct, n_elements(tstring))
      tstruct.text = tstring
   endif

   if (keyword_set(rstruct)) then begin
      ;----------
      ; Get the universal time (UT) from the header in the format '12:34',
      ; then convert to Mountain standard time (MST), which is 7 (or 17)
      ; hours different from UT.
      ; Instead, only store the TAI time which we will convert to UT time
      ; on-the-fly.

;      ut = strmid(sxpar(hdr, 'TAIHMS'),0,5)
;      mst = string((long(strmid(ut,0,2))+17) MOD 24,format='(i2.2)' ) $
;      + strmid(ut,2,3)
      tai = sxpar(hdr, 'TAI')

      ; The following prevents a crash in MWRFITS.
      if (NOT keyword_set(shortplugfile)) then shortplugfile = ' '

      rstruct = create_struct('FILENAME', string(filename), $
                              'PLUGFILE', string(shortplugfile), $
                              'MJD', long(mjd), $
                              'PLATE', long(plate), $
                              'EXPNUM', long(filee), $
                              'EXPTIME', double(sxpar(hdr, 'EXPTIME')), $
                              'FLAVOR', string(flavor), $
                              'CAMERA', string(camnames[icam]), $
                              'TAI', double(tai), $
                              rstruct )
   endif

   if (keyword_set(tstruct) OR keyword_set(rstruct)) then begin
      splog, 'Appending to FITS log file '+logfile
      apo_appendlog, logfile, rstruct, tstruct
      splog, 'Done with append'
   endif

   ;----------
   ; After being passed any 'r2' frame, and if it was reduced,
   ; we re-generate the HTML file for all plates and the S/N plot for
   ; this plate.
   ; Optionally copy it to the directory specified by COPYDIR.
   ; Make an additional copy of the HTML file called 'logsheet-current.html'.

;   if (camnames[icam] EQ 'r2' AND keyword_set(rstruct)) then begin
; Instead, create the HTML file after any reduced frame.
   if (keyword_set(rstruct) OR keyword_set(tstruct)) then begin

      if (myflavor EQ 'science') then begin
         plotfile = filepath('snplot-'+mjdstr+'-'+platestr+'.ps', $
          root_dir=outdir)
         splog, 'Generating S/N plot '+plotfile
         apo_plotsn, logfile, plate, plugdir=plugdir, plotfile=plotfile
         splog, 'Done generating plot'
      endif

      splog, 'Generating HTML file '+htmlfile
      apo_log2html, logfile, htmlfile
      splog, 'Done generating HTML file'

      ; Generate a copy of the HTML file, 'logsheet-current.html',
      ; that includes the Java script to auto-load the page every 60 seconds.

      squote = "\'"
      addstring = $
       '<BODY ONLOAD=\"timerID=setTimeout(' $
       +squote+'location.reload(true)'+squote+',60000)\">'
      sedcommand = '"s/<\/HEAD>/<\/HEAD>'+addstring+'/g"'
      setenv, 'SHELL=bash'
      spawn, 'sed -e ' + sedcommand + ' ' + htmlfile + ' > ' + currentfile

       ; scp1 does not work on sos.apo, let's switch to scp

      if (keyword_set(copydir)) then begin
         splog, 'Copying files to ', copydir
         spawn, 'scp ' + htmlfile + ' ' + copydir
         spawn, 'scp ' + currentfile + ' ' + copydir
         spawn, 'scp ' + logfile  + ' ' + copydir
         if (keyword_set(plotfile)) then $
          spawn, 'scp ' + plotfile + ' ' + copydir
         splog, 'Done.'
      endif
   endif

   ;----------
   ; Close splog file

   splog, 'Elapsed time = ', systime(1)-t0
   splog, 'Finished at ', systime()
   splog, /close

   return
end
;------------------------------------------------------------------------------
