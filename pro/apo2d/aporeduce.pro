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
;   minexp     - Minimum exposure time for science frames; default to 61 sec
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
;
; PROCEDURES CALLED:
;   apo_appendlog
;   apo_log2html
;   apo_plotsn
;   fits_wait()
;   quickextract()
;   quicktrace()
;   quickwave()
;   splog
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel & S. Burles, APO
;-
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
   if (NOT keyword_set(minexp)) then minexp = 61

   filer = strmid(filename,0,3)  ; root 'sdR'
   filec = strmid(filename,4,2)  ; camera name
   filee = strmid(filename,7,8)  ; exposure number

   camnames = ['b1','b2','r1','r2']

   icam = (where(filec EQ camnames))[0]

   if (filer NE 'sdR' OR icam EQ -1) then begin
      splog, 'Cannot parse FILENAME '+filename
      return
   endif

   fullname = (findfile(filepath(filename, root_dir=indir), count=ct))[0]
   if (ct NE 1) then begin
      splog, 'Found '+string(ct)+' instead of 1'
      print, fullname
      return
   endif

   ;----------
   ; Open the log file to catch WARNINGs and ABORTs.

   splgfile = filepath('splog-'+filec+'-'+filee+'.log', root_dir=outdir)
   splog, filename=splgfile, prelog=filename
   splog, 'Started at ', systime()

   ;----------
   ; Wait for a file to be fully written to disk, and exit if that doesn't
   ; happen within 3 minutes.

   if (fits_wait(fullname, deltat=10, tmax=180) EQ 0) then begin
      splog, 'File never fully written to disk: '+ fullname
      splog, /close
      return
   endif

   ;----------
   ; Find flavor, plate and MJD

   hdr = sdsshead(fullname)

   flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)

   plate = sxpar(hdr, 'PLATEID')
   platestr = string(plate, format='(i4.4)')
   mjd = sxpar(hdr, 'MJD')
   mjdstr = strtrim(string(mjd),2)

   ;----------
   ; Determine names for the FITS and HTML output log files

   logfile = filepath('logfile-' + mjdstr + '.fits', root_dir=outdir)
   htmlfile = filepath('logfile-' + mjdstr + '.html', root_dir=outdir)
    
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
;         if (NOT flatexist AND plugexist) then $ ; Only reduce 1 flat/camera?
         if (plugexist) then $
          rstruct = quicktrace(fullname, tsetfile1, fullplugfile) $
         else $
          splog, 'Unable to reduce this flat exposure'
      end

      'arc' : begin
;         if (flatexist AND (NOT arcexist)) then $ ; Only reduce 1 arc/camera?
         if (flatexist) then $
          rstruct = quickwave(fullname, tsetfile_last, wsetfile1, fflatfile1) $
          else $
           splog, 'Unable to reduce this arc exposure'
      end

      'science': begin
          exptime = sxpar(hdr, 'EXPTIME')
          outsci = filepath('sci-'+platestr+'-'+filec+'-'+filee+'.fits',$
                 root_dir=outdir)

          if (flatexist AND arcexist AND exptime GT minexp) then $
           rstruct = quickextract(tsetfile_last, wsetfile_last, fflatfile_last, $
            fullname, outsci) $
          else $
           splog, 'Unable to reduce this science exposure'
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

   if (keyword_set(rstruct)) then begin

      ;----------
      ; Find WARNINGs and ABORTs from splog file.  Recast them as string
      ; arrays that are not empty (e.g., ''), or MWRFITS will fail.

      spawn, 'grep WARNING '+splgfile, warnings
      spawn, 'grep ABORT '+splgfile, aborts

      if (warnings[0] NE '') then begin
         warnings = strtrim([warnings],2)
      endif else begin
         warnings = [' ']
      endelse

      if (aborts[0] NE '') then begin
         aborts = strtrim([aborts],2)
      endif else begin
         aborts = [' ']
      endelse

      ;----------
      ; Get the universal time (UT) from the header in the format '12:34',
      ; then convert to Mountain standard time (MST), which is 7 (or 17)
      ; hours different from UT.

      ut = strmid(sxpar(hdr, 'TAIHMS'),0,5)
      mst = string((long(strmid(ut,0,2))+17) MOD 24,format='(i2.2)' ) $
       + strmid(ut,2,3)

      ; The following prevents a crash in MWRFITS.
      if (NOT keyword_set(shortplugfile)) then shortplugfile = ' '

      rstruct = create_struct('FILENAME', filename, $
                              'PLUGFILE', shortplugfile, $
                              'MJD', mjd, $
                              'PLATE', plate, $
                              'EXPNUM', filee, $
                              'EXPTIME', sxpar(hdr, 'EXPTIME'), $
                              'FLAVOR', flavor, $
                              'CAMERA', camnames[icam], $
                              'MST', mst, $
                              rstruct, $
                              'WARNINGS', warnings, $
                              'ABORTS', aborts )

      apo_appendlog, logfile, rstruct
   endif

   ;----------
   ; After being passed any 'r2' frame, and if it was reduced,
   ; we re-generate the HTML file for all plates and the S/N plot for
   ; this plate.
   ; Optionally copy it to the directory specified by COPYDIR.
   ; Make an additional copy of the HTML file called 'logsheet-current.html'.

;   if (camnames[icam] EQ 'r2' AND keyword_set(rstruct)) then begin
; Instead, create the HTML file after any reduced frame.
   if (keyword_set(rstruct)) then begin

      if (myflavor EQ 'science') then begin
         plotfile = filepath('snplot-'+mjdstr+'-'+platestr+'.ps', $
          root_dir=outdir)
         splog, 'Generating S/N plot '+plotfile
         apo_plotsn, logfile, plate, plugdir=plugdir, plotfile=plotfile
      endif

      splog, 'Generating HTML file '+htmlfile
      apo_log2html, logfile, htmlfile

      if (keyword_set(copydir)) then begin
         splog, 'Copying files to ', copydir
         spawn, 'scp1 ' + htmlfile + ' ' + copydir
         spawn, 'scp1 ' + htmlfile + ' ' + $
          filepath('logfile-current.html', root_dir=copydir)
         spawn, 'scp1 ' + plotfile + ' ' + copydir
         spawn, 'scp1 ' + logfile  + ' ' + copydir
         splog, 'Done.'
      endif
   endif

   ;----------
   ; Close splog file

   splog, 'Finished at ', systime()
   splog, /close

   return
end
;------------------------------------------------------------------------------
