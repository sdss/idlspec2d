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
;   plugdir    - Input directory for PLUGFILE; default to '.'
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
;   apo_log2html
;   apo_plotsn
;   djs_lockfile()
;   djs_unlockfile
;   mwrfits
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
   endif

   if (NOT keyword_set(indir)) then indir = './'
   if (NOT keyword_set(plugdir)) then plugdir = './'
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
   ; Find flavor, plate and MJD

   hdr = headfits(fullname)

   flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)

   plate = sxpar(hdr, 'PLATEID')
   platestr = string(plate, format='(i4.4)')
   mjd = sxpar(hdr, 'MJD')
   mjdstr = strtrim(string(mjd),2)

   ;----------
   ; Fix up some header information from early data

   if (flavor EQ 'target') then flavor = 'science'

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
   ; Determine if a flat or arc for this plate has already been reduced,
   ; and test if the plugmap file exists.

   tsetfile = filepath('tset-'+platestr+'-'+filec+'.fits', root_dir=outdir)
   wsetfile = filepath('wset-'+platestr+'-'+filec+'.fits', root_dir=outdir)
   fflatfile = filepath('fflat-'+platestr+'-'+filec+'.fits', root_dir=outdir)

   plugexist = keyword_set(fullplugfile)
   flatexist = keyword_set( findfile(tsetfile) )
   arcexist = keyword_set( findfile(wsetfile) )

   ;----------
   ; Reduce file depending on its flavor: flat, arc, or science

   rstruct = 0
   case flavor of
      'flat' : begin
;         if (NOT flatexist AND plugexist) then $ ; Only reduce 1 flat/camera?
         if (plugexist) then $
          rstruct = quicktrace(fullname, tsetfile, fullplugfile) $
         else $
          splog, 'Unable to reduce this flat exposure'
      end

      'arc' : begin
;         if (flatexist AND (NOT arcexist)) then $ ; Only reduce 1 arc/camera?
         if (flatexist) then $
          rstruct = quickwave(fullname, tsetfile, wsetfile, fflatfile) $
          else $
           splog, 'Unable to reduce this arc exposure'
      end

      'science': begin
          exptime = sxpar(hdr, 'EXPTIME')
          outsci = filepath('sci-'+platestr+'-'+filec+'-'+filee+'.fits',$
                 root_dir=outdir)

          if (flatexist AND arcexist AND exptime GT minexp) then $
           rstruct = quickextract(tsetfile, wsetfile, fflatfile, $
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

   if (keyword_set(rstruct)) then begin
      while(djs_lockfile(logfile) EQ 0) do wait, 1
      rstruct = create_struct('MJD', mjd, $
                              'PLATE', plate, $
                              'EXPNUM', filee, $
                              rstruct )
      mwrfits, rstruct, logfile
      djs_unlockfile, logfile
   endif

   ;----------
   ; After being passed any 'r2' frame, and if it was reduced,
   ; we re-generate the HTML file for all plates and the S/N plot for
   ; this plate.
   ; Optionally copy it to the directory specified by COPYDIR.
   ; Make an additional copy of the HTML file called 'logsheet-current.html'.

   if (camnames[icam] EQ 'r2' AND keyword_set(rstruct)) then begin
      wait, 10

      plotfile = filepath('snplot-'+mjdstr+'-'+platestr+'.ps', root_dir=outdir)
      apo_plotsn, logfile, plate, plugfile=fullplugfile, plotfile=plotfile

      apo_log2html, logfile, htmlfile

      if (keyword_set(copydir)) then begin
         splog
         splog, 'Copying files to ', copydir
         spawn, 'scp1 ' + htmlfile + ' ' + copydir
         spawn, 'scp1 ' + htmlfile + ' ' + $
          filepath('logfile-current.html', root_dir=copydir)
         spawn, 'scp1 ' + plotfile + ' ' + copydir
         spawn, 'scp1 ' + logfile  + ' ' + copydir
         splog, 'Done.'
      endif
   endif

   return
end
;------------------------------------------------------------------------------
