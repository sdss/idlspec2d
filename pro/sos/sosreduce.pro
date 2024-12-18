;+
; NAME:
;   sosreduce
;
; PURPOSE:
;   Quick on-the-mountain reduction pipeline for 1 file at a time.
;
; CALLING SEQUENCE:
;   sosreduce, filename, [ indir=, outdir=, $
;    plugfile=, plugdir=, minexp=, $
;    copydir=, /no_diskcheck, /no_lock ]
;
; INPUTS:
;   filename   - Raw spectroscopic image file name(s) of any flavor; this
;                can be an array of file names, but cannot include wildcards
;
; OPTIONAL INPUTS:
;   indir      - Input directory for FILENAME; default to './';
;                conventionally will explicitly contain $MJD in the name
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
;                "scp" copy (not "scp1" any longer).  Make an additional copy
;                of the HTML file called 'logsheet-current.html'.
;   no_diskcheck- If set, then do not do the check for filling input or
;                output disks.  (This option is always set by the APOALL proc).
;   no_lock    - If set, then do not create lock files for the input files;
;                this option is useful for calls from APOALL.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   After reducing any 'r2' frame, we re-generate the HTML file and optionally
;   copy it to the file specified by COPYDIR.
;
;   The copy of the file "logsheet-current.html" also has a line of Java-script
;   added that does an auto-refresh every 60 seconds.
;
; EXAMPLES:
;
; BUGS:
;   scp1 does not exist on sos.apo.nmsu.edu, reverted to scp
;
; INTERNAL SUPPORT ROUTINES:
;   sos_diskcheck
;
; PROCEDURES CALLED:
;   sos_appendlog
;   sos_log2html
;   sos_plotsn
;   djs_filepath()
;   fits_wait()
;   get_tai
;   idlspec2d_version()
;   idlutils_version()
;   quickbias()
;   quickextract()
;   quicktrace()
;   quickwave()
;   splog
;   tai2airmass()
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel & S. Burles, APO
;-
;------------------------------------------------------------------------------
function latest_flat, flatlist, this_expid = this_expid
  if flatlist[0] EQ '' then return, flatlist
  expids=strarr(n_elements(flatlist))
  foreach flat, flatlist, i do begin
    expids[i]=(STRSPLIT(flat, '-', /EXTRACT))[-2]
  endforeach
  if not keyword_set(this_expid) then begin
    return, flatlist[where(expids EQ max(expids))]
  endif else begin
    mindiff = min(abs((long64(expids) - long64(this_expid))), idxmin)
    return, flatlist[idxmin]
  endelse
end


;------------------------------------------------------------------------------
pro find_cal, fieldstr, mjdstr, filee, filec, tsetfile_last = tsetfile_last,$
             wsetfile_last = wsetfile_last, fflatfile_last = fflatfile_last,$
             outdir = outdir
   ; Determine if there is a reduced flat or arc for this field. If so use the latest
   ; for that field, else uses the cloest in expid
   
   tsetfiles = findfile(filepath( $
        'tset-'+mjdstr+'-'+fieldstr+'-*-'+filec+'.fits', root_dir=outdir))
   if (n_elements(tsetfiles) gt 0) and (tsetfiles[0] ne '') then begin
      tsetfile_last = latest_flat(tsetfiles)
   endif else begin
      tsetfiles = findfile(filepath('tset-'+mjdstr+'-*-*-'+filec+'.fits', root_dir=outdir))
      tsetfile_last = latest_flat(tsetfiles, this_expid=filee)
   endelse
   
   wsetfile_last = 0
   fflatfile_last = 0
   if not keyword_set(nocal) then begin
        wsetfiles = findfile(filepath( $
            'wset-'+mjdstr+'-'+fieldstr+'-*-'+filec+'.fits',root_dir=outdir))
        if (n_elements(wsetfiles) gt 0)  and (wsetfiles[0] ne '') then begin
            wsetfile_last = max(wsetfiles)
        endif else wsetfile_last = 0
        
        fflatfiles = findfile(filepath( $
            'fflat-'+mjdstr+'-'+fieldstr+'-*-'+filec+'.fits', root_dir=outdir))
        if (n_elements(fflatfiles) gt 0) and (fflatfiles[0] ne '') then begin
            fflatfile_last = max(fflatfiles)
        endif else fflatfile_last = 0
   endif

   if not keyword_set(wsetfile_last) then begin
        wsetfiles = findfile(filepath( $
            'wset-'+mjdstr+'-*-*-'+filec+'.fits', root_dir=outdir))
        wsetfile_last = latest_flat(wsetfiles, this_expid=filee)
   endif
   if not keyword_set(fflatfile_last) then begin
        fflatfiles = findfile(filepath( $
            'fflat-'+mjdstr+'-*-*-'+filec+'.fits', root_dir=outdir))
        fflatfile_last = latest_flat(fflatfiles, this_expid=filee)
   endif
   return
end


;------------------------------------------------------------------------------

; Check disk space on the input or output disk.
pro sos_diskcheck, dirname

   if (NOT keyword_set(dirname)) then return

   spawn, 'df -k '+dirname, dfout
   if (size(dfout,/tname) EQ 'STRING') then begin
      dfout_entry = dfout[n_elements(dfout)-1]
      if (dfout_entry NE '') then begin
         perc  = str_sep(dfout_entry,'%')
         percentfull = long(strmid(perc[0],strpos(perc[0],' ',/reverse_search)))
         if (percentfull GT 95) then $
          splog, 'WARNING: SOS disk '+dirname+' is ' $
           +strtrim(string(percentfull),2)+'% full'
      endif else splog, 'Warning: Could not check disk space on '+dirname
   endif else splog, 'Warning: Could not check disk space on '+dirname

   return
end
;------------------------------------------------------------------------------
pro sosreduce, filename, indir=indir, outdir=outdir, $
 plugfile=plugfile, plugdir=plugdir, minexp=minexp, nocal=nocal,$
 copydir=copydir,  no_diskcheck=no_diskcheck, no_lock=no_lock, $
 fps=fps, noreject=noreject, sdssv_sn2=sdssv_sn2, sn2_15=sn2_15,$
 arc2trace=arc2trace, forcea2t=forcea2t
   if (n_params() LT 1) then begin
      doc_library, 'sosreduce'
      return
   endif

   if (size(filename, /tname) NE 'STRING') then begin
      splog, 'FILENAME is not a string'
      return
   endif

   ;----------
   ; Create the output directory if it does not exist

   if (keyword_set(outdir)) then begin
      if (file_test(outdir, /directory) EQ 0) then begin
         spawn, '\mkdir -p '+outdir
      endif
      if (file_test(outdir, /directory, /write) EQ 0) then begin
         splog, 'OUTDIR not a writeable directory '+outdir
         return
      endif
   endif

   ;----------
   ; If multiple file names are passed, then call this script recursively
   ; for each file.

   if (n_elements(filename) GT 1) then begin
      for ifile=0, n_elements(filename)-1 do $
       sosreduce, filename[ifile], indir=indir, outdir=outdir, $
       plugfile=plugfile, plugdir=plugdir, minexp=minexp, $
       copydir=copydir, no_diskcheck=no_diskcheck, no_lock=no_lock, $
       fps=fps, noreject=noreject
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

   if getenv('OBSERVATORY') eq 'LCO' then begin
	camnames = ['b2','r2'] 
	lco = 1
   endif else begin
	camnames = ['b1','r1']
	lco = 0
   endelse
   icam = (where(filec EQ camnames))[0]

   if (filer NE 'sdR' OR icam EQ -1) then begin
      splog, 'Cannot parse FILENAME '+filename
      return
   endif

   if keyword_set(forcea2t) then arc2trace = 1

   ;----------
   ; Open the log file to catch WARNINGs and ABORTs.

   splgfile = filepath('splog-'+filec+'-'+filee+'.log', root_dir=outdir)
   splog, filename=splgfile, prelog=filename
   splog, 'Log file ' + splgfile + ' opened ' + systime()
   t0 = systime(1)

   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   splog, 'DISPLAY=' + getenv('DISPLAY')

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()

   ;----------
   ; Check disk space on both the input and the output disk.

   if (NOT keyword_set(no_diskcheck)) then begin
      sos_diskcheck, indir
      sos_diskcheck, outdir
   endif

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

   do_lock = keyword_set(no_lock) EQ 0

   ;----------
   ; Find flavor, hartmann status, plate and MJD

   splog, 'Using SDSSHEAD() to read FITS header'
   hdr = sdsshead(fullname, do_lock=do_lock)
   hartmann=strtrim(sxpar(hdr,'HARTMANN'),2)
   flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)
   if (NOT keyword_set(fps)) then begin
       ;The configuration id is set as the plateid, this is only for the sdss-v plate program
       config = sxpar(hdr, 'PLATEID')
       confstr = config_to_string(config)
   endif else begin
       config = sxpar(hdr, 'CONFID')
       confstr = strtrim(config,2)
       fieldid = sxpar(hdr, 'FIELDID')
       camera = strtrim(sxpar(hdr, 'CAMERAS'),2)
   endelse
   splog,config
   confstr = config_to_string(config)
   cartid = sxpar(hdr, 'CARTID')
   mjd = sxpar(hdr, 'MJD')
   mjdstr = strtrim(string(mjd),2)
   exposure = long( sxpar(hdr, 'EXPOSURE') )
   threshold=2000.

   if (NOT keyword_set(fps)) then begin
       splog, 'FLAVOR=', flavor, ' PLATEID=', config, ' MJD=', mjd
   endif else begin
       splog, 'FLAVOR=', flavor, ' CONFID=', config, ' MJD=', mjd
   endelse

   platetype0='BHM';'BOSSHALF'
   platetype1='BHM&MWM';'APOGEE-BOSS'
;; Modified by Vivek for omitting the manga plates. 
;; Modified by HJIM for the SDSS-V. 
   if (flavor NE 'dark') then begin; or flavor NE 'bias') then begin 
   if (flavor NE 'bias') then begin
   ; Check CONFIGTYP for BOSS or EBOSS (e.g. not MANGA)
   ; If keyword is missing (older data), assume this is BOSS
   platetype = sxpar(hdr, 'PLATETYP', count=nhdr)
       if (nhdr GT 0) then begin
         platetype = strupcase(strtrim(platetype,2))
         if (platetype NE platetype0) && (platetype NE platetype1) then begin
            splog, 'Skipping ' + platetype + ' plate ', plateid,' exposure ', exposure
            flavor = 'unknown'
         endif
      endif
   endif
   endif
   ;----------
   ; Determine names for the FITS and HTML output log files

   logfile = filepath('logfile-' + mjdstr + '.fits', root_dir=outdir)
   htmlfile = filepath('logfile-' + mjdstr + '.html', root_dir=outdir)
   currentfile = filepath('logfile-current.html', root_dir=outdir)
    

   ;----------
   ; Find the full name of the plugmap file

   if (NOT keyword_set(plugfile)) then begin
       ; This string should contain PLATE-MJD-PLUGID, but it may not
       ; in some of the early data, in which case we're search using wildcards
       if (NOT keyword_set(fps)) then begin
           name = strtrim(sxpar(hdr,'NAME'),2)
           if (strlen(name) LT 13) then name = '*' + name + '*'
           plugfile = 'plPlugMapM-'+name+'.par'
       endif else begin
           name = strtrim(sxpar(hdr,'CONFID'),2)
           confile = (findfile(filepath('confSummaryF-' + name + '.par',$
                                                root_dir=plugdir, subdir=['*','*']), count=ct))[0]
           if ct ne 0 then plugfile = 'confSummaryF-'+name+'.par' $
                      else plugfile = 'confSummary-'+name+'.par'
       endelse
   endif
   fullplugfile = findfile( filepath(plugfile, root_dir=plugdir), count = ct)
   if ct ne 0 then begin
        ; If we found several plugmap files (using wildcards), take the most
        ; recent as determined by simply doing an ASCII sort of the file names.
        if (n_elements(fullplugfile) EQ 1) then fullplugfile = fullplugfile[0] $
        else fullplugfile = fullplugfile[ (reverse(sort(fullplugfile)))[0] ]
   endif else Undefine,fullplugfile
   spd1=1

   if (NOT keyword_set(fps)) then begin
     fieldid = long(config)
     fieldstr= confstr
   endif else begin
     if flavor NE 'unknown' then begin
        fieldstr=field_to_string(fieldid)
     endif else begin
        fieldid = long(config)
        fieldstr= confstr
     endelse
   endelse
   
   ;----------
   ; Construct the names of the flat and arc output files if we generate
   ; them from this exposure.

   tsetfile1 = filepath( $
      'tset-'+mjdstr+'-'+fieldstr+'-'+filee+'-'+filec+'.fits', root_dir=outdir)
   wsetfile1 = filepath( $
      'wset-'+mjdstr+'-'+fieldstr+'-'+filee+'-'+filec+'.fits', root_dir=outdir)
   fflatfile1 = filepath( $
      'fflat-'+mjdstr+'-'+fieldstr+'-'+filee+'-'+filec+'.fits', root_dir=outdir)

   

   find_cal, fieldstr, mjdstr, filee, filec, tsetfile_last = tsetfile_last, $
             wsetfile_last = wsetfile_last, fflatfile_last = fflatfile_last, $
             outdir = outdir

   splog, 'TSETFILE = ' + tsetfile_last
   splog, 'WSETFILE = ' + wsetfile_last
   splog, 'FFLATFILE = ' + fflatfile_last

   plugexist = keyword_set(fullplugfile)
   flatexist = keyword_set(tsetfile_last) AND keyword_set( findfile(tsetfile_last) )
   arcexist = keyword_set(wsetfile_last) AND keyword_set( findfile(wsetfile_last) )
   splog, 'PLUGEXIST = ', plugexist
   splog, 'FLATEXIST = ', flatexist
   splog, 'ARCEXIST = ', arcexist

   ;----------
   ; Reduce file depending on its flavor: bias/dark, flat, arc, or science/smear
   ; Report to log file if exposure is HARTMANN left of right

   rstruct = 0
   myflavor = flavor
   if (myflavor EQ 'smear') then myflavor = 'science'
   if (myflavor EQ 'dark') then myflavor = 'bias'
   if ((hartmann EQ 'Left') OR (hartmann EQ 'Right')) then begin
      myflavor = 'hartmann'
   endif

   case myflavor of
      'bias' : begin
         rstruct = quickbias(fullname, do_lock=do_lock)
      end

      'flat' : begin
         if (plugexist) then begin
            if keyword_set(fps) then begin
               rstruct = quicktrace(fullname, tsetfile1, plugmapfile=fullplugfile, $
                                    do_lock=do_lock, fps=fps, plugdir=outdir, noreject=noreject)
            endif else begin
               rstruct = quicktrace(fullname, tsetfile1, plugmapfile=fullplugfile, $
                                    do_lock=do_lock, fps=fps, noreject=noreject)
            endelse
         endif else begin
            if keyword_set(fps) then begin
                rstruct = quicktrace(fullname, tsetfile1, do_lock=do_lock, fps=fps, $
                                    plugdir=outdir, noreject=noreject)
            endif else begin
                splog, 'ABORT: Unable to reduce this flat exposure (need plug-map)'
            endelse
         endelse
      end

      'hartmann': begin
         splog, 'Skipping Hartmann exposure'
         if not strmatch(strtrim(sxpar(hdr, 'FLAVOR'),2), 'arc*', /fold_case) then splog, 'INFO: Skipping exposure with Hartmann ', hartmann
      end

      'arc' : begin
         if (flatexist) then begin
            rstruct = quickwave(fullname, tsetfile_last, wsetfile1, noreject=noreject,$
                   fflatfile1, lco = lco, do_lock=do_lock, nocal=nocal)
            if keyword_set(arc2trace) then begin
                setup_arc2trace, tsetfile_last, fflatfile1, fullname, indir, outdir, mjd, camnames[icam], fieldstr

            endif
         endif else begin
             splog, 'INFO: Arc exposure, waiting for flat before reducing'
         endelse
      end

      'science': begin
          exptime = sxpar(hdr, 'EXPTIME')
          outsci = filepath('sci-'+confstr+'-'+filec+'-'+filee+'.fits',root_dir=outdir)
          ;Added the keyword 'splitsky'.- vivek
          if (camnames[icam] eq 'r1') or (camnames[icam] eq 'r2')  then splitsky = 1B else splitsky = 0B
          if (flatexist AND arcexist AND exptime GE minexp) then begin
            rstruct = quickextract(tsetfile_last, wsetfile_last, $
                fflatfile_last, fullname, outsci, fullplugfile, outdir, mjd,$
                splitsky=splitsky, do_lock=do_lock,threshold=threshold,$
                sdssv_sn2=sdssv_sn2,sn2_15=sn2_15,arc2trace=arc2trace,forcea2t=forcea2t)
          endif else begin
             if (NOT keyword_set(flatexist)) then $
                splog, 'ABORT: Unable to reduce this science exposure (need flat)'
             if (NOT keyword_set(arcexist)) then $
                splog, 'ABORT: Unable to reduce this science exposure (need arc)'
             if (exptime LT minexp) then $
                splog, 'ABORT: Exposure time = ' + string(exptime) + ' < ' + string(minexp)
          endelse
       end

       else : begin
          splog, 'Unknown flavor: ', flavor
       end
   endcase

   ;----------
   ; Append to binary FITS log file a structure with info for this frame
   ; Lock the file to do this.

   if keyword_set(fullplugfile) then begin
    i = strpos(fullplugfile,'/',/reverse_search)
    if (i[0] EQ -1) then shortplugfile = fullplugfile $
        else shortplugfile = strmid(fullplugfile,i+1)
   endif
   ;----------
   ; Find WARNINGs and ABORTs from splog file.  Recast them as string
   ; arrays that are not empty (e.g., ''), or MWRFITS will fail.
   spawn, 'grep -e WARNING -e ABORT -e INFO '+splgfile, tstring
   if (keyword_set(tstring)) then begin
      tstruct = create_struct('FILENAME', filename, $
                              'MJD', mjd, $
                              'CONFIG', config, $
                              'FIELD', fieldid, $
                              'CARTID', cartid, $
                              'EXPNUM', filee, $
                              'CAMERA', camnames[icam], $
                              'TEXT', '' )
      tstruct = replicate(tstruct, n_elements(tstring))
      tstruct.text = tstring
   endif else begin
      tstruct = create_struct('FILENAME', filename, $
                              'MJD', mjd, $
                              'CONFIG', config, $
                              'FIELD', fieldid, $
                              'CARTID', cartid, $
                              'EXPNUM', filee, $
                              'CAMERA', camnames[icam], $
                              'TEXT', '' )
   endelse      
;   if keyword_set(tstruct) then print, tstruct
   if (keyword_set(rstruct)) then begin

      ; Get the time in TAI, which we convert on-the-fly to UT when needed.
      get_tai, hdr, tai_beg, tai_mid, tai_end

      ; Get telescope position, which is used to compute airmass.
      radeg = sxpar(hdr,'RADEG')
      decdeg = sxpar(hdr,'DECDEG')

      ; Get the CCD temperatures.
      ; Note that b1=01, b2=03, r1=04, r2=02
      ;cardname = (['TEMP01', 'TEMP03', 'TEMP04', 'TEMP02'])[icam]
      ; Note that b1=01, r1=02
      cardname = (['TEMP01', 'TEMP02'])[icam]
      ccdtemp = float(sxpar(hdr,cardname))

      airtemp = float(sxpar(hdr,'AIRTEMP', count=ct))
      if (ct EQ 0) then begin
         case strmid(camnames[icam],1,1) of
         '1': airtemp = float(sxpar(hdr,'MC1TEMDN', count=ct))
         '2': airtemp = float(sxpar(hdr,'MC2TEMDN', count=ct))
         endcase
      endif

      ; The following prevents a crash in MWRFITS.
      if (NOT keyword_set(shortplugfile)) then shortplugfile = ' '
      rstruct = create_struct('FILENAME', string(filename), $
                              'PLUGFILE', string(shortplugfile), $
                              'MJD', long(mjd), $
                              'CONFIG', long(config), $
                              'FIELD', long(fieldid), $
                              'CARTID', strtrim(sxpar(hdr,'CARTID'),2), $
                              'DESIGNID',strtrim(sxpar(hdr,'DESIGNID'),2),$
                              'EXPNUM', long(filee), $
                              'EXPTIME', float(sxpar(hdr, 'EXPTIME')), $
                              'FLAVOR', string(flavor), $
                              'CAMERA', string(camnames[icam]), $
                              'TAI', double(tai_mid), $
                              'AIRTEMP', float(airtemp), $
                              'CCDTEMP', float(ccdtemp), $
                              'QUALITY', string(sxpar(hdr,'QUALITY')), $
                              'NAME', string(sxpar(hdr,'NAME')), $
                              'OBSCOMM', string(sxpar(hdr,'OBSCOMM')), $
                              'RADEG', float(radeg), $
                              'DECDEG', float(decdeg), $
                              'AIRMASS', float(tai2airmass(radeg,decdeg,tai=tai_mid,site=getenv('OBSERVATORY'))), $
                              rstruct )
   endif

   if (keyword_set(tstruct) OR keyword_set(rstruct)) then begin
      splog, 'Appending to FITS log file '+logfile
      sos_appendlog, logfile, rstruct, tstruct
      splog, 'Done with append'
   endif

   ;----------
   ; After being passed any 'r2' frame, and if it was reduced,
   ; we re-generate the HTML file for all plates and the S/N plot for
   ; this plate.
   ; Optionally copy it to the directory specified by COPYDIR.
   ; Make an additional copy of the HTML file called 'logsheet-current.html'.

   ; Instead, create the HTML file after any reduced frame.
   if (keyword_set(rstruct) OR keyword_set(tstruct)) then begin

      if (myflavor EQ 'science') then begin
         ; Generate the added S/N^2 for this one exposure only
         plotfile1 = filepath('snplot-'+mjdstr+'-'+confstr+'-'+filee+'.ps', root_dir=outdir)
         jpegfiletmp1 = filepath('snplot-'+mjdstr+'-'+confstr+'-'+filee+'-'+filec+'.jpeg', root_dir=outdir)
         jpegfile1 = filepath('snplot-'+mjdstr+'-'+confstr+'-'+filee+'.jpeg', root_dir=outdir)
         splog, 'Generating S/N plot '+plotfile1
         sos_plotsn, logfile, config, expnum=long(filee), plugdir=plugdir, plotfile=plotfile1, fps=fps, ccd=string(camnames[icam])
         cmd = '/usr/bin/convert '+plotfile1+' '+jpegfiletmp1+' ; \mv '+jpegfiletmp1+' '+jpegfile1+' &'
         splog, 'SPAWN '+cmd, sh_out, sh_err
         spawn, cmd
         splog, 'SPAWN out=', sh_out
         splog, 'SPAWN err=', sh_err
         splog, 'Done generating plot'

         ; Generate the added S/N^2 for all exposures on this plate
         plotfile = filepath('snplot-'+mjdstr+'-'+confstr+'.ps', root_dir=outdir)
         jpegfile = filepath('snplot-'+mjdstr+'-'+confstr+'.jpeg', root_dir=outdir)
         jpegfiletmp = filepath('snplot-'+mjdstr+'-'+confstr+'-'+filec+'.jpeg', root_dir=outdir)
         splog, 'Generating S/N plot '+plotfile
         sos_plotsn, logfile, config, plugdir=plugdir, plotfile=plotfile, fps=fps, ccd=string(camnames[icam])
         cmd = '/usr/bin/convert '+plotfile+' '+jpegfiletmp+' ; \mv '+jpegfiletmp+' '+jpegfile+' &'
         splog, 'SPAWN '+cmd, sh_out, sh_err
         spawn, cmd
         splog, 'SPAWN out=', sh_out
         splog, 'SPAWN err=', sh_err
         splog, 'Done generating plot'
         
         if keyword_set(sdssv_sn2) then begin
             ; Generate the added S/N^2 for this one exposure only
             plotfile1_v2 = filepath('snplot-sdssv-'+mjdstr+'-'+confstr+'-'+filee+'.ps', root_dir=outdir)
             jpegfiletmp1_v2 = filepath('snplot-sdssv-'+mjdstr+'-'+confstr+'-'+filee+'-'+filec+'.jpeg', root_dir=outdir)
             jpegfile1_v2 = filepath('snplot-sdssv-'+mjdstr+'-'+confstr+'-'+filee+'.jpeg', root_dir=outdir)
             splog, 'Generating SDSS-V S/N plot '+plotfile1_v2
             sos_plotsn, logfile, config, expnum=long(filee), plugdir=plugdir,$
                         plotfile=plotfile1_v2, fps=fps,sdssv_sn2=sdssv_sn2, ccd=string(camnames[icam])
             cmd = '/usr/bin/convert '+plotfile1_v2+' '+jpegfiletmp1_v2+' ; \mv '+jpegfiletmp1_v2+' '+jpegfile1_v2+' &'
             splog, 'SPAWN '+cmd, sh_out, sh_err
             spawn, cmd
             splog, 'SPAWN out=', sh_out
             splog, 'SPAWN err=', sh_err
             splog, 'Done generating plot'

             ; Generate the added S/N^2 for all exposures on this plate
             plotfile_v2 = filepath('snplot-sdssv-'+mjdstr+'-'+confstr+'.ps', root_dir=outdir)
             jpegfile_v2 = filepath('snplot-sdssv-'+mjdstr+'-'+confstr+'.jpeg', root_dir=outdir)
             jpegfiletmp_v2 = filepath('snplot-sdssv-'+mjdstr+'-'+confstr+'-'+filec+'.jpeg', root_dir=outdir)
             splog, 'Generating  SDSS-V S/N plot '+plotfile_v2
             sos_plotsn, logfile, config, plugdir=plugdir, plotfile=plotfile_v2, fps=fps,sdssv_sn2=sdssv_sn2, ccd=string(camnames[icam])
             cmd = '/usr/bin/convert '+plotfile_v2+' '+jpegfiletmp_v2+' ; \mv '+jpegfiletmp_v2+' '+jpegfile_v2+' &'
             splog, 'SPAWN '+cmd, sh_out, sh_err
             spawn, cmd
             splog, 'SPAWN out=', sh_out
             splog, 'SPAWN err=', sh_err
             splog, 'Done generating SDSS-V plot'
         endif
         if keyword_set(sn2_15) then begin
             ; Generate the added S/N^2 for this one exposure only
             plotfile1 = filepath('snplot-sdssv15-'+mjdstr+'-'+confstr+'-'+filee+'.ps', root_dir=outdir)
             jpegfiletmp1 = filepath('snplot-sdssv15-'+mjdstr+'-'+confstr+'-'+filee+'-'+filec+'.jpeg', root_dir=outdir)
             jpegfile1 = filepath('snplot-sdssv15-'+mjdstr+'-'+confstr+'-'+filee+'.jpeg', root_dir=outdir)
             splog, 'Generating SDSS-V S/N plot '+plotfile1
             sos_plotsn, logfile, config, expnum=long(filee), plugdir=plugdir,$
                         plotfile=plotfile1, fps=fps,sn2_15=sn2_15, ccd=string(camnames[icam])
             cmd = '/usr/bin/convert '+plotfile1+' '+jpegfiletmp1+' ; \mv '+jpegfiletmp1+' '+jpegfile1+' &'
             splog, 'SPAWN '+cmd, sh_out, sh_err
             spawn, cmd
             splog, 'SPAWN out=', sh_out
             splog, 'SPAWN err=', sh_err
             splog, 'Done generating plot'

             ; Generate the added S/N^2 for all exposures on this plate
             plotfile = filepath('snplot-sdssv15-'+mjdstr+'-'+confstr+'.ps', root_dir=outdir)
             jpegfile = filepath('snplot-sdssv15-'+mjdstr+'-'+confstr+'.jpeg', root_dir=outdir)
             jpegfiletmp = filepath('snplot-sdssv15-'+mjdstr+'-'+confstr+'-'+filec+'.jpeg', root_dir=outdir)
             splog, 'Generating  SDSS-V Mag 15 S/N plot '+plotfile
             sos_plotsn, logfile, config, plugdir=plugdir, plotfile=plotfile, fps=fps,sn2_15=sn2_15, ccd=string(camnames[icam])
             cmd = '/usr/bin/convert '+plotfile+' '+jpegfiletmp+' ; \mv '+jpegfiletmp+' '+jpegfile+' &'
             splog, 'SPAWN '+cmd, sh_out, sh_err
             spawn, cmd
             splog, 'SPAWN out=', sh_out
             splog, 'SPAWN err=', sh_err
             splog, 'Done generating SDSS-V Mag 15 plot'
         endif
      endif

      splog, 'Generating HTML file '+htmlfile
      sos_log2html, logfile, htmlfile, fps=fps, sn2_15=sn2_15;, sdssv_sn2=sdssv_sn2
      splog, 'Done generating HTML file'

      ; Generate a copy of the HTML file, 'logsheet-current.html',
      ; that includes the Java script to auto-load the page every 60 seconds.

      squote = "\'"
      addstring = $
       '<BODY ONLOAD=\"timerID=setTimeout(' $
       +squote+'location.reload(true)'+squote+',60000)\">'
      sedcommand = '-e "s/<\/HEAD>/<\/HEAD>'+addstring+'/g"'
      sedcommand = sedcommand + ' -e "s/BOSS Spectro/BOSS Spectro (Current)/g"'
      setenv, 'SHELL=bash'
      spawn, 'sed ' + sedcommand + ' ' + htmlfile + ' > ' + currentfile

      if (keyword_set(copydir)) then begin
         FILE_MKDIR, copydir
         splog, 'Copying files to ', copydir
         spawn, 'scp ' + htmlfile + ' ' + copydir
         spawn, 'scp ' + currentfile + ' ' + copydir
         htmlfile_c = djs_filepath(file_basename(htmlfile), root_dir=copydir)
         currentfile_c = djs_filepath(file_basename(currentfile), root_dir=copydir)

         yesterday = strtrim((long(mjd)-1),2)
         sedcommand = ' -e "s/Yesterday: <A HREF=..\/'+yesterday+'\//Yesterday: <A HREF=/g"'
         tomorrow = strtrim((long(mjd)+1),2)
         sedcommand = sedcommand + ' -e "s/Tomorrow: <A HREF=..\/'+tomorrow+'\//Tomorrow: <A HREF=/g"'
         spawn, 'sed ' + sedcommand + ' ' + htmlfile + ' > ' + htmlfile_c
         spawn, 'sed ' + sedcommand + ' ' + currentfile + ' > ' + currentfile_c

         
;         spawn, 'scp ' + logfile  + ' ' + copydir
;         if (keyword_set(plotfile)) then $
;          spawn, 'scp ' + plotfile + ' ' + plotfile1 $
;           + ' ' + jpegfile + ' ' + jpegfile1 + ' ' + copydir
;         if (keyword_set(plotfile_v2)) then $
;          spawn, 'scp ' + plotfile_v2 + ' ' + plotfile1_v2 $
;           + ' ' + jpegfile_v2 + ' ' + jpegfile1_v2 + ' ' + copydir
         splog, 'Done.'
      endif
   endif

   ;----------
   ; Close splog file

   splog, 'Elapsed time = ', systime(1)-t0
   splog, 'Finished at ', systime()
   splog, /close
   ;asdsa
   return
end
;------------------------------------------------------------------------------
