;+
; NAME:
;   apofix
;
; PURPOSE:
;   Add line to sdHdrFix file to denote change in FITS header for sdR files.
;
; CALLING SEQUENCE:
;   apofix, expnum, [ card, value, camera=, copydir=, /bad, /test, /excellent, /not_sos ]
;
; INPUTS:
;   expnum     - Exposure number
;
; OPTIONAL INPUTS:
;   card       - FITS header keyword to change; this is case-insensitive,
;                so that 'exptime' is the same as 'EXPTIME'.
;   value      - New value for FITS header keyword.
;   camera     - Camera name in which to change values, e.g. 'b1', 'r1',
;                'b2' or 'r2'.  A '?' can be used as a wildcard, for example
;                '?2' to denote a change to both 'b2' and 'r2'.  Default to
;                '??' to denote a change to sdR files for all 4 cameras.
;   copydir    - If set, then copy the output log files to this directory, default is
;                to copy to /data/boss/sos/combined
;   bad        - If set, then declare the specified exposure number to be bad.
;                This is equivalent to setting QUALITY='bad'.
;   test       - If set, then declare the specified exposure number to be test.
;                This is equivalent to setting QUALITY='test'.
;   excellent  - If set, then declare the specified exposure number to be excellent.
;                This is equivalent to setting QUALITY='excellent'.
;   not_sos    - This keyword can be set to run this proc on a machine
;                that is not named "sos" or "sos3".  This would only be done for
;                testing purposes, or if Son-of-Spectro has been moved
;                to another machine.
;   sos_dir    - This keyword sets the location of the sos output files
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only the following keywords can be modified with this procedure:
;     CAMERAS, FLAVOR, PLATEID, NAME, EXPTIME, TAI-BEG, TAI-END, TAI,
;     FFS, FF, NE, HGCD, OBSCOMM, QUALITY
;   The AIRMASS is not read from the header, but computed from RADEG,DECDEG
;   and the TAI-BEG,TAI-END keywords.
;   Refer to the Son-of-Spectro documentation for the specifics of
;   valid values for each keyword.
;
;   Note that string values must be enclosed in single- or double-quotes.
;   Numeric values should not be in quotes.  Double-precision numbers should
;   be written with "d" notation, for example 3.14d7 instead of 3.14e7.
;
; EXAMPLES:
;   Fix the exposure time for exposure #1234 to be 900 sec:
;     IDL> apofix, 1234, 'exptime', 900
;
;   Fix the TAI time for exposure #1234 to be 4.443852968d+09
;   (use the "d" notation for double-precision, even though it
;   will appear in the sdHdrFix file with an "e"):
;     IDL> apofix, 1234, 'TAI', 4.443852968d+09
;
;   Declare exposure #1234 as bad:
;     IDL> apofix, 1234, /bad
;   or equivalently:
;     IDL> apofix, 1234, 'quality', 'bad'
;
;   Only declare the 'b1' camera bad for exposure number 1234:
;     IDL> apofix, 1234, /bad, camera='b1'
;
;   The wrong NAME is in the header, which is necessary to identify
;   the proper plug-map file.  If the correct plug-map file
;   is plPlugMapM-0328-52277-01, then edit as follows:
;     IDL> apofix, 1234, 'name', '0328-52277-01'
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_lockfile()
;   djs_modfits
;   djs_unlockfile
;   fileandpath()
;   fits_wait
;   headfits()
;   struct_append
;   sxpar()
;   yanny_read
;
; REVISION HISTORY:
;   22-Apr-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro apofix, expnum, card, newval, copydir=copydir, camera=camera, bad=bad, test=test, $
 excellent=excellent, not_sos=not_sos, sos_dir=sos_dir

   common apofix_com, apo_uname

   if (n_params() LT 1) then begin
      doc_library, 'apofix'
      return
   end

   !quiet = 1

   ;----------
   ; Insist that this proc only run on machines named "sos"

   if (NOT keyword_set(apo_uname)) then begin
      spawn, 'uname -n', uname_string
      apo_uname = (strsplit(uname_string[0], '.', /extract))[0]
   endif

   if (apo_uname NE 'sos' AND apo_uname NE 'sos3' AND keyword_set(not_sos) EQ 0) then begin
      print, 'This procedure can only be run on the machine sos.apo.nmsu.edu'
      print, 'If Son-of-Spectro is now running on a different machine, you'
      print, 'can call this routine with /not_sos'
      return
   endif

   if keyword_set(not_sos) eq 1 then begin
      print,'You are not running this machine on sos or sos3'
      print,'Make sure the svn repository has been updated before committing any changes'
   endif

   ;----------
   ; Set input directory for sdR files and sdHdrFix file

   rawdata_dir = getenv('BOSS_SPECTRO_DATA')
   if (NOT keyword_set(rawdata_dir)) then $
    rawdata_dir = '/data/spectro'

   speclog_dir = getenv('SPECLOG_DIR')
   if (NOT keyword_set(speclog_dir)) then $
    speclog_dir = '/home/sdss3/products/NULL/speclog/trunk'

   if (NOT keyword_set(copydir)) then $
    copydir = '/data/boss/sos/combined/'


   ;----------
   ; Create output structure

   fstruct = create_struct( $
    name = 'OPHDRFIX', $
    'fileroot', '', $
    'keyword' , '', $
    'value'   , '' )

   ;----------
   ; Sanity checks on EXPNUM, CAMERA, CARD, VALUE

   expnum = long(expnum)
   if (expnum LE 0 OR expnum GT 99999999L OR n_elements(expnum) NE 1) then begin
      print, 'EXPNUM must be a number between 1 and 99999999'
      return
   endif

   if (NOT keyword_set(camera)) then camera = '??'
   c1 = strmid(camera,0,1)
   c2 = strmid(camera,1,1)
   if (size(camera, /tname) NE 'STRING' $
    OR strlen(camera) NE 2 $
    OR n_elements(camera) NE 1 $
    OR (c1 NE '?' AND c1 NE 'b' AND c1 NE 'r') $
    OR (c2 NE '?' AND c2 NE '1' AND c2 NE '2') ) then begin
      print, 'CAMERA must be a 2-character string'
      return
   endif

   if (keyword_set(bad) AND keyword_set(test)) then begin
      print, 'Invalid to set both the /BAD and /TEST keywords. Aborting!'
      return
   endif

   if (keyword_set(bad) AND keyword_set(excellent)) then begin
      print, 'Invalid to set both the /BAD and /EXCELLENT keywords. Aborting!'
      return
   endif

   if (keyword_set(excellent) AND keyword_set(test)) then begin
      print, 'Invalid to set both the /EXCELLENT and /TEST keywords. Aborting!'
      return
   endif

   if (keyword_set(bad)) then begin
      if (keyword_set(card) OR keyword_set(value)) then begin
         print, 'Invalid to set the /BAD flag along with CARD or VALUE. Aborting!'
         return
      endif
      card = 'QUALITY'
      newval = 'bad'
   endif

   if (keyword_set(test)) then begin
      if (keyword_set(card) OR keyword_set(value)) then begin
         print, 'Invalid to set the /TEST flag along with CARD or VALUE. Aborting!'
         return
      endif
      card = 'QUALITY'
      newval = 'test'
   endif

   if (keyword_set(excellent)) then begin
      if (keyword_set(card) OR keyword_set(value)) then begin
         print, 'Invalid to set the /EXCELLENT flag along with CARD or VALUE. Aborting!'
         return
      endif
      card = 'QUALITY'
      newval = 'excellent'
   endif

   if (size(card, /tname) NE 'STRING' $
    OR n_elements(card) NE 1 $
    OR strlen(card) EQ 0 OR strlen(card) GT 8) then begin
      print, 'CARD must be a string of 1 to 8 characters'
      return
   endif

   if (n_elements(newval) NE 1) then begin
      print, 'VALUE must be specified (and a scalar)'
      return
   endif
   if (size(newval, /tname) EQ 'DOUBLE') then format='(e17.10)' $
    else format=''
   strval = strtrim(string(newval,format=format),2)
   if (strpos(strval,'"') NE -1 OR strpos(strval,"'") NE -1) then begin
      print, 'VALUE cannot contain single or double-quotes'
      return
   endif
   ; Enclose any string in single-quotes, so that the SPHDRFIX routine
   ; can desriminate between strings and numbers.
   if (size(newval, /tname) EQ 'STRING') then strval = "'" + strval + "'"

   ;----------
   ; Explicitly test the values for each possible keyword name.
   ; This is a different test for each keyword.

   qstring = size(newval, /tname) EQ 'STRING'
   case strupcase(card) of
   'EXPOSURE': begin
      print, 'The EXPOSURE number is always set to that in the file name. Aborting!'
      return
      end
   'CAMERAS': begin
      possible = ['b1','b2','r1','r2']
      if (total(newval EQ possible) EQ 0) then begin
         print, 'Valid values = ', "'"+possible+"'"
         print, 'Invalid value for CAMERAS. Aborting!'
         return
      endif
      end
   'FLAVOR': begin
      if (newval EQ 'unknown') then begin
         print, 'Please set QUALITY to bad instead of setting FLAVOR=unknown. Aborting!'
         return
      endif
      possible = ['bias','dark','flat','arc','science','smear']
      if (total(newval EQ possible) EQ 0) then begin
         print, 'Valid values = ', "'"+possible+"'"
         print, 'Invalid value for CAMERAS. Aborting!'
         return
      endif
      end
   'MJD': begin
      print, 'Not possible to change MJD. This is simply big trouble. Aborting!'
      return
      end
   'PLATEID': begin
      if (long(newval) LE 0 OR long(newval) GT 9999) then begin
         print, 'PLATEID must be between 1 and 9999. Quitting!'
         return
      endif
      end
   'NAME': begin
      if (strmatch(newval,'[0-9][0-9][0-9][0-9]-[5-9][0-9][0-9][0-9][0-9]-[0-9][0-9]')) then begin
         print, 'NAME must be of the form "????-?????-??" [all digits]. Aborting!'
         return
      endif
      end
   'EXPTIME': begin
      if (long(newval) LT 0 OR long(newval) GT 3600 OR qstring) then begin
         print, 'Valid exposure times must be between 0 and 3600 sec. Aborting!'
         return
      endif
      end
   'TAI-BEG': begin
      if (double(newval) LT 4d9 OR double(newval) GT 6d9) then begin
         print, 'Valid times are between 4d9 and 6d9. Aborting!'
         return
      endif
      end
   'TAI-END': begin
      if (double(newval) LT 4d9 OR double(newval) GT 6d9) then begin
         print, 'Valid times are between 4d9 and 6d9. Aborting!'
         return
      endif
      end
   'TAI': begin
      if (double(newval) LT 4d9 OR double(newval) GT 6d9) then begin
         print, 'Valid times are between 4d9 and 6d9. Aborting!'
         return
      endif
      end
   'FFS': begin
      if (strmatch(newval,'[01] [01] [01] [01] [01] [01] [01] [01]')) then begin
         print, 'FFS must be of the form "? ? ? ? ? ? ? ?" [all 0 or 1]. Aborting!'
         return
      endif
      end
   'FF': begin
      if (strmatch(newval,'[01] [01] [01] [01]')) then begin
         print, 'FF must be of the form "? ? ? ?" [all 0 or 1]. Aborting!'
         return
      endif
      end
   'NE': begin
      if (strmatch(newval,'[01] [01] [01] [01]')) then begin
         print, 'NE must be of the form "? ? ? ?" [all 0 or 1]. Aborting!'
         return
      endif
      end
   'HGCD': begin
      if (strmatch(newval,'[01] [01] [01] [01]')) then begin
         print, 'HGCD must be of the form "? ? ? ?" [all 0 or 1]. Aborting!'
         return
      endif
      end
   'OBSCOMM': begin
      possible = ['{dithered flats-flat}', '{dithered flats-arc}', $
       '{focus, hartmann l}', '{focus, hartmann r}']
      if (total(newval EQ possible) EQ 0) then begin
         print, 'Valid values = ', "'"+possible+"'"
         print, 'Invalid value for OBSCOMM. Aborting!'
         return
      endif
      end
   'QUALITY': begin
      possible = ['excellent','test','bad']
      if (total(newval EQ possible) EQ 0) then begin
         print, 'Valid values = ', "'"+possible+"'"
         print, 'Invalid value for QUALITY. Aborting!'
         return
      endif
      end
   'TILEID': begin
      if (long(newval) LE 0 OR long(newval) GT 9999) then begin
         print, 'TILEID must be between 1 and 9999. Quitting!'
         return
      endif
      end
   'RA': begin
      if (double(newval) LT 0d OR double(newval) GT 360d OR qstring) then begin
         print, 'Valid RA are between 0d and 360d. Aborting!'
         return
      endif
      end
   'DEC': begin
      if (double(newval) LT -90d OR double(newval) GT 90d OR qstring) then begin
         print, 'Valid DEC are between -90d and 90d. Aborting!'
         return
      endif
      end
   'RADEG': begin
      if (double(newval) LT 0d OR double(newval) GT 360d OR qstring) then begin
         print, 'Valid RADEG are between 0d and 360d. Aborting!'
         return
      endif
      end
   'DECDEG': begin
      if (double(newval) LT -90d OR double(newval) GT 90d OR qstring) then begin
         print, 'Valid DECDEG are between -90d and 90d. Aborting!'
         return
      endif
      end
   'AIRTEMP': begin
      if (float(newval) LT -40 OR double(newval) GT 40 OR qstring) then begin
         print, 'Valid AIRTEMP are between -40 and 40. Aborting!'
         return
      endif
      end
   else: begin
      print, 'This keyword is not of interest. Aborting!'
      return
      end
   endcase

   ;----------
   ; Test that sdR files exist that correspond to the exposure number
   ; and camera(s) specified.

   fileroot = string(camera, expnum, format='("sdR-",a2,"-",i8.8)')

   filename = findfile(filepath(fileroot+'.fit*', $
    root_dir=rawdata_dir, subdir='*'), count=nfile)
   if (nfile EQ 0) then begin
      print, 'File=' + fileroot + '.fit not found in dir=' + rawdata_dir
      return
   endif

   ;----------
   ; Read the values of the cards in these sdR files
   ; Also, compare values in sdR headers to requested values.

   for ifile=0, nfile-1 do begin
      qdone = fits_wait(filename[ifile], deltat=2, tmax=10, /header_only)
      if (qdone) then begin
         thishdr = headfits(filename[ifile])
         thismjd = sxpar(thishdr, 'MJD')
         oldval = sxpar(thishdr, card)
         print, strmid(fileandpath(filename[ifile]),0,15), strupcase(card), $
          strtrim(string(oldval,format=format),2), strval, $
          format='(a, 1x, a8, "=[", a, "] -> [", a, "]")'
      endif
   endfor

   ;----------
   ; Construct the output structure

   fstruct.fileroot = fileroot
   fstruct.keyword = strupcase(card)
   fstruct.value = strval

   ;----------
   ; Read the sdHdrFix file

   if (NOT keyword_set(thismjd)) then begin
      print, 'MJD could not be determined from the FITS headers'
      return
   endif

   mjdstr = string(thismjd, format='(i5.5)')
   sdfixname = filepath('sdHdrFix-' + mjdstr + '.par', root_dir=speclog_dir, $
    subdir=mjdstr)

   while (djs_lockfile(sdfixname) EQ 0) do wait, 2
   yanny_read, sdfixname, pdata, hdr=hdr, enums=enums, structs=structs, $
    stnames=stnames, /anonymous
   djs_unlockfile, sdfixname

   ;----------
   ; Use an explicitly-defined structs such that the fTCL Yanny-readers
   ; will still work.

   structs = ['typedef struct {', $
              '  char fileroot[20]; # Root of file name, without any ".fit" suffix', $
              '  char keyword[9]; # Keyword name', $
              '  char value[80]; # Keyword value (as a string)', $
              '} OPHDRFIX;']

   ;----------
   ; Append to the existing data structures in the sdHdrFix file,
   ; or create a new one if the file does not yet exist.

   if (NOT keyword_set(pdata)) then begin
      while (djs_lockfile(sdfixname) EQ 0) do wait, 2
      yanny_write, sdfixname, ptr_new(fstruct), structs=structs
      djs_unlockfile, sdfixname
      ncorr = 1
   endif else begin
      ; Append data to the relevant data structure.
      i = (where(stnames EQ tag_names(fstruct, /structure_name)))[0]
      if (i[0] EQ -1) then begin
         print, 'The file ' + sdfixname + ' appears to be invalid'
         return
      endif
      pdata[i] = ptr_new(struct_append(*pdata[i], fstruct))

      while (djs_lockfile(sdfixname) EQ 0) do wait, 2
      yanny_write, sdfixname, pdata, hdr=hdr, enums=enums, structs=structs, $
       stnames=stnames
      djs_unlockfile, sdfixname
      ncorr = n_elements(*pdata[i])
   endelse

   yanny_free, pdata

   print, 'File ' + fileandpath(sdfixname) + ' contains ' $
    + strtrim(string(ncorr),2) + ' declared changes.'

   ;----------
   ; If QUALITY keyword is changed, and the routine was run from sos or sos3, then edit the APO logfile.

   if (strupcase(card) EQ 'QUALITY' AND keyword_set(not_sos) eq 0) then begin

      if (NOT keyword_set(sos_dir)) then $
      sos_dir = '/data/boss/sos/'
      log_root = 'logfile-' + mjdstr + '.fits'
;      logfile = filepath(logfile, root_dir=sos_dir, subdir=mjdstr)
      logfile=findfile(filepath(log_root+'*',root_dir=sos_dir,subdir=mjdstr))
print,'logfile = ',logfile
      if (NOT keyword_set(findfile(logfile))) then begin
         splog, 'Unable to find logfile '+logfile
         return
      endif
      htmlfile = filepath('logfile-' + mjdstr + '.html', root_dir=sos_dir, subdir=mjdstr)
      currentfile = filepath('logfile-current.html', root_dir=sos_dir, subdir=mjdstr)

      splog, 'Trying to lock the logfile  ' + logfile
      while(djs_lockfile(logfile) EQ 0) do wait, 1

      ;----------
      ; Loop through each of the first 4 HDUs in the log file.

      splog, 'Reading the logfile ' + logfile
      for thishdu=1, 4 do begin
         rstruct = mrdfits(logfile, thishdu, /silent)
         nstruct = n_elements(rstruct) * (keyword_set(rstruct))

         qchange = 0
         for i=0, nstruct-1 do begin
            if (expnum EQ rstruct[i].expnum $
             AND strmatch(rstruct[i].camera,camera)) then begin
               rstruct[i].quality = newval
               qchange = 1
            endif
         endfor

         if (qchange) then begin
            djs_modfits, logfile, rstruct, exten_no=thishdu
         endif
      endfor

      djs_unlockfile, logfile

      junk='stringplaceholder'
      while(djs_lockfile(copydir+'logfile-current.html') EQ 0) do wait, 1
      openr,101,copydir+'logfile-current.html'
      readf,101,junk,format='(a)'
      readf,101,junk,current_mjd,format='(a40,i5)'
      close,101
      if current_mjd le thismjd then begin
         ; Now convert log fits file into html webpage and update plateDB accordingly
         apo_log2html, logfile, htmlfile
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
   
         ; now copy files over to directory that is displayed on the web browser
         spawn, 'cp ' + htmlfile + ' ' + copydir
         spawn, 'cp ' + currentfile + ' ' + copydir
      endif

      spawn, 'cp ' + logfile  + ' ' + copydir
      djs_unlockfile,copydir+'logfile-current.html'

      name = strtrim(sxpar(thishdr,'NAME'),2)
      plateid = strtrim(sxpar(thishdr,'PLATEID'),2)

      fileroot = string(plateid,camera, expnum, format='("sci-",a4,"-",a2,"-",i8.8)')
      sciname = findfile(filepath(fileroot+'.fit*', root_dir=sos_dir, subdir=mjdstr), count=nscience)

      ; This string should contain PLATE-MJD-PLUGID, but it may not
      ; in some of the early data, in which case we're search using wildcards
      if (strlen(name) LT 13) then name = '*' + name + '*'
      plugfile = 'plPlugMapM-'+name+'.par'
      fullplugfile = findfile( filepath(plugfile, root_dir=speclog_dir, subdir=mjdstr) )
 
      for i=0,nscience-1 do begin
         while(djs_lockfile(sciname[i]) EQ 0) do wait, 1
         rstruct = mrdfits(sciname[i],0, /silent)
         hdr=headfits(sciname[i])
         sxaddpar,hdr,'QUALITY',newval
         djs_modfits,sciname[i],rstruct,hdr
         djs_unlockfile, sciname[i]
         spawn,'loadSN2Value --update -v ' + sciname[i] + ' ' + fullplugfile
      endfor
      spawn,'svn add ' + sdfixname

   endif else begin
      print, 'Changes will not appear unless you run sos_redo.'
      print, 'Equivalently, you could type the following:'
      print,'sos_aporedo -g -m '+string(thismjd)+' -z '+string(expnum)
;      print, '  remove2redo, mjd=' + mjdstr $
;       + ', expnum=' + strtrim(string(expnum),2)
      print, ''
   endelse

   if keyword_set(not_sos) eq 1 then begin
      print,'You are not running this machine from APO'
      print,'You must commit the changes you have made to the hdrfixfile manually in order'
      print,'to have the flags appear in the svn repository'
   endif


   return
end
;------------------------------------------------------------------------------
