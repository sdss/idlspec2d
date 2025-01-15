;+
; NAME:
;   sos_log2html
;
; PURPOSE:
;   Convert output FITS file from APOREDUCE to HTML format.
;
; CALLING SEQUENCE:
;   sos_log2html, logfile, [ htmlfile ]
;
; INPUTS:
;   logfile    - Input log file as a FITS binary file with an extension
;                for each frame reduced; this file is written by APOREDUCE.
;
; OPTIONAL INPUTS:
;   htmlfile   - Output log file in HTML format; default to the same name
;                as LOGFILE, replacing the '.fits' extension with '.html'
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   sos_checklimits()
;   copy_struct_inx
;   djs_filepath()
;   djs_findfile()
;   djs_lockfile()
;   djs_unlockfile
;   fileandpath()
;   headfits()
;   mrdfits
;   repstr()
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   sos_color2hex()
;   sos_log_header()
;   sos_log_endfile()
;   sos_log_tableline()
;   sos_log_beginplate()
;   sos_log_endplate()
;   sos_log_fields()
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
function sos_log_header, title1

   ; Include a Java script to auto-load this page every 60 seconds
   ; --> No.  Disable this now.  A copy of the file is made from APOREDUCE
   ; that includes this Java script.

   textout = '<HTML>'
   textout = [textout, '<HEAD><TITLE>' + title1 + '</TITLE></HEAD>']
   squote = "'"
;   textout = [textout, $
;    '<BODY ONLOAD="timerID=setTimeout('+squote+'location.reload(true)'+squote+',60000)">']

   return, textout
end

;------------------------------------------------------------------------------
function sos_log_endfile, title

   textout = '</BODY></HTML>'

   return, textout
end

;------------------------------------------------------------------------------
function sos_log_tableline, ncams

   textout = ''

; The following would generate a bunch of horizontal rules in each field...
;   rowsep = ' <TR> <TD> <HR> '
;   colsep = ' <TD> <HR> '
;
;   textout = rowsep + colsep + colsep
;   for icam=0, ncams-1 do $
;    textout = textout + colsep

   return, textout
end

;------------------------------------------------------------------------------
function sos_log_beginplate, platenum, cartid, mjd, fieldid, camnames, designMode, outdir=outdir, fps=fps

   rowsep = ' <TR> <TH> '
   colsep = ' <TH> '

   if (NOT keyword_set(fps)) then begin
      var_str='Plate'
      var_str1='PLATE'
   endif else begin
      var_str='Configuration'
      var_str1='CONFIG'
   endelse
   
   ncams = n_elements(camnames)

   mjdstr = strtrim(string(mjd),2)
   platestr = strtrim(string(platenum),2)+" ("+strtrim(designMode,2)+")"
   cartstr = strtrim(string(cartid),2)
   fieldstr = strtrim(string(fieldid),2)
   ;platestr4 = plate_to_string(platenum)
   platestr4 = config_to_string(platenum)
   plotfile = 'snplot-'+mjdstr+'-'+platestr4+'.ps'
   jpegfile = 'snplot-'+mjdstr+'-'+platestr4+'.jpeg'

   textout = ['<A NAME="'+var_str1 + platestr + '">']
   textout = [textout, '<TABLE BORDER=1 CELLPADDING=3>']
   textout = [textout, sos_log_tableline(ncams)]
   
   if keyword_set(fps) then begin
     nextline = '<CAPTION><FONT SIZE="+2"><B> '+var_str+' ' + platestr $
      + ' on Cart ' + cartstr + ' on Field #' + fieldstr + '</B></FONT>'
   endif else begin
     nextline = '<CAPTION><FONT SIZE="+2"><B> '+var_str+' ' + platestr $
      + ' on Cart ' + cartstr + '</B></FONT>' 
   endelse
   textout = [textout, nextline]
   nextline = rowsep + colsep
   for icam=0, ncams-1 do $
    nextline = nextline + colsep + camnames[icam]
   nextline = nextline + colsep + 'EXPTIME' + colsep + 'TEMP' $
    + colsep + 'UT' + colsep + 'QUALITY'
   textout = [textout, nextline]

   textout = [textout, sos_log_tableline(ncams)]

   return, textout
end

;------------------------------------------------------------------------------
function sos_log_endplate

   textout = ['</TABLE>']

   return, textout
end

;------------------------------------------------------------------------------
function sos_log_fields, pp, fields, printnames=printnames, formats=formats

   common com_sos_log, camnames

   rowsep = ' <TR> <TH> '
   colsep = ' <TD ALIGN=RIGHT> '

   ncams = n_elements(pp)
   igood = where(strtrim(pp.flavor,2) NE '')
   if (igood[0] EQ -1) then return, ''

   flavor = strtrim(pp[igood[0]].flavor,2)

   ; Print the exposure number as long as EXPNUM is set, which is always
   ; case except for the TOTAL S/N^2 row.
   if (keyword_set(pp[igood[0]].expnum)) then begin
      expstring = string(pp[igood[0]].expnum, format='(i8.8)')
   endif else begin
      expstring = ''
   endelse

   ; Print the UT time as long as TAI is set, which is always the
   ; case except for the TOTAL S/N^2 row.
   ; Note that we always get EXPTIME,AIRTEMP,UT from the first camera listed
   ; in the pp structure.
   if (pp[igood[0]].tai NE 0) then begin
      jd = 2400000.5D + pp[igood[0]].tai / (24.D*3600.D)
      caldat, jd, jd_month, jd_day, jd_year, jd_hr, jd_min, jd_sec
      utstring = string(jd_hr, jd_min, format='(i2.2,":",i2.2," Z")')

      airtempstring = string(pp[igood[0]].airtemp, format='(f6.1)')
      exptimestring = sos_checklimits(pp[igood[0]].flavor, 'EXPTIME', $
       pp[igood[0]].camera, pp[igood[0]].exptime, /html) $
       + string(pp[igood[0]].exptime, format='(f8.1)')
      qualstring = sos_checklimits(pp[igood[0]].flavor, 'QUALITY', $
       pp[igood[0]].camera, pp[igood[0]].quality, /html) $
       + pp[igood[0]].quality
   endif else begin
      utstring = ''
      airtempstring = ''
      exptimestring = ''
      qualstring = ''
   endelse

   tags = tag_names(pp[igood[0]])
   for ifield=0, n_elements(fields)-1 do begin
      itag = (where(fields[ifield] EQ tags))[0]
      if (keyword_set(printnames)) then nextline = colsep + printnames[ifield] $
       else nextline = colsep + fields[ifield]
      if (keyword_set(formats)) then format = formats[ifield]
      for icam=0, ncams-1 do begin
         value = ' '
         if (keyword_set(pp[icam].flavor)) then begin
            tmpval = pp[icam].(itag)
            if (keyword_set(tmpval)) then $
             value = string(tmpval, format=format)
            if strcmp(strtrim(value,2),'NaN') then value = '-'
            value = sos_checklimits(flavor, fields[ifield], $
             camnames[icam], tmpval, /html) + value
         endif
         nextline = nextline + colsep + value
      endfor
      if (ifield EQ 0) then begin
       nextline = nextline + colsep
       textout = rowsep + strupcase(flavor) + '-' + expstring $
                + nextline + exptimestring + colsep $
                + airtempstring + colsep + utstring + colsep + qualstring
      endif else $
        textout = [textout, rowsep + nextline]
   endfor

   textout = [textout, sos_log_tableline(ncams)]

   return, textout
end

;------------------------------------------------------------------------------
pro sos_log2html, logfile, htmlfile, fps=fps, sdssv_sn2=sdssv_sn2, obs=obs, sn2_15=sn2_15, brightsn2=brightsn2
    
   common com_sos_log, camnames

   if (n_params() EQ 0) then begin
      print, 'Syntax: sos_log2html, logfile, [ htmlfile ]'
      return
   endif else if (n_params() EQ 1) then begin
      thisfile = fileandpath(logfile, path=thispath)
      ipos = strpos(thisfile, '.fits', /reverse_search)
      if (ipos EQ -1) then ipos = strlen(thisfile)
      htmlfile = djs_filepath(strmid(thisfile, 0, ipos) + '.html', $
       root_dir=thispath)
   endif

   if (NOT keyword_set(fps)) then begin
      var_str='Plate'
      var_str1='PLATE'
   endif else begin
      var_str='Configuration'
      var_str1='CONFIG'
   endelse
   
   junk = fileandpath(htmlfile, path=outdir)
   ;camnames = ['b1', 'r1', 'b2', 'r2']
   if not keyword_set(obs) then obs = strtrim(STRLOWCASE(getenv('OBSERVATORY')),2)
   if strmatch(obs, 'apo',/fold_case) eq 1 then begin
        camnames = ['b1', 'r1']
   endif else begin
        camnames = ['b2', 'r2']
   endelse
   ncams = n_elements(camnames)

   ; Lock the files.
   ;print, htmlfile
   ;print, logfile
   while(djs_lockfile(htmlfile, lun=html_lun) EQ 0) do wait, 5
   while(djs_lockfile(logfile) EQ 0) do wait, 5

   ; Read the 0th header to get the version of the code
   hdr = headfits(logfile)
   vers2d = sxpar(hdr, 'VERS2D')
   run2d = sxpar(hdr,'RUN2D')

   ; Read in all the HDU's in the log file as structures
   PPBIAS = mrdfits(logfile, 1)
   PPFLAT = mrdfits(logfile, 2)
   PPARC = mrdfits(logfile, 3)
   PPSCIENCE = mrdfits(logfile, 4)
   PPTEXT = mrdfits(logfile, 5)
   djs_unlockfile, logfile
   if (NOT keyword_set(PPBIAS) AND NOT keyword_set(PPFLAT) $
    AND NOT keyword_set(PPTEXT)) AND NOT keyword_set(PPSCIENCE) then begin
      djs_unlockfile, htmlfile, lun=html_lun
      return
   endif

   allplates = [0]
   allcarts = ['']
   allModes = ['']
   if keyword_set(fps) then allfields = [0]
   if (keyword_set(PPBIAS)) then begin
      allplates = [allplates, PPBIAS.config]
      allcarts = [allcarts, PPBIAS.cartid]
      allModes = [allModes, PPBIAS.allModes]
      if keyword_set(fps) then allfields = [allfields, PPBIAS.field]
      thismjd = PPBIAS[0].mjd
   endif
   if (keyword_set(PPFLAT)) then begin
      allplates = [allplates, PPFLAT.config]
      allcarts = [allcarts, PPFLAT.cartid]
      allModes = [allModes, PPFLAT.allModes]
      if keyword_set(fps) then allfields = [allfields, PPFLAT.field]
      thismjd = PPFLAT[0].mjd
   endif
   if (keyword_set(PPTEXT)) then begin
      allplates = [allplates, PPTEXT.config]
      allcarts = [allcarts, PPTEXT.cartid]
      allModes = [allModes, PPTEXT.allModes]
      if keyword_set(fps) then allfields = [allfields, PPTEXT.field]
      thismjd = PPTEXT[0].mjd
   endif
   if (keyword_set(PPSCIENCE)) then begin
      allplates = [allplates, PPSCIENCE.config]
      allcarts = [allcarts, PPSCIENCE.cartid]
      allModes = [allModes, PPSCIENCE.allModes]
      if keyword_set(fps) then allfields = [allfields, PPSCIENCE.field]
      thismjd = PPSCIENCE[0].mjd
   endif
   if (keyword_set(PPARC)) then begin
      allplates = [allplates, PPARC.config]
      allcarts = [allcarts, PPARC.cartid]
      allModes = [allModes, PPARC.allModes]
      if keyword_set(fps) then allfields = [allfields, PPARC.field]
      thismjd = PPARC[0].mjd
   endif

   allplates = allplates[1:n_elements(allplates)-1]
   allcarts = allcarts[1:n_elements(allcarts)-1]
   allModes = allModes[1:n_elements(allModes)-1]
   if keyword_set(fps) then allfields = allfields[1:n_elements(allfields)-1]
   indx = uniq(allplates, sort(allplates))
   allplates = allplates[indx]
   allcarts = allcarts[indx]
   allModes = allModes[indx]
   if keyword_set(fps) then allfields = allfields[indx]
   nplates = n_elements(allplates)
   mjdstr = strtrim(thismjd,2)

   ;----------
   ; Consruct the header of the output text

   title1 = strupcase(obs)+' BOSS Spectro MJD=' + mjdstr + ' '+var_str+'='
   platelist = var_str+'='
   for iplate=0, nplates-1 do begin
      platestr = strtrim(string(allplates[iplate]),2)
      title1 = title1 + platestr
      platelist = platelist + '<A HREF="#'+var_str1 + platestr + '">' + platestr + '</A>'
      if (iplate NE nplates-1) then begin
         title1 = title1 + ','
         platelist = platelist + ', '
      endif
   endfor
   textout = sos_log_header(title1)
   
   flags = []
   if keyword_set(fps) then flags = [flags, " /fps"]
   if keyword_set(sdssv_sn2) then flags = [flags, " /sdssv_sn2"]
   if keyword_set(sn2_15) then flags = [flags, " /sn2_15"]
   if keyword_set(brightsn2) then flags = [flags, " /brightsn2"]
   if keyword_set(flags) then begin
    flags = strjoin(flags,',')
    textout = [textout, "<!-- Flags:  "+flags+" -->"]
   endif
   
;   textout = [textout, '<FONT SIZE="+4">']
   prevmjd = string(thismjd-1,format='(i5.5)')
   nextmjd = string(thismjd+1,format='(i5.5)')
   prevfile = 'logfile-' + prevmjd + '.html'
   nextfile = 'logfile-' + nextmjd + '.html'
   textout = [textout, $
    '<TABLE CELLSPACING=0 CELLPADDING=0><TR>']
   textout = [textout, $
    '<TD WIDTH="33%" ALIGN="LEFT">Yesterday: ' $
    + '<A HREF=../'+prevmjd+'/'+prevfile+'>MJD='+prevmjd+'</A></TD>']
   textout = [textout, $
    '<TD WIDTH="34%" ALIGN="CENTER"><B><FONT SIZE="+4">BOSS Spectro MJD '+mjdstr+'</FONT></B></TD>']
   textout = [textout, $
    '<TD WIDTH="33%" ALIGN="RIGHT">Tomorrow: ' $
    + '<A HREF=../'+nextmjd+'/'+nextfile+'>MJD='+nextmjd+'</A></TD></TR>']
   textout = [textout, $
    '<TR><TD></TD><TD WIDTH="100%" ALIGN="CENTER"><B><FONT SIZE="+2"><A HREF=../'$
    +mjdstr+'/Summary_'+mjdstr+'.html>SOS Summary Plots</A></FONT></B></TD><TD></TD></TR>']
    
   outf = fileandpath(htmlfile, path=outp)
   arc_html = '../'+mjdstr+'/trace/'+mjdstr+'/arcs_'+mjdstr+'_'+obs+'.html'
    if file_test(FILE_DIRNAME(djs_filepath(arc_html, root_dir=outp)), /DIRECTORY) then begin
        textout = [textout, $
            '<TR><TD></TD><TD WIDTH="100%" ALIGN="CENTER"><B><FONT SIZE="+2">'$
                +'<A HREF='+arc_html+'>Arc Shift Plots</A>'$
                +'</FONT></B></TD><TD></TD></TR>']
   endif
   textout = [textout, $
    '<TR><TD></TD><TD WIDTH="100%" ALIGN="CENTER"><FONT SIZE="+2">'+platelist+'</FONT></TD><TD></TD></TR>']
   textout = [textout, $
    '</TABLE>']

   textout = [textout, $
    '<P>IDLSPEC2D version ' + vers2d ];+ ' ('
    $+ '<A HREF="https://sdss-idlspec2d.readthedocs.io/en/latest/sos.html">documentation</A>).']
   if keyword_set(run2d) then textout = [textout, '<BR> RUN2D '+run2d ]
   if (!version.release LT '5.4') then $
    textout = [textout, $
     '<BR>This page last updated <B>'+systime()+' local time</B>.<P>'] $
   else $
    textout = [textout, $
     '<BR>This page last updated <B>'+systime(/ut)+' UT</B>.<P>']

   ;---------------------------------------------------------------------------
   ; Loop over each plate
   ;---------------------------------------------------------------------------
   disk_warnings = ''
   for iplate=0, nplates-1 do begin

      ;----------
      ; Find which structures correspond to this plate

      thisplate = allplates[iplate]
      thiscart = allcarts[iplate]
      designMode = allModes[iplate]
      if keyword_set(fps) then thisfield = allfields[iplate] else thisfield = thisplate
      textout = [textout, $
       sos_log_beginplate(thisplate, thiscart, thismjd, thisfield, camnames, designMode, outdir=outdir, fps=fps)]

      ;----------
      ; Find all biases and loop over each exposure number with any

      if (keyword_set(PPBIAS)) then ii = where(PPBIAS.config EQ thisplate) $
       else ii = -1
      if (ii[0] NE -1) then begin
         allexp = PPBIAS[ii].expnum
         allexp = allexp[ uniq(allexp, sort(allexp)) ]
         nexp = n_elements(allexp)
         onebias = create_struct(PPBIAS[0], 'PERCENTILE98', 0.0)
         struct_assign, {junk:0}, onebias ; Zero-out all elements
         for iexp=0, nexp-1 do begin
            pbias = replicate(onebias, ncams)
            for icam=0, ncams-1 do begin
               jj = (where(PPBIAS.config EQ thisplate $
                AND PPBIAS.camera EQ camnames[icam] $
                AND PPBIAS.expnum EQ allexp[iexp]))[0]
               if (jj NE -1) then begin
                  copy_struct_inx, PPBIAS[jj], pbias, index_to=icam
                  pbias[icam].percentile98 = PPBIAS[jj].percentile[97]
               endif
            endfor

            mjdstr = strtrim(string(thismjd),2)
            expstring = string(allexp[iexp], format='(i8.8)')
            figl = STRLOWCASE(pbias[0].flavor)+'Plot-'+expstring[0]+'.jpeg'
            figl = ['<A HREF="../'+mjdstr+'/'+ figl + '">PERCENTILE98</A>']

            ; Output table line for this one bias exposure
            fields = ['PERCENTILE98']
            formats = ['(i4)', '(f7.1)', '(f7.1)']
            textout = [ textout, $
             sos_log_fields(pbias, fields, printnames=figl, formats=formats) ]

         endfor
      endif

      ;----------
      ; Find all flats and loop over each exposure number with any

      if (keyword_set(PPFLAT)) then ii = where(PPFLAT.config EQ thisplate) $
       else ii = -1
      if (ii[0] NE -1) then begin
         allexp = PPFLAT[ii].expnum
         allexp = allexp[ uniq(allexp, sort(allexp)) ]
         nexp = n_elements(allexp)
         oneflat = PPFLAT[0]
         struct_assign, {junk:0}, oneflat ; Zero-out all elements
         for iexp=0, nexp-1 do begin
            pflats = replicate(oneflat, ncams)
            for icam=0, ncams-1 do begin
               jj = (where(PPFLAT.config EQ thisplate $
                AND PPFLAT.camera EQ camnames[icam] $
                AND PPFLAT.expnum EQ allexp[iexp]))[0]
               if (jj NE -1) then $
                copy_struct_inx, PPFLAT[jj], pflats, index_to=icam
            endfor

            ; Output table line for this one flat exposure
            fields = ['NGOODFIBER', 'XMID', 'XSIGMA']
            formats = ['(i4)', '(f7.1)', '(f5.2)']
            textout = [ textout, $
             sos_log_fields(pflats, fields, formats=formats) ]

         endfor
      endif

      ;----------
      ; Find all arcs and loop over each exposure number with any

      if (keyword_set(PPARC)) then ii = where(PPARC.config EQ thisplate) $
       else ii = -1
      if (ii[0] NE -1) then begin
         allexp = PPARC[ii].expnum
         allexp = allexp[ uniq(allexp, sort(allexp)) ]
         nexp = n_elements(allexp)
         onearc = PPARC[0]
         struct_assign, {junk:0}, onearc ; Zero-out all elements
         for iexp=0, nexp-1 do begin
            parcs = replicate(onearc, ncams)
            for icam=0, ncams-1 do begin
               jj = (where(PPARC.config EQ thisplate $
                AND PPARC.camera EQ camnames[icam] $
                AND PPARC.expnum EQ allexp[iexp]))[0]
               if (jj NE -1) then $
                copy_struct_inx, PPARC[jj], parcs, index_to=icam
            endfor

            formats = ['(f7.1)', '(f4.2)', '(i)', '(f5.2)']
            fields = ['WAVEMID', 'BESTCORR', 'NLAMPS', 'WSIGMA']
            textout = [ textout, $
             sos_log_fields(parcs, fields, formats=formats) ]
         endfor
      endif

      ;----------
      ; Find all science exposures and collect them into one structure

      ; Now find all unique science exposure numbers for this plate
      if (keyword_set(PPSCIENCE)) then ii = where(PPSCIENCE.config EQ thisplate) $
       else ii = -1
      if (ii[0] NE -1) then begin
         allexp = PPSCIENCE[ii].expnum
         allexp = allexp[ uniq(allexp, sort(allexp)) ]
         nexp = n_elements(allexp)
         onescience = PPSCIENCE[0]
         struct_assign, {junk:0}, onescience ; Zero-out all elements
         pscience = replicate(onescience, ncams, nexp)
         for iexp=0, nexp-1 do begin
            for icam=0, ncams-1 do begin
               jj = (where(PPSCIENCE.config EQ thisplate $
                AND PPSCIENCE.camera EQ camnames[icam] $
                AND PPSCIENCE.expnum EQ allexp[iexp]))[0]
               if (jj NE -1) then $
                copy_struct_inx, PPSCIENCE[jj], pscience, $
                 index_to=icam+iexp*ncams
            endfor
         endfor

         ;----------
         ; Output SKYPERSEC for science exposures

         for iexp=0, nexp-1 do begin
            textout = [ textout, $
             sos_log_fields(pscience[*,iexp], 'SKYPERSEC', $
              printnames='SKY/SEC', formats='(f8.2)') ]
         endfor

         ;----------
         ; Output SN2 for science exposures

         if keyword_set(brightsn2) then begin
            if keyword_set(fps) then this_sn2_15 = sn2_15 else this_sn2_15 = 0
         endif else begin
            if keyword_set(fps) and (long(thisfield) lt 100000) then this_sn2_15 = sn2_15 else this_sn2_15 = 0
         endelse
         this_sdssv_sn2 = 0 ; sdssv_sn2

         for iexp=0, nexp-1 do begin
            mjdstr = strtrim(string(thismjd),2)
            platestr4 = config_to_string(thisplate)
            expstring = string(pscience[*,iexp].expnum, format='(i8.8)')
            
            jpegfile1 = 'snplot-'+mjdstr+'-'+platestr4+'-'+expstring[0]+'.jpeg'
            jpegfile_15= 'snplot-sdssv15-'+mjdstr+'-'+platestr4+'-'+expstring[0]+'.jpeg'
            jpegfile_v2= 'snplot-sdssv-'+mjdstr+'-'+platestr4+'-'+expstring[0]+'.jpeg'
            sn2_lab = ['SN2']
            sn2_pnames = ['<A HREF="../'+mjdstr+'/'+ jpegfile1 + '">(S/N)^2</A>']
            sn2_formats = ['(f7.1)']
            
            if keyword_set(this_sdssv_sn2) then begin
                sn2_lab = [sn2_lab,'SN2_V2']
                sn2_pnames = [sn2_pnames, '<A HREF="../'+mjdstr+'/'+  jpegfile_v2 + '">v2 (S/N)^2</A>']
                sn2_formats = [sn2_formats,'(f7.1)']
            endif

            if keyword_set(this_sn2_15) then begin
                sn2_lab = [sn2_lab,'SN2_15']
                sn2_pnames = [sn2_pnames, '<A HREF="../'+mjdstr+'/'+  jpegfile_15 + '">Mag15 (S/N)^2</A>']
                sn2_formats = [sn2_formats,'(f7.1)']
            endif

            textout = [ textout, $
                        sos_log_fields(pscience[*,iexp], sn2_lab, $
                                       printnames=sn2_pnames,$
                                       formats=sn2_formats) ]

         endfor

         ;----------
         ; Output TOTAL-SN2

         rstruct = create_struct('MJD', 0L, $
                                 'CONFIG', 0L, $
                                 'EXPNUM', '', $
                                 'TAI', '', $
                                 'FLAVOR', 'TOTAL', $
                                 'CAMERA', '', $
                                 'TOTALSN2', 0.0 )
        if keyword_set(this_sdssv_sn2) then begin
            rstruct = struct_addtags(rstruct, create_struct('TOTALSN2_v2', 0.0 ))
        endif
        
        if keyword_set(this_sn2_15) then begin
            rstruct = struct_addtags(rstruct, create_struct('TOTALSN2_15', 0.0 ))
        endif
        
        
         ptotal = replicate(rstruct, ncams)
         for icam=0, ncams-1 do begin
            for iexp=0, nexp-1 do begin
               ; Only add if a 'science' exposure, not a 'smear',
               ; and (S/N)^2 is not flagged as anything bad
               ; in the opLimits file (currently anything < 2.0 is bad).
               if (pscience[icam,iexp].flavor EQ 'science' $
                 AND strmatch(pscience[icam,iexp].quality, 'excellent') $
                 AND sos_checklimits('science', 'SN2', $
                      pscience[icam,iexp].camera, $
                      pscience[icam,iexp].sn2) EQ '') then begin
;                 AND pscience[icam,iexp].sn2 GE 2.0) then begin
                  ptotal[icam].totalsn2 = ptotal[icam].totalsn2 + $
                   pscience[icam,iexp].sn2
               endif
               if keyword_set(this_sdssv_sn2) then begin
                   if (pscience[icam,iexp].flavor EQ 'science' $
                     AND strmatch(pscience[icam,iexp].quality, 'excellent') $
                     AND sos_checklimits('science', 'SN2', $
                          pscience[icam,iexp].camera, $
                          pscience[icam,iexp].sn2_v2) EQ '') then begin
                      ptotal[icam].TOTALSN2_v2 = ptotal[icam].TOTALSN2_v2 + $
                       pscience[icam,iexp].sn2_v2
                   endif
               endif
               if keyword_set(this_sn2_15) then begin
                   if (pscience[icam,iexp].flavor EQ 'science' $
                     AND strmatch(pscience[icam,iexp].quality, 'excellent') $
                     AND sos_checklimits('science', 'SN2', $
                          pscience[icam,iexp].camera, $
                          pscience[icam,iexp].SN2_15) EQ '') then begin
                      ptotal[icam].TOTALSN2_15 = ptotal[icam].TOTALSN2_15 + $
                       pscience[icam,iexp].SN2_15
                   endif
               endif
            endfor
         endfor
         mjdstr = strtrim(string(thismjd),2)
         platestr4 = config_to_string(thisplate)
        jpegfile1 = 'snplot-'+mjdstr+'-'+platestr4+'.jpeg'
        jpegfile_15= 'snplot-sdssv15-'+mjdstr+'-'+platestr4+'.jpeg'
        jpegfile_v2= 'snplot-sdssv-'+mjdstr+'-'+platestr4+'.jpeg'

        sn2_lab = ['TOTALSN2']
        sn2_pnames = ['<A HREF="../'+mjdstr+'/'+ jpegfile1 + '">TOTAL (S/N)^2</A>']
        sn2_formats = ['(f7.1)']
            
        if keyword_set(this_sdssv_sn2) then begin
            sn2_lab = [sn2_lab,'TOTALSN2_V2']
            sn2_pnames = [sn2_pnames, '<A HREF="../'+mjdstr+'/'+  jpegfile_v2 + '">Total v2 (S/N)^2</A>']
            sn2_formats = [sn2_formats,'(f7.1)']
        endif

        if keyword_set(this_sn2_15) then begin
            sn2_lab = [sn2_lab,'TOTALSN2_15']
            sn2_pnames = [sn2_pnames, '<A HREF="../'+mjdstr+'/'+  jpegfile_15 + '">Total Mag15 (S/N)^2</A>']
            sn2_formats = [sn2_formats,'(f7.1)']
        endif
         
        textout = [ textout, $
                    sos_log_fields(ptotal, sn2_lab, $
                                   printnames=sn2_pnames,$
                                   formats=sn2_formats) ]
      endif

      textout = [textout, sos_log_endplate()]

      ;----------
      ; Print all WARNINGs and ABORTs for this plate

      if (keyword_set(PPTEXT)) then ii = where(PPTEXT.config EQ thisplate) $
       else ii = -1
      if (ii[0] NE -1) then begin
         ; Remove leading+trailing spaces
         addtext_temp = strtrim(PPTEXT[ii].text, 2)
         addtext = ''
         ; Remove the first word from each line (which is the name of the
         ; IDL proc that generated the warning or abort message)
         for jj=0, n_elements(addtext_temp)-1 do $
            addtext_temp[jj] = strmid( addtext_temp[jj], strpos(addtext_temp[jj],' ')+1 )
          
         foreach at, addtext_temp do begin
            if strmatch(at, '*SOS disk*', /fold_case) eq 0 then begin
                ;if keyword_set(addtext) then
                addtext = [addtext, at]
            endif else begin
                at_sub = strjoin((strsplit(at,/extract))[-6:*],' ')
                junk = where(strmatch(disk_warnings, at_sub, /fold_case), ct)
                if ct eq 0 then begin
                    addtext = [addtext, at]
                    disk_warnings = [disk_warnings, at_sub]
                endif
            endelse
         endforeach
         
         addtext = repstr(addtext, 'WARNING', $
          '<B><FONT COLOR="' + sos_color2hex('YELLOW') + '">WARNING</FONT></B>')
         addtext = repstr(addtext, 'ABORT', $
          '<B><FONT COLOR="' + sos_color2hex('RED') + '">ABORT</FONT></B>')
         textout = [textout, '<PRE>', addtext, '</PRE>']
      endif
   endfor

   textout = [textout, sos_log_endfile()]

   for i=0, n_elements(textout)-1 do $
    printf, html_lun, textout[i]

   ; Now unlock the HTML file.
   djs_unlockfile, htmlfile, lun=html_lun

   return
end
;------------------------------------------------------------------------------
