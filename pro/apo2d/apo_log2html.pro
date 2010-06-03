;+
; NAME:
;   apo_log2html
;
; PURPOSE:
;   Convert output FITS file from APOREDUCE to HTML format.
;
; CALLING SEQUENCE:
;   apo_log2html, logfile, [ htmlfile ]
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
;   apo_checklimits()
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
;   apo_color2hex()
;   apo_log_header()
;   apo_log_endfile()
;   apo_log_tableline()
;   apo_log_beginplate()
;   apo_log_endplate()
;   apo_log_fields()
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
function apo_log_header, title1

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
function apo_log_endfile, title

   textout = '</BODY></HTML>'

   return, textout
end

;------------------------------------------------------------------------------
function apo_log_tableline, ncams

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
function apo_log_beginplate, platenum, cartid, mjd, camnames, outdir=outdir

   rowsep = ' <TR> <TH> '
   colsep = ' <TH> '

   ncams = n_elements(camnames)

   mjdstr = strtrim(string(mjd),2)
   platestr = strtrim(string(platenum),2)
   cartstr = strtrim(string(cartid),2)
   platestr4 = string(platenum, format='(i4.4)')
   plotfile = 'snplot-'+mjdstr+'-'+platestr4+'.ps'
   jpegfile = 'snplot-'+mjdstr+'-'+platestr4+'.jpeg'

   textout = ['<A NAME="PLATE' + platestr + '">']
   textout = [textout, '<TABLE BORDER=1 CELLPADDING=3>']
   textout = [textout, apo_log_tableline(ncams)]
   nextline = '<CAPTION><FONT SIZE="+2"><B> Plate ' + platestr $
    + ' on Cart #' + cartstr + '</B></FONT>'
   textout = [textout, nextline]
   nextline = rowsep + colsep
   for icam=0, ncams-1 do $
    nextline = nextline + colsep + camnames[icam]
   nextline = nextline + colsep + 'EXPTIME' + colsep + 'TEMP' $
    + colsep + 'UT' + colsep + 'QUALITY'
   textout = [textout, nextline]

   textout = [textout, apo_log_tableline(ncams)]

   return, textout
end

;------------------------------------------------------------------------------
function apo_log_endplate

   textout = ['</TABLE>']

   return, textout
end

;------------------------------------------------------------------------------
function apo_log_fields, pp, fields, printnames=printnames, formats=formats

   common com_apo_log, camnames

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
      exptimestring = apo_checklimits(pp[igood[0]].flavor, 'EXPTIME', $
       pp[igood[0]].camera, pp[igood[0]].exptime, /html) $
       + string(pp[igood[0]].exptime, format='(f8.1)')
      qualstring = apo_checklimits(pp[igood[0]].flavor, 'QUALITY', $
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
            value = apo_checklimits(flavor, fields[ifield], $
             camnames[icam], tmpval, /html) + value
         endif
         nextline = nextline + colsep + value
      endfor
      nextline = nextline + colsep

      if (ifield EQ 0) then $
       textout = rowsep + strupcase(flavor) + '-' + expstring $
        + nextline + exptimestring + colsep $
        + airtempstring + colsep + utstring + colsep + qualstring + colsep  $
      else $
       textout = [textout, rowsep + nextline]
   endfor

   textout = [textout, apo_log_tableline(ncams)]

   return, textout
end

;------------------------------------------------------------------------------
pro apo_log2html, logfile, htmlfile

   common com_apo_log, camnames

   if (n_params() EQ 0) then begin
      print, 'Syntax: apo_log2html, logfile, [ htmlfile ]'
      return
   endif else if (n_params() EQ 1) then begin
      thisfile = fileandpath(logfile, path=thispath)
      ipos = strpos(thisfile, '.fits', /reverse_search)
      if (ipos EQ -1) then ipos = strlen(thisfile)
      htmlfile = djs_filepath(strmid(thisfile, 0, ipos) + '.html', $
       root_dir=thispath)
   endif

   junk = fileandpath(htmlfile, path=outdir)
   camnames = ['b1', 'r1', 'b2', 'r2']
   ncams = n_elements(camnames)

   ; Lock the files.
   while(djs_lockfile(htmlfile, lun=html_lun) EQ 0) do wait, 5
   while(djs_lockfile(logfile) EQ 0) do wait, 5

   ; Read the 0th header to get the version of the code
   hdr = headfits(logfile)
   vers2d = sxpar(hdr, 'VERS2D')

   ; Read in all the HDU's in the log file as structures
   PPBIAS = mrdfits(logfile, 1)
   PPFLAT = mrdfits(logfile, 2)
   PPARC = mrdfits(logfile, 3)
   PPSCIENCE = mrdfits(logfile, 4)
   PPTEXT = mrdfits(logfile, 5)
   djs_unlockfile, logfile
   if (NOT keyword_set(PPBIAS) AND NOT keyword_set(PPFLAT) $
    AND NOT keyword_set(PPTEXT)) then begin
      djs_unlockfile, htmlfile, lun=html_lun
      return
   endif

   allplates = [0]
   allcarts = [0]
   if (keyword_set(PPBIAS)) then begin
      allplates = [allplates, PPBIAS.plate]
      allcarts = [allcarts, PPBIAS.cartid]
      thismjd = PPBIAS[0].mjd
   endif
   if (keyword_set(PPFLAT)) then begin
      allplates = [allplates, PPFLAT.plate]
      allcarts = [allcarts, PPFLAT.cartid]
      thismjd = PPFLAT[0].mjd
   endif
   if (keyword_set(PPTEXT)) then begin
      allplates = [allplates, PPTEXT.plate]
      allcarts = [allcarts, PPTEXT.cartid]
      thismjd = PPTEXT[0].mjd
   endif
   allplates = allplates[1:n_elements(allplates)-1]
   allcarts = allcarts[1:n_elements(allcarts)-1]
   indx = uniq(allplates, sort(allplates))
   allplates = allplates[indx]
   allcarts = allcarts[indx]
   nplates = n_elements(allplates)
   mjdstr = strtrim(thismjd,2)

   ;----------
   ; Consruct the header of the output text

   title1 = 'BOSS Spectro MJD=' + mjdstr + ' Plate='
   platelist = 'Plate='
   for iplate=0, nplates-1 do begin
      platestr = strtrim(string(allplates[iplate]),2)
      title1 = title1 + platestr
      platelist = platelist + '<A HREF="#PLATE' + platestr + '">' + platestr + '</A>'
      if (iplate NE nplates-1) then begin
         title1 = title1 + ','
         platelist = platelist + ','
      endif
   endfor
   textout = apo_log_header(title1)

;   textout = [textout, '<FONT SIZE="+4">']
   prevmjd = string(thismjd-1,format='(i5.5)')
   nextmjd = string(thismjd+1,format='(i5.5)')
   prevfile = 'logfile-' + prevmjd + '.html'
   nextfile = 'logfile-' + nextmjd + '.html'
   textout = [textout, $
    '<TABLE CELLSPACING=0 CELLPADDING=0><TR>']
   textout = [textout, $
    '<TD WIDTH="33%" ALIGN="LEFT">Yesterday: ' $
    + '<A HREF='+prevfile+'>MJD='+prevmjd+'</A></TD>']
   textout = [textout, $
    '<TD WIDTH="34%" ALIGN="CENTER"><B><FONT SIZE="+4">BOSS Spectro MJD '+mjdstr+'</FONT></B></TD>']
   textout = [textout, $
    '<TD WIDTH="33%" ALIGN="RIGHT">Tomorrow: ' $
    + '<A HREF='+nextfile+'>MJD='+nextmjd+'</A></TD></TR>']
   textout = [textout, $
    '<TR><TD></TD><TD WIDTH="100%" ALIGN="CENTER"><FONT SIZE="+2">'+platelist+'</FONT></TD><TD></TD></TR>']
   textout = [textout, $
    '</TABLE>']
;   textout = [textout, '<ALIGN=CENTER>' + platelist + '</ALIGN>']
;   textout = [textout, '<FONT SIZE="+0">']

   textout = [textout, $
    '<P>IDLSPEC2D version ' + vers2d + ' (' $
    + '<A HREF="http://sdsshost2.apo.nmsu.edu/doc/idlspec2d/spectroSOS.html">documentation</A>).']
   if (!version.release LT '5.4') then $
    textout = [textout, $
     '<BR>This page last updated <B>'+systime()+' local time</B>.<P>'] $
   else $
    textout = [textout, $
     '<BR>This page last updated <B>'+systime(/ut)+' UT</B>.<P>']

   ;---------------------------------------------------------------------------
   ; Loop over each plate
   ;---------------------------------------------------------------------------

   for iplate=0, nplates-1 do begin

      ;----------
      ; Find which structures correspond to this plate

      thisplate = allplates[iplate]
      thiscart = allcarts[iplate]

      textout = [textout, $
       apo_log_beginplate(thisplate, thiscart, thismjd, camnames, outdir=outdir)]

      ;----------
      ; Find all biases and loop over each exposure number with any

      if (keyword_set(PPBIAS)) then ii = where(PPBIAS.plate EQ thisplate) $
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
               jj = (where(PPBIAS.plate EQ thisplate $
                AND PPBIAS.camera EQ camnames[icam] $
                AND PPBIAS.expnum EQ allexp[iexp]))[0]
               if (jj NE -1) then begin
                  copy_struct_inx, PPBIAS[jj], pbias, index_to=icam
                  pbias[icam].percentile98 = PPBIAS[jj].percentile[97]
               endif
            endfor

            ; Output table line for this one bias exposure
            fields = ['PERCENTILE98']
            formats = ['(i4)', '(f7.1)', '(f7.1)']
            textout = [ textout, $
             apo_log_fields(pbias, fields, formats=formats) ]

         endfor
      endif

      ;----------
      ; Find all flats and loop over each exposure number with any

      if (keyword_set(PPFLAT)) then ii = where(PPFLAT.plate EQ thisplate) $
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
               jj = (where(PPFLAT.plate EQ thisplate $
                AND PPFLAT.camera EQ camnames[icam] $
                AND PPFLAT.expnum EQ allexp[iexp]))[0]
               if (jj NE -1) then $
                copy_struct_inx, PPFLAT[jj], pflats, index_to=icam
            endfor

            ; Output table line for this one flat exposure
            fields = ['NGOODFIBER', 'XMID', 'XSIGMA']
            formats = ['(i4)', '(f7.1)', '(f5.2)']
            textout = [ textout, $
             apo_log_fields(pflats, fields, formats=formats) ]

         endfor
      endif

      ;----------
      ; Find all arcs and loop over each exposure number with any

      if (keyword_set(PPARC)) then ii = where(PPARC.plate EQ thisplate) $
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
               jj = (where(PPARC.plate EQ thisplate $
                AND PPARC.camera EQ camnames[icam] $
                AND PPARC.expnum EQ allexp[iexp]))[0]
               if (jj NE -1) then $
                copy_struct_inx, PPARC[jj], parcs, index_to=icam
            endfor

            formats = ['(f7.1)', '(f4.2)', '(i)', '(f5.2)']
            fields = ['WAVEMID', 'BESTCORR', 'NLAMPS', 'WSIGMA']
            textout = [ textout, $
             apo_log_fields(parcs, fields, formats=formats) ]
         endfor
      endif

      ;----------
      ; Find all science exposures and collect them into one structure

      ; Now find all unique science exposure numbers for this plate
      if (keyword_set(PPSCIENCE)) then ii = where(PPSCIENCE.plate EQ thisplate) $
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
               jj = (where(PPSCIENCE.plate EQ thisplate $
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
             apo_log_fields(pscience[*,iexp], 'SKYPERSEC', $
              printnames='SKY/SEC', formats='(f8.2)') ]
         endfor

         ;----------
         ; Output SN2 for science exposures

         for iexp=0, nexp-1 do begin
            mjdstr = strtrim(string(thismjd),2)
            platestr4 = string(thisplate, format='(i4.4)')
            expstring = string(pscience[*,iexp].expnum, format='(i8.8)')
            jpegfile1 = 'snplot-'+mjdstr+'-'+platestr4+'-'+expstring+'.jpeg'
            printnames = '<A HREF="' + jpegfile1 + '">(S/N)^2</A>'
            textout = [ textout, $
             apo_log_fields(pscience[*,iexp], 'SN2', $
             printnames=printnames, formats='(f7.1)') ]
         endfor

         ;----------
         ; Output TOTAL-SN2

         rstruct = create_struct('MJD', 0L, $
                                 'PLATE', 0L, $
                                 'EXPNUM', '', $
                                 'TAI', '', $
                                 'FLAVOR', 'TOTAL', $
                                 'CAMERA', '', $
                                 'TOTALSN2', 0.0 )
         ptotal = replicate(rstruct, ncams)
         for icam=0, ncams-1 do begin
            for iexp=0, nexp-1 do begin
               ; Only add if a 'science' exposure, not a 'smear',
               ; and (S/N)^2 is not flagged as anything bad
               ; in the opLimits file (currently anything < 2.0 is bad).
               if (pscience[icam,iexp].flavor EQ 'science' $
                 AND strmatch(pscience[icam,iexp].quality, 'excellent') $
                 AND apo_checklimits('science', 'SN2', $
                      pscience[icam,iexp].camera, $
                      pscience[icam,iexp].sn2) EQ '') then begin
;                 AND pscience[icam,iexp].sn2 GE 2.0) then begin
                  ptotal[icam].totalsn2 = ptotal[icam].totalsn2 + $
                   pscience[icam,iexp].sn2
               endif
            endfor
         endfor
         mjdstr = strtrim(string(thismjd),2)
         platestr4 = string(thisplate, format='(i4.4)')
         jpegfile = 'snplot-'+mjdstr+'-'+platestr4+'.jpeg'
         printnames = '<A HREF="' + jpegfile + '">TOTAL (S/N)^2</A>'
         textout = [ textout, $
          apo_log_fields(ptotal, 'TOTALSN2', $
           printnames=printnames, formats='(f7.1)') ]
      endif

      textout = [textout, apo_log_endplate()]

      ;----------
      ; Print all WARNINGs and ABORTs for this plate

      if (keyword_set(PPTEXT)) then ii = where(PPTEXT.plate EQ thisplate) $
       else ii = -1
      if (ii[0] NE -1) then begin
         ; Remove leading+trailing spaces
         addtext = strtrim(PPTEXT[ii].text, 2)
         ; Remove the first word from each line (which is the name of the
         ; IDL proc that generated the warning or abort message)
         for jj=0, n_elements(addtext)-1 do $
          addtext[jj] = strmid( addtext[jj], strpos(addtext[jj],' ')+1 )
         addtext = repstr(addtext, 'WARNING', $
          '<B><FONT COLOR="' + apo_color2hex('YELLOW') + '">WARNING</FONT></B>')
         addtext = repstr(addtext, 'ABORT', $
          '<B><FONT COLOR="' + apo_color2hex('RED') + '">ABORT</FONT></B>')
         textout = [textout, '<PRE>', addtext, '</PRE>']
      endif
   endfor

   textout = [textout, apo_log_endfile()]

   for i=0, n_elements(textout)-1 do $
    printf, html_lun, textout[i]

   ; Now unlock the HTML file.
   djs_unlockfile, htmlfile, lun=html_lun

   return
end
;------------------------------------------------------------------------------
