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
;   djs_lockfile()
;   djs_unlockfile
;
; INTERNAL SUPPORT ROUTINES:
;   apo_color2hex()
;   apo_stringreplace()
;   apo_log_header()
;   apo_log_endfile()
;   apo_log_tableline()
;   apo_log_beginplate()
;   apo_log_endplate()
;   apo_log_fields()
;
; REVISION HISTORY:
;
; REVISION HISTORY:
;   30-Apr-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
function apo_color2hex, colorname

   case strupcase(strtrim(colorname,2)) of
   'RED': hexname = '#FF0000'
   'YELLOW': hexname = '#9F9F00'
   endcase

   return, hexname
end

;------------------------------------------------------------------------------
; The following replaces only the first instance of STRING1 within TEXT
; with STRING2.
function apo_stringreplace, text, string1, string2

   newtext = ''

   i = strpos(text, string1)
   j = strlen(string1)
   k = strlen(text)
   if (i[0] NE -1) then begin
      if (i GT 0) then newtext = newtext + strmid(text,0,i)
      newtext = newtext + string2
      if (i+j LT k) then newtext = newtext + strmid(text,i+j,k-i-j)
   endif

   return, newtext
end
;------------------------------------------------------------------------------
function apo_checklimits, field, camera, value

   common apo_limits, slimits

   markstring = ''
   if (NOT keyword_set(value)) then return, markstring

   ; Read this Yanny file only the first time this routine is called,
   ; then save the limits in a common block.
   if (NOT keyword_set(slimits)) then begin
      limitfile = filepath('opLimits.par', root_dir=getenv('IDLSPEC2D_DIR'), $
       subdirectory='examples')
      yanny_read, limitfile, pdata
      slimits = *pdata[0]
      yanny_free, pdata
   endif

   indx = where(slimits.field EQ field AND slimits.camera EQ camera, nlim)
   for ilim=0, nlim-1 do begin
      if (value GE slimits[indx[ilim]].lovalue $
       AND value LE slimits[indx[ilim]].hivalue) then $
       markstring = '<B><FONT COLOR="' $
        + apo_color2hex(slimits[indx[ilim]].color) + '">'
   endfor

   return, markstring
end

;------------------------------------------------------------------------------
function apo_log_header, title1, title2

   ; Include a Java script to auto-load this page every 60 seconds

   textout = '<HTML>'
   textout = [textout, '<HEAD><TITLE>' + title1 + '</TITLE></HEAD>']
   textout = [textout, '<H1 ALIGN=CENTER>' + title2 + '</H1>']
   squote = "'"
   textout = [textout, $
    '<BODY ONLOAD="timerID=setTimeout('+squote+'location.reload(true)'+squote+',60000)">']

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
function apo_log_beginplate, platenum, mjd, camnames

   rowsep = ' <TR> <TH> '
   colsep = ' <TH> '

   ncams = n_elements(camnames)

   mjdstr = strtrim(string(mjd),2)
   platestr = strtrim(string(platenum),2)
   platestr4 = string(platenum, format='(i4.4)')
   plotfile = 'snplot-'+mjdstr+'-'+platestr4+'.ps'

   textout = ['<A NAME="PLATE' + platestr + '">']
   textout = [textout, '<TABLE BORDER=1 CELLPADDING=3>']
   textout = [textout, apo_log_tableline(ncams)]
   textout = [textout, $
    '<CAPTION><FONT SIZE="+3"><B> PLATE ' + platestr + '</B></FONT>' $
    + ' - <A HREF="' + plotfile + '">S/N FIGURE</A>' + '</CAPTION>' ]
   nextline = rowsep + colsep
   for icam=0, ncams-1 do $
    nextline = nextline + colsep + camnames[icam]
   nextline = nextline + colsep + 'MST'
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

   flavor = pp[igood[0]].flavor
   expstring = strtrim(string( pp[igood[0]].expnum ),2)
   mststring = strtrim(string( pp[igood[0]].mst ),2)
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
            value = apo_checklimits(fields[ifield], camnames[icam], tmpval) $
             + value
         endif
         nextline = nextline + colsep + value
      endfor
      nextline = nextline + colsep

      if (ifield EQ 0) then $
       textout = rowsep + strupcase(flavor) + '-' + expstring $
        + nextline + mststring + colsep $
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
      ipos = rstrpos(logfile, '.fits')
      if (ipos EQ -1) then ipos = strlen(logfile)
      htmlfile = strmid(logfile, 0, ipos) + '.html'
   endif

   camnames = ['b1', 'r1', 'b2', 'r2']
   ncams = n_elements(camnames)

   ; Lock the file to do this.
   while(djs_lockfile(htmlfile, lun=html_lun) EQ 0) do wait, 5

   ; Read in all the HDU's in the log file as structures
   PPFLAT = mrdfits(logfile, 1)
   PPARC = mrdfits(logfile, 2)
   PPSCIENCE = mrdfits(logfile, 3)
   if (NOT keyword_set(PPFLAT)) then begin
      djs_unlockfile, htmlfile, lun=html_lun
      return
   endif

   allplates = PPFLAT[ uniq(PPFLAT.plate, sort(PPFLAT.plate)) ].plate
   nplates = n_elements(allplates)
   mjdstr = strtrim(string(PPFLAT[0].mjd),2)

   ;----------
   ; Consruct the header of the output text

   title1 = 'APO SPECTRO LOGSHEET MJD=' + mjdstr + ' PLATE='
   title2 = 'APO SPECTRO LOGSHEET MJD=' + mjdstr + '<BR> PLATE='
   for iplate=0, nplates-1 do begin
      platestr = strtrim(string(allplates[iplate]),2)
      title1 = title1 + platestr
      title2 = title2 + '<A HREF="#PLATE' + platestr + '">' + platestr + '</A>'
      if (iplate NE nplates-1) then begin
         title1 = title1 + ','
         title2 = title2 + ','
      endif
   endfor
   textout = apo_log_header(title1, title2)

   ;---------------------------------------------------------------------------
   ; Loop over each plate
   ;---------------------------------------------------------------------------

   for iplate=0, nplates-1 do begin

      ;----------
      ; Find which structures correspond to this plate

      thisplate = allplates[iplate]

      textout = [textout, $
       apo_log_beginplate(thisplate, PPFLAT[0].mjd, camnames)]

      ;----------
      ; Append all WARNINGs and ABORTs for this plate to the following

      warnings = ''
      aborts = ''

      ;----------
      ; Find all flats and loop over each exposure number with any

      if (keyword_set(PPFLAT)) then ii = where(PPFLAT.plate EQ thisplate) $
       else ii = -1
      if (ii[0] NE -1) then begin
         warnings = [warnings, PPFLAT[ii].warnings]
         aborts = [aborts, PPFLAT[ii].aborts]

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
            fields = ['NGOODFIBER', 'XMIN', 'XMAX']
            formats = ['(i4)', '(f7.1)', '(f7.1)']
            textout = [ textout, $
             apo_log_fields(pflats, fields, formats=formats) ]

         endfor
      endif

      ;----------
      ; Find all arcs and loop over each exposure number with any

      if (keyword_set(PPARC)) then ii = where(PPARC.plate EQ thisplate) $
       else ii = -1
      if (ii[0] NE -1) then begin
         warnings = [warnings, PPARC[ii].warnings]
         aborts = [aborts, PPARC[ii].aborts]

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

            formats = ['(f7.1)', '(f7.1)', '(f4.2)', '(i)']
            fields = ['WAVEMIN', 'WAVEMAX', 'BESTCORR', 'NLAMPS']
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
         warnings = [warnings, PPSCIENCE[ii].warnings]
         aborts = [aborts, PPSCIENCE[ii].aborts]

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
            textout = [ textout, $
             apo_log_fields(pscience[*,iexp], 'SN2', $
              printnames='(S/N)^2', formats='(f7.1)') ]
         endfor

         ;----------
         ; Output TOTAL-SN2

         rstruct = create_struct('MJD', 0L, $
                                 'PLATE', 0L, $
                                 'EXPNUM', '', $
                                 'MST', '', $
                                 'FLAVOR', 'TOTAL', $
                                 'CAMERA', '', $
                                 'TOTALSN2', 0.0 )
         ptotal = replicate(rstruct, ncams)
         for icam=0, ncams-1 do begin
            for iexp=0, nexp-1 do begin
               if (keyword_set(pscience[icam,iexp])) then begin
                  ptotal[icam].totalsn2 = ptotal[icam].totalsn2 + $
                   pscience[icam,iexp].sn2
               endif
            endfor
         endfor
         textout = [ textout, $
          apo_log_fields(ptotal, 'TOTALSN2', $
           printnames='TOTAL (S/N)^2', formats='(f7.1)') ]
      endif

      textout = [textout, apo_log_endplate()]

      ;----------
      ; Print all WARNINGs and ABORTs for this plate

      for j=0, n_elements(warnings)-1 do $
       warnings[j] = apo_stringreplace(warnings[j], 'WARNING', $
        '<B><FONT COLOR="' + apo_color2hex('YELLOW') + '">WARNING</FONT></B>')

      for j=0, n_elements(aborts)-1 do $
       aborts[j] = apo_stringreplace(aborts[j], 'ABORT', $
        '<B><FONT COLOR="' + apo_color2hex('RED') + '">ABORT</FONT></B>')

      j = where(warnings NE '')
      if (j[0] NE -1) then $
       textout = [textout, '<PRE>', warnings[j], '</PRE>']

      j = where(aborts NE '')
      if (j[0] NE -1) then $
       textout = [textout, '<PRE>', aborts[j], '</PRE>']
   endfor

   textout = [textout, apo_log_endfile()]

   for i=0, n_elements(textout)-1 do $
    printf, html_lun, textout[i]

   ; Now unlock the HTML file.
   djs_unlockfile, htmlfile, lun=html_lun

   return
end
;------------------------------------------------------------------------------
