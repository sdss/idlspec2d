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

   jj = where(slimits.field EQ field AND slimits.camera EQ camera)
   if (jj[0] NE -1) then begin
      if (value LT slimits[jj].lovalue $
       OR value GT slimits[jj].hivalue) then $
       markstring = '<B><FONT COLOR="#FF0000">'
   endif

   return, markstring
end

;------------------------------------------------------------------------------
function apo_log_header, title

   ; Include a Java script to auto-load this page every 30 seconds

   textout = '<HTML>'
   textout = [textout, $
'<SCRIPT LANGUAGE="JavaScript"> function reload(){ history.go(0) } </SCRIPT>']
   textout = [textout, '<HEAD><TITLE>' + title + '</TITLE></HEAD>']
   textout = [textout, '<H2 ALIGN=CENTER>' + title + '</H2>']
   squote = "'"
   textout = [textout, $
    '<BODY ONLOAD="timerID=setTimeout('+squote+'reload()'+squote+',30000)">']

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
   platestr = string(platenum, format='(i4.4)')
   plotfile = 'snplot-'+mjdstr+'-'+platestr+'.ps'

   textout = ['<TABLE BORDER=1 CELLPADDING=3>']
   textout = [textout, apo_log_tableline(ncams)]
   textout = [textout, $
    '<CAPTION><B> PLATE ' + strtrim(string(platenum),2) + '</B>' $
    + ' - <A HREF="' + plotfile + '">S/N FIGURE</A>' + '</CAPTION>' ]
   nextline = rowsep + colsep
   for icam=0, ncams-1 do $
    nextline = nextline + colsep + camnames[icam]
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
   igood = where(pp NE ptr_new())
   if (igood[0] EQ -1) then return, ''

   flavor = (*pp[igood[0]]).flavor
   expstring = strtrim(string( (*pp[igood[0]]).expnum ),2)
   tags = tag_names(*pp[igood[0]])

   for ifield=0, n_elements(fields)-1 do begin
      itag = (where(fields[ifield] EQ tags))[0]
      if (keyword_set(printnames)) then nextline = colsep + printnames[ifield] $
       else nextline = colsep + fields[ifield]
      if (keyword_set(formats)) then format = formats[ifield]
      for icam=0, ncams-1 do begin
         value = ' '
         if (keyword_set(pp[icam])) then begin
            tmpval = (*pp[icam]).(itag)
            if (keyword_set(tmpval)) then $
             value = string(tmpval, format=format)
            value = apo_checklimits(fields[ifield], camnames[icam], tmpval) $
             + value
         endif
         nextline = nextline + colsep + value
      endfor
      nextline = nextline + colsep

      if (ifield EQ 0) then $
       textout = rowsep + strupcase(flavor) + '-' + expstring + nextline $
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

   camnames = ['b1', 'r1', 'b2', 'r2', 'ALL']
   ncams = n_elements(camnames)

   ; Lock the file to do this.
   while(djs_lockfile(htmlfile, lun=html_lun) EQ 0) do wait, 5

   ; Read in all the HDU's in the log file as structures
   pstruct = apo_readlog(logfile)
   if (NOT keyword_set(pstruct)) then return
   nstruct = n_elements(pstruct)

   mjd = lonarr(nstruct)
   plate = lonarr(nstruct)
   expnum = lonarr(nstruct)
   flavor = strarr(nstruct)
   camera = strarr(nstruct)
   for ii=0, nstruct-1 do begin
      mjd[ii] = (*pstruct[ii]).mjd
      plate[ii] = (*pstruct[ii]).plate
      expnum[ii] = (*pstruct[ii]).expnum
      flavor[ii] = (*pstruct[ii]).flavor
      camera[ii] = (*pstruct[ii]).camera
   endfor

   allplates = plate[ uniq(plate, sort(plate)) ]
   nplates = n_elements(allplates)
   mjdstr = strtrim(string(mjd[0]),2)

   ;----------
   ; Consruct the header of the output text

   title = 'APO SPECTRO LOGSHEET MJD=' + mjdstr + ' PLATE='
   for iplate=0, nplates-1 do begin
      title = title + strtrim(string(allplates[iplate]),2)
      if (iplate NE nplates-1) then title=title+','
   endfor
   textout = apo_log_header(title)

   ;---------------------------------------------------------------------------
   ; Loop over each plate
   ;---------------------------------------------------------------------------

   for iplate=0, nplates-1 do begin

      ;----------
      ; Find which structures correspond to this plate

      thisplate = allplates[iplate]

      textout = [textout, apo_log_beginplate(thisplate, mjd[0], camnames)]

      ;----------
      ; Find all flats and loop over each exposure number with any

      ii = where(plate EQ thisplate AND flavor EQ 'flat')
      if (ii[0] NE -1) then begin
         allexp = expnum[ii[ uniq(expnum[ii], sort(expnum[ii])) ]]
         nexp = n_elements(allexp)
         pflats = replicate(ptr_new(), ncams)
         for iexp=0, nexp-1 do begin
            for icam=0, ncams-1 do begin
               jj = where(plate EQ thisplate AND flavor EQ 'flat' $
                AND camera EQ camnames[icam] AND expnum EQ allexp[iexp])
               if (jj[0] NE -1) then pflats[icam] = pstruct[jj[0]] $
                else pflats[icam] = ptr_new()
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

      ii = where(plate EQ thisplate AND flavor EQ 'arc')
      if (ii[0] NE -1) then begin
         allexp = expnum[ii[ uniq(expnum[ii], sort(expnum[ii])) ]]
         nexp = n_elements(allexp)
         parcs = replicate(ptr_new(), ncams)
         for iexp=0, nexp-1 do begin
            for icam=0, ncams-1 do begin
               jj = where(plate EQ thisplate AND flavor EQ 'arc' $
                AND camera EQ camnames[icam] AND expnum EQ allexp[iexp])
               if (jj[0] NE -1) then parcs[icam] = pstruct[jj[0]] $
                else parcs[icam] = ptr_new()
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
      ii = where(plate EQ thisplate AND flavor EQ 'science')
      if (ii[0] NE -1) then begin
         allexp = expnum[ii[ uniq(expnum[ii], sort(expnum[ii])) ]]
         nexp = n_elements(allexp)

         pscience = replicate(ptr_new(), ncams, nexp)
         for iexp=0, nexp-1 do begin
            for icam=0, ncams-1 do begin
               jj = where(plate EQ thisplate AND flavor EQ 'science' $
                AND camera EQ camnames[icam] AND expnum EQ allexp[iexp])
               if (jj[0] NE -1) then pscience[icam,iexp] = pstruct[jj[0]] $
                else pscience[icam,iexp] = ptr_new()
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
              printnames='(S/N)^2', formats='(f6.0)') ]
         endfor

         ;----------
         ; Output TOTAL-SN2

         rstruct = create_struct('MJD', 0L, $
                                 'PLATE', 0L, $
                                 'EXPNUM', '', $
                                 'FLAVOR', '', $
                                 'CAMERA', '', $
                                 'TOTALSN2', 0.0 )
         ptotal = replicate(ptr_new(), ncams)
         for icam=0, ncams-2 do begin
            ptotal[icam] = ptr_new(rstruct)
            for iexp=0, nexp-1 do begin
               if (keyword_set(pscience[icam,iexp])) then begin
                  (*ptotal[icam]).totalsn2 = (*ptotal[icam]).totalsn2 + $
                   (*pscience[icam,iexp]).sn2
               endif
            endfor
         endfor
         ptotal[ncams-1] = ptr_new(rstruct) ; 'ALL' camera
         (*ptotal[ncams-1]).totalsn2 = 9999
         for icam=0, ncams-2 do $
          (*ptotal[ncams-1]).totalsn2 = min([ (*ptotal[ncams-1]).totalsn2, $
           (*ptotal[icam]).totalsn2 ]) > 1.0e-10 ; Set to non-zero so that this
                                                 ; field will be printed.
         textout = [ textout, $
          apo_log_fields(ptotal, 'TOTALSN2', $
           printnames='TOTAL (S/N)^2', formats='(f6.0)') ]
      endif

      textout = [textout, apo_log_endplate()]

      ;----------
      ; Append all WARNINGs and ABORTs for this plate

      ii = where(plate EQ thisplate)
      warnings = ''
      aborts = ''
      for j=0, n_elements(ii)-1 do begin
         warnings = [warnings, strtrim((*pstruct[ii[j]]).warnings,2)]
         aborts = [aborts, strtrim((*pstruct[ii[j]]).aborts,2)]
      endfor
      j = where(warnings NE '')
      if (j[0] NE -1) then warnings = warnings[j] $
       else warnings = ''
      j = where(aborts NE '')
      if (j[0] NE -1) then aborts = aborts[j] $
       else aborts = ''

      textout = [textout, '<PRE>', warnings, '</PRE>']
      textout = [textout, '<PRE>', aborts, '</PRE>']
   endfor

   textout = [textout, apo_log_endfile()]

   for i=0, n_elements(textout)-1 do $
    printf, html_lun, textout[i]

   ; Now unlock the HTML file.
   djs_unlockfile, htmlfile, lun=html_lun

   ; Free pointers
   for i=0, nstruct-1 do $
    if (keyword_set(pstruct[i])) then ptr_free, pstruct[i]
   heap_gc ; Just in case we missed some!

   return
end
;------------------------------------------------------------------------------
