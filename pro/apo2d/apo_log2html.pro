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
;   apo_log_science()
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

   textout = '<HTML>'
   textout = [textout, '<HEAD><TITLE>' + title + '</TITLE></HEAD>']
   textout = [textout, '<H2 ALIGN=CENTER>' + title + '</H2>']

   return, textout
end

;------------------------------------------------------------------------------
function apo_log_endfile, title

   textout = '</HTML>'

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
function apo_log_beginplate, platenum, camnames

   rowsep = ' <TR> <TH> '
   colsep = ' <TH> '

   ncams = n_elements(camnames)

   textout = ['<TABLE BORDER=1 CELLPADDING=3>']
   textout = [textout, apo_log_tableline(ncams)]
   textout = [textout, $
    '<CAPTION><H3> PLATE ' + strtrim(string(platenum),2) + '</H3></CAPTION>' ]
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
function apo_log_fields, pp, fields, formats=formats

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
      nextline = colsep + fields[ifield]
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
   while(djs_lockfile(htmlfile, lun=html_lun) EQ 0) do wait, 1

   ; Read in all the HDU's in the log file as structures

   ihdu = 1
   pp = 1
   while (keyword_set(pp)) do begin
      pp = mrdfits(logfile, ihdu)
      if (keyword_set(pp)) then begin
         if (NOT keyword_set(pstruct)) then pstruct = ptr_new(pp) $
          else pstruct = [pstruct, ptr_new(pp)]
      endif
      ihdu = ihdu + 1
   endwhile
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

   ;----------
   ; Consruct the header of the output text

   title = 'APO SPECTRO LOGSHEET MJD=' + strtrim(string(mjd[0]),2) + ' PLATE='
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

      textout = [textout, apo_log_beginplate(thisplate, camnames)]

      ;----------
      ; Find the first (and presumably only) FLAVOR=flat for each camera

      pflats = replicate(ptr_new(), ncams)
      for icam=0, ncams-1 do begin
         ii = where(plate EQ thisplate AND flavor EQ 'flat' $
                    AND camera EQ camnames[icam])
         if (ii[0] NE -1) then pflats[icam] = pstruct[ii[0]]
      endfor
      fields = ['NGOODFIBER', 'XMIN', 'XMAX']
      formats = ['(i4)', '(f7.1)', '(f7.1)']
      textout = [ textout, apo_log_fields(pflats, fields, formats=formats) ]

      ;----------
      ; Find the first (and presumably only) FLAVOR=arc for each camera

      parcs = replicate(ptr_new(), ncams)
      for icam=0, ncams-1 do begin
         ii = where(plate EQ thisplate AND flavor EQ 'arc' $
                    AND camera EQ camnames[icam])
         if (ii[0] NE -1) then parcs[icam] = pstruct[ii[0]]
      endfor
      formats = ['(f7.1)', '(f7.1)', '(f4.2)', '(i)']
      fields = ['WAVEMIN', 'WAVEMAX', 'BESTCORR', 'NLAMPS']
      textout = [ textout, apo_log_fields(parcs, fields, formats=formats) ]

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
               if (jj[0] NE -1) then pscience[icam,iexp] = pstruct[jj[0]]
            endfor
         endfor

         ;----------
         ; Output SKYPERSEC for science exposures

         for iexp=0, nexp-1 do begin
            textout = [ textout, $
             apo_log_fields(pscience[*,iexp], 'SKYPERSEC', formats='(f8.2)') ]
         endfor

         ;----------
         ; Output SN2 for science exposures

         for iexp=0, nexp-1 do begin
            textout = [ textout, $
             apo_log_fields(pscience[*,iexp], 'SN2', formats='(i6)') ]
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
         (*ptotal[ncams-1]).totalsn2 = 1.0e-10 ; Set to non-zero so that this
                                               ; field will be printed.
         for icam=0, ncams-2 do $
          (*ptotal[ncams-1]).totalsn2 = min([ (*ptotal[ncams-1]).totalsn2, $
           (*ptotal[icam]).totalsn2 ]) > 1.0e-10
         textout = [ textout, $
          apo_log_fields(ptotal, 'TOTALSN2', formats='(i6)') ]

      endif

      textout = [textout, apo_log_endplate()]
   endfor

   textout = [textout, apo_log_endfile()]

   for i=0, n_elements(textout)-1 do $
    printf, html_lun, textout[i]

   ; Now unlock the HTML file.
   djs_unlockfile, htmlfile, lun=html_lun

   ; Free pointers
   for i=0, nstruct-1 do $
    if (keyword_set(pstruct[i])) then ptr_free, pstruct[i]

   return
end
;------------------------------------------------------------------------------
