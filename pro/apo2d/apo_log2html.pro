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
function apo_log_header, title

   textout = '<HTML>'
   textout = [textout, '<HEAD><TITLE>' + title + '</TITLE></HEAD>']
   textout = [textout, '<H1 ALIGN="center">' + title + '</H1>']

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

   textout = ['<TABLE BORDER=1>']
   textout = [textout, apo_log_tableline(ncams)]
   textout = [textout, $
    '<CAPTION><B> PLATE ' + strtrim(string(platenum),2) + '</B></CAPTION>' ]
   nextline = rowsep + colsep
   for icam=0, ncams-1 do $
    nextline = nextline + colsep + camnames[icam]
   nextline = nextline + colsep + 'ALL'
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
function apo_log_fields, pp, fields

   rowsep = ' <TR> <TH> '
   colsep = ' <TD> '

   ncams = n_elements(pp)
   igood = where(keyword_set(pp))
   if (igood[0] EQ -1) then return, ''

   flavor = (*pp[igood[0]]).flavor
   tags = tag_names(*pp[igood[0]])

   for ifield=0, n_elements(fields)-1 do begin
      itag = (where(fields[ifield] EQ tags))[0]
      nextline = colsep + fields[ifield]
      for icam=0, ncams-1 do begin
         if (keyword_set(pp[icam])) then value = (*pp[icam]).(itag) $
          else value = ' '
         nextline = nextline + colsep + string(value)
      endfor
      nextline = nextline + colsep

      if (ifield EQ 0) then textout = rowsep + strupcase(flavor) + nextline $
       else textout = [textout, rowsep + nextline]
   endfor

   textout = [textout, apo_log_tableline(ncams)]

   return, textout
end

;------------------------------------------------------------------------------
function apo_log_science, pp, fieldname, total=total

   rowsep = ' <TR> <TH> '
   colsep = ' <TD> '

   ndim = size(pp, /n_dimen)
   dims = size(pp, /dimens)
   ncams = dims[0]
   if (ndim EQ 1) then nexp = 1 $
    else nexp = dims[1]

   igood = where(keyword_set(pp))
   if (igood[0] EQ -1) then return, ''

   flavor = (*pp[igood[0]]).flavor
   tags = tag_names(*pp[igood[0]])

   totals = fltarr(ncams)

   for iexp=0, nexp-1 do begin
      ; Get the exposure number
      igood = where(keyword_set(pp[*,iexp]))
      if (igood[0] NE -1) then expnum = (*pp[igood[0],iexp]).expnum $
       else expnum = 0

      itag = (where(fieldname EQ tags))[0]
      nextline = strupcase(flavor) + '-' $
       + strtrim(string(long(expnum)),2) $
       + colsep + fieldname
      for icam=0, ncams-1 do begin
         if (keyword_set(pp[icam,iexp])) then value = (*pp[icam,iexp]).(itag) $
          else value = ''
         nextline = nextline + colsep + string(value)
         totals[icam] = totals[icam] + float(value)
      endfor
      nextline = nextline + colsep

      if (iexp EQ 0) then textout = rowsep + nextline $
       else textout = [textout, rowsep + nextline]
   endfor

   textout = [textout, apo_log_tableline(ncams)]

   if (keyword_set(total)) then begin
      nextline = strupcase(flavor) + '-TOTAL' $
       + colsep + fieldname
      for icam=0, ncams-1 do $
       nextline = nextline + colsep + string(totals[icam])
      nextline = nextline + colsep + string(min(totals))
      textout = [textout, rowsep + nextline]

      textout = [textout, apo_log_tableline(ncams)]
   endif

   return, textout
end

;------------------------------------------------------------------------------
pro apo_log2html, logfile, htmlfile

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
      fields = ['XMIN', 'XMAX']
      textout = [ textout, apo_log_fields(pflats, fields) ]

      ;----------
      ; Find the first (and presumably only) FLAVOR=arc for each camera

      parcs = replicate(ptr_new(), ncams)
      for icam=0, ncams-1 do begin
         ii = where(plate EQ thisplate AND flavor EQ 'arc' $
                    AND camera EQ camnames[icam])
         if (ii[0] NE -1) then parcs[icam] = pstruct[ii[0]]
      endfor
      fields = ['WAVEMIN', 'WAVEMAX', 'BESTCORR', 'NLAMPS']
      textout = [ textout, apo_log_fields(parcs, fields) ]

      ;----------
      ; Find all science exposures

      ; Now find all unique science exposure numbers for this plate
      ii = where(plate EQ thisplate AND flavor EQ 'science')
      allexp = expnum[ii[ uniq(expnum[ii], sort(expnum[ii])) ]]
      nexp = n_elements(allexp)

      pscience = replicate(ptr_new(), ncams, nexp)
      for icam=0, ncams-1 do begin
         for iexp=0, nexp-1 do begin
            ii = where(plate EQ thisplate AND flavor EQ 'science' $
                       AND camera EQ camnames[icam] AND expnum EQ allexp[iexp])
            if (ii[0] NE -1) then pscience[icam, iexp] = pstruct[ii[0]]
         endfor
      endfor
      textout = [ textout, apo_log_science(pscience, 'SKYPERSEC') ]
      textout = [ textout, apo_log_science(pscience, 'SN2', /total) ]

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
