;+
; NAME:
;   guidermonfile
;
; PURPOSE:
;   Create guiderMon file using proc-gimg files.
;
; CALLING SEQUENCE:
;   guidermonfile, [ mjd=, mjstart=, mjend=, /clobber ]
;
; INPUTS:
;   mjd        - MJD(s)
;   mjstart    - Starting MJD
;   mjend      - Ending MJD
;
; OPTIONAL INPUTS:
;   clobber    - If set, then over-write existing file, otherwise do not.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The GUIDE_DIR, SPECLOG_DIR environment variables must be set.
;   Data is input from
;     $GUIDE_DIR/$MJD/proc-gimg-????.fits
;   Output file is
;     $SPECLOG_DIR/$MJD/guiderMon-$MJD.par
;
; EXAMPLES:
;
; BUGS:
;
; REVISION HISTORY:
;   13-Jan-2010  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro guidermonfile, mjd=mjd, mjstart=mjstart, mjend=mjend, clobber=clobber

   speclog_dir = getenv('SPECLOG_DIR')
   if (NOT keyword_set(speclog_dir)) then $
    message, 'Must set environment variable SPECLOG_DIR'
   guide_dir = getenv('GUIDE_DIR')
   if (NOT keyword_set(guide_dir)) then $
    message, 'Must set environment variable GUIDE_DIR'

   ;----------
   ; Create a list of the MJD directories (as strings)

   mjdlist = get_mjd_dir(speclog_dir, mjd=mjd, mjstart=mjstart, mjend=mjend)
   nmjd = n_elements(mjdlist)

   if (nmjd GT 1) then begin
      splog, 'Number of MJDs = ', nmjd
      for i=0L, nmjd-1L do guidermonfile, mjd=mjdlist[i], clobber=clobber
      return
   endif

   outfile = filepath('guiderMon-'+mjdlist+'.par', root_dir=speclog_dir, $
    subdir=mjdlist)
   if (keyword_set(clobber) EQ 0) then begin
      oldfile = findfile(outfile, count=ct)
      if (ct GT 0) then begin
         splog, 'Skipping existing file '+outfile
         return
      endif
   endif

   infile = filepath('proc-gimg-????.fits', root_dir=guide_dir, $
    subdir=mjdlist)
   infile = findfile(infile, count=nfile)

   for i=0L, nfile-1L do begin
      hdr = headfits(infile[i])
      gdat1 = mrdfits(infile[i], 6, /silent)
      nfib1 = n_elements(gdat1)
      if (keyword_set(gdat1)) then begin
         if (keyword_set(outdat) EQ 0) then begin
            outdat1 = gdat1[0]
            if (tag_exist(gdat1, 'timestamp') EQ 0) then $
             outdat1 = struct_addtags({timestamp: 0LL}, outdat1)
            if (tag_exist(gdat1, 'fiberid') EQ 0) then $
             outdat1 = struct_addtags({fiberid: 0L}, outdat1)
            struct_assign, {junk: 0}, outdat1
            outdat = replicate(outdat1, nfile*17)
         endif
         index_to=17*i+lindgen(nfib1)
         copy_struct_inx, gdat1, outdat, index_to=index_to
         if (tag_exist(gdat1, 'timestamp') EQ 0) then $
          outdat[index_to].timestamp = $
           long64(date_conv(sxpar(hdr,'DATE-OBS'), 'MODIFIED') * 24.D * 3600.D)
         if (tag_exist(gdat1, 'fiberid') EQ 0) then $
          outdat[index_to].fiberid = 1 + indgen(nfib1)
      endif
   endfor

   ; Trim to non-zero entries
   if (NOT keyword_set(outdat)) then return
   indx = where(outdat.timestamp GT 0, ct)
   if (ct EQ 0) then return
   outdat = outdat[indx]

   outhdr = 'mjd '+mjdlist

   ; Create the output directory if it does not exist already
   junk = fileandpath(outfile, path=outdir)
   if (keyword_set(find_all_dir(outdir)) EQ 0) then begin
      spawn, 'mkdir -p '+outdir
      if (keyword_set(find_all_dir(outdir)) EQ 0) then $
       message, 'Unable to create output directory '+outdir
   endif

   splog, 'Writing file '+outfile
   yanny_write, outfile, ptr_new(outdat), hdr=outhdr, stnames='GUIDEOBJ', $
    /align

   return
end
;------------------------------------------------------------------------------
