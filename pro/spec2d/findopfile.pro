;+
; NAME:
;   findopfile
;
; PURPOSE:
;   Find the op file corresponding to a specified MJD.
;
; CALLING SEQUENCE:
;   filename = findopfile( expres, mjd, [ indir, /abort_notfound ] )
;
; INPUTS:
;   expres     - Op file names to match which may include any wildcards
;                other other expressions valid for the function FINDFILE().
;   mjd        - Number (MJD) for matching corresponding op file.
;
; OPTIONAL INPUT:
;   indir      - Input directory
;   abort_notfound - If set and no files are found to match the expression,
;                    then abort with the MESSAGE command; otherwise return ''.
;
; OUTPUTS:
;   filename   - Matched file name without path information
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The MJD for the op file is determined from the file name itself
;   by looking at whatever number is immediately after the first '-'.
;   For example, the file 'opBC-52000.par' is interpreted as having
;   MJD=52000.
;
;   The matching op file is the one with the same MJD as that requested,
;   or the last one preceding the requested MJD.  If the requested MJD
;   precedes any for the existing op files, then return the op file
;   with the lowest MJD.  For example, if we have two opBC files
;   'opBC-50000.par' and 'opBC-52000.par', then the first is returned
;   for all MJDs up to 51999, and the latter used for MJD=52000 and later.
;
; EXAMPLES:
;   Find the opBC file in the directory $IDLSPEC2D_DIR/examples
;   that is appropriate for data taken on MJD=52000:
;     IDL> indir = getenv('IDLSPEC2D_DIR')+'/examples'
;     IDL> filename = findopfile('opBC*.par', 52000, indir)
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_filepath()
;   fileandpath()
;   splog
;
; INTERNAL SUPPORT ROUTINES:
;
; DATA FILES:
;
; REVISION HISTORY:
;   27-Feb-2002  Written by Scott Burles & David Schlegel.
;                Broken out from an internal function of SDSSPROC.
;-
;------------------------------------------------------------------------------
function findopfile, expres, mjd, indir, abort_notfound=abort_notfound

   if (n_params() LT 2) then begin
      print, 'Syntax - filename = findopfile( expres, mjd, [ indir, /abort_notfound ] )'
      return, ''
   endif

   files = findfile(djs_filepath(expres, root_dir=indir), count=nfile)
   if (nfile EQ 0) then begin
      if (keyword_set(abort_notfound)) then $
       message, 'Cannot find opFile '+expres $
      else $
       return, ''
   endif

   mjdlist = lonarr(nfile)
   for i=0, nfile-1 do begin
      thisfile = fileandpath(files[i])
; Change below to select which file to use based upon the MJD in
; the file name rather than the MJD in the header of a Yanny param file.
;      if (strmid(thisfile, strpos(thisfile, '.')) EQ '.par') then begin
;         ; Get the MJD from the header of the Yanny file
;         yanny_read, files[i], hdr=hdr
;         mjdlist[i] = long(yanny_par(hdr, 'mjd'))
;      endif else begin
         ; Get the MJD from the file name; use the 5 digits after
         ; the first minus sign.  Prepend the '0' below to avoid triggering
         ; an IDL "Type conversion error".
         mjdlist[i] = long( '0'+strmid(thisfile, strpos(thisfile,'-')+1,5) )
;      endelse
   endfor

   ; Sort the files by MJD in descending order
   isort = reverse(sort(mjdlist))
   mjdlist = mjdlist[isort]
   files = files[isort]

   ; Select the op file with the same MJD or the closest MJD preceding
   ; this one.
   ibest = (where(mjdlist LE mjd))[0]
   if (ibest[0] EQ -1) then begin
      splog, 'WARNING: No ' + expres + ' op files appear to have an MJD <= ' $
       + strtrim(string(mjd),2)
      ibest = 0
   endif
   selectfile = fileandpath(files[ibest])
   splog, 'Selecting op file ' + selectfile

   return, selectfile
end
;------------------------------------------------------------------------------
