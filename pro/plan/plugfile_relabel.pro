;+
; NAME:
;   plugfile_relabel
;
; PURPOSE:
;   Relabel the OBJTYPE's for specified objects in all plPlugMapM files.
;
; CALLING SEQUENCE:
;   plufile_relabel, [ infile, objtype= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   infile     - ASCII file with plate, mjd, fiberid, as written
;                by BADFSTARS (for example); default to 'badfstars.log'
;   objtype    - New OBJTYPE for these objects; default to 'SERENDIPITY_MANUAL'
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This procedure was written to get rid of objects from the plPlugMapM
;   files that were incorrectly being used as F star calibrators
;   (either SPECTROPHOTO_STD or REDDEN_STD targets).  For a specified list
;   of objects, this proc changes the OBJTYPE to 'SERENDIPITY_MANUAL'.
;   Note that the target flags are not changed, since they are not used
;   anyway by Spectro-2D.  This means that we will still have the 2 or 32
;   values set in SECTARGET.
;
;   For the input list of objects, I use READSPEC to find the coordinates
;   of the object, and then the object is changed in all the plPlugMapM
;   files for that plate.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   readspec
;   splog
;   yanny_read
;   yanny_write
;
; REVISION HISTORY:
;   06-Feb-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro plugfile_relabel, infile, objtype=objtype

   if (NOT keyword_set(infile)) then infile = 'badfstars.log'
   if (NOT keyword_set(objtype)) then objtype = 'SERENDIPITY_MANUAL'

   ;----------
   ; Insist on modifying a CVS version of the speclog product, and
   ; not some tagged version.

   spawn, 'speclog_version', logvers, /noshell
   print, 'SPECLOG_VERSION= ' + logvers
   if (strmatch(logvers, 'NOCVS*') EQ 0) then begin
      print, 'The speclog version is not a CVS version, so quitting'
      return
   endif
   logdir = getenv('SPECLOG_DIR')

   ;----------
   ; Read the input file

   readcol, infile, plate, mjd, fiberid, format='(L,L,L)'
   if (NOT keyword_set(plate)) then begin
      print, 'No objects read from file ' + infile
      return
   endif

   ;----------
   ; Loop over each object

   for iobj=0L, n_elements(plate)-1 do begin
      ; Get the RA,DEC for this object
      readspec, plate[iobj], mjd=mjd[iobj], fiberid[iobj], plug=thisplug
      thisra = thisplug.ra
      thisdec = thisplug.dec

      ; Find all plPlugMapM files for this plate
      platestr = string(plate[iobj], format='(I4.4)')
      mapfile = 'plPlugMapM-' + platestr + '-*.par'
      mapfile = findfile(filepath(mapfile, root_dir=logdir, subdir='?????'), $
       count=nfile)
      splog, 'Object #', iobj, ' has ', nfile, ' plug-map files'

      ; Loop over each plPlugMapM file, changing this object as appropriate
      for ifile=0L, nfile-1 do begin
         yanny_read, mapfile[ifile], pp, hdr=hdr, enums=enums, $
          structs=structs, stnames=stnames, /anonymous
         ip = (where(stnames EQ 'PLUGMAPOBJ'))[0]
         if (ip EQ -1) then $
          message, 'Invalid plug-map file ' + mapfile[ifile]
         adiff = djs_diff_angle((*(pp[ip])).ra, (*(pp[ip])).dec, $
          thisplug.ra, thisplug.dec)
         imatch = where(adiff LT 1./3600 $
          AND strmatch((*(pp[ip])).holetype,'OBJECT*'), nmatch)

         ; Change the OBJTYPE for any matching objects, and over-write the file
         if (nmatch GT 0) then begin
            (*(pp[ip]))[imatch].objtype = objtype
            splog, 'Modifying ' + mapfile[ifile]
            yanny_write, mapfile[ifile], pp, hdr=hdr, enums=enums, $
             structs=structs, stnames=stnames
         endif
      endfor
   endfor

   return
end
;------------------------------------------------------------------------------
