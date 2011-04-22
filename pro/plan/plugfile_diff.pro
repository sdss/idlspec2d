;+
; NAME:
;   plugfile_diff
;
; PURPOSE:
;   Find instances where the plPlugMapM file in the speclog file differs
;   from what was used in the reductions.
;
; CALLING SEQUENCE:
;   plugfile_diff, [ run2d=, run1d= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   run2d       - Limit the platelist to this version of the 2D reductions
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only plates that completed reductions are considered.
;   Generate a FITS file 'plugfile_diff.fits' containing all of objects
;   that have had their plugging changed.  File contains:
;     HDU #1 - ZANS structure for those objects
;     HDU #2 - PLUGMAP info used for reductions
;     HDU #3 - PLUGMAP info currently in speclog product
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   platelist
;   readspec
;   splog
;
; REVISION HISTORY:
;   21-Apr-2011  Written by David Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro plugfile_diff, run2d=run2d, run1d=run1d

   platelist, plist=plist
   if (NOT keyword_set(plist)) then begin
      splog, 'No plates found'
      return
   endif
;   qtrim = strmatch(plist.platequality,'good*')
   qtrim = strmatch(plist.status1d,'Done*')
   if (keyword_set(run2d)) then $
    qtrim *= strtrim(plist.run2d) EQ run2d
   if (keyword_set(run1d)) then $
    qtrim *= strtrim(plist.run1d) EQ run1d
   itrim = where(qtrim, nplate)
   if (nplate EQ 0) then begin
      splog, 'No plates found'
      return
   endif
   plist = plist[itrim]

   for i=0L, nplate-1L do begin
      print, 'File ', i, ' of ', nplate
      readspec, plist[i].plate, mjd=plist[i].mjd, plug=plug1, objhdr=hdr, $
       zans=zans1, run2d=plist[i].run2d, run1d=plist[i].run1d, /silent
      if (keyword_set(plug1[0]) * keyword_set(hdr) EQ 0) then $
       message, 'Error reading plate '+string(plist[i].plate)
      plugfile = 'plPlugMapM-' + sxpar(hdr,'NAME') + '.par'
      thismjd = sxpar(hdr,'MJD')
      plugdir = getenv('SPECLOG_DIR') + '/' $
       + string(thismjd,format='(i5.5)')
      plug2 = readplugmap(plugfile, plugdir=plugdir)

      ; Match the FIBERID from the two lists.  This matters if only
      ; half the fibers exist in some reductions.
      ig1 = where(plug1.fiberid GT 0)
      ig2 = where(plug2.fiberid GT 0)
      match, plug1[ig1].fiberid, plug2[ig2].fiberid, i1, i2
      plug1 = plug1[ig1[i1]]
      zans1 = zans1[ig1[i1]]
      plug2 = plug2[ig2[i2]]

      adist = djs_diff_angle(plug1.ra, plug1.dec, plug2.ra, plug2.dec)
      idiff = where(adist GT 2./3600, ndiff)
      splog, 'Change ', ndiff, ' fibers ', plugfile
      if (ndiff GT 0) then begin
         zanslist = struct_append(zanslist, zans1[idiff])
         badlist1 = struct_append(badlist1, plug1[idiff])
         badlist2 = struct_append(badlist2, plug2[idiff])
      endif
   endfor

   outfile = 'plugfile_diff.fits'
   mwrfits, zanslist, outfile, /create
   mwrfits, badlist1, outfile
   mwrfits, badlist2, outfile

   return
end
;------------------------------------------------------------------------------
