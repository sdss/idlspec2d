;+
; NAME:
;   plugfile_unmapped
;
; PURPOSE:
;   List which fibers in a plPlugMapM file are unmapped, and the info
;   that might help decide where those fibers should be.
;
; CALLING SEQUENCE:
;   plugfile_unmapped, filename
;
; INPUTS:
;   filename    - Look for a plPlugMapM file in the speclog directory
;                 matching $SPECLOG_DIR/?????/FILENAME.
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   There are two lists of information provided: First, the list of fiber
;   numbers that are not plugged, and the RA,DEC ranges of the other 19 fibers
;   in that bundle on this plate.  Second, a list of the object information
;   for the unplugged objects, including their RA,DEC position and the fiber
;   numbers of the nearest plugged fibers.
;
;   Presumably, those unplugged objects will be plugged by one of the fibers
;   whose neighboring fiber numbers are nearby.  This isn't necessarily true,
;   for example if the pluggers simply skipped one fiber entirely while
;   plugging, and the unplugged object is in a completely un-reachable
;   position elsewhere on the plate.
;
;   A companion procedure to this is FIND_UNPLUGGED.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_diff_angle()
;   splog
;   yanny_readone()
;
; REVISION HISTORY:
;   10-Sep-2003  Written by David Schlegel, Princeton (not checked in then)
;-
;------------------------------------------------------------------------------
pro plugfile_unmapped, plugfile

   topdir = getenv('SPECLOG_DIR')
   filename = findfile(filepath(plugfile, root_dir=topdir, $
    subdir='?????'))
   filename = filename[0] ; Select the first file that matches
   if (NOT keyword_set(filename)) then begin
      print, 'File not found'
      return
   endif
   plug = yanny_readone(filename)

   plug = plug[where(strmatch(plug.holetype,'OBJECT*'))]
   qgood = lonarr(640)
   indx = where(plug.fiberid GE 1)
   qgood[plug[indx].fiberid-1] = 1
   ibad = where(qgood EQ 0, nbad)
   if (nbad EQ 0) then begin
      splog, 'No unmapped fibers'
   endif else begin
      splog, 'Unmapped fiber IDs: ', ibad+1
      for j=0L, nbad-1 do begin
         bundlenum = fix(ibad[j]/20)
         inbundle = where(plug.fiberid-1 GE bundlenum*20 $
          AND plug.fiberid LE (bundlenum+1)*20)
         rarange = minmax(plug[inbundle].ra)
         ; Fix wrapping of RA about 360 deg.
         if ((rarange[1] - rarange[0]) GT 180) then $
          rarange = [min(plug[inbundle].ra + (plug[inbundle].ra LE 180)*360), $
                     max(plug[inbundle].ra - (plug[inbundle].ra GE 180)*360) ]
         decrange = minmax(plug[inbundle].dec)
         splog, 'Fiber# ', ibad[j]+1, $
          '   Bundle at RA=', string(rarange,format='(2f7.2)'), $
          ' DEC=', string(decrange,format='(2f7.2)')
      endfor

      jndx = where(plug.fiberid LE 0, nbad)
      for j=0L, nbad-1 do begin
         ; Find the nearest mapped fibers to this unmapped hole
         adist = djs_diff_angle(plug[jndx[j]].ra, plug[jndx[j]].dec, $
          plug[indx].ra, plug[indx].dec)
         inearest = (plug[indx[sort(adist)]].fiberid)[0:9]
         splog, ' '
         splog, string(plug[jndx[j]].objtype+'         ', format='(a10)'), $
          ' MAG=' + string(plug[jndx[j]].mag, format='(5f6.2)'), $
          '  RA=' + string(plug[jndx[j]].ra, format='(f11.6)'), $
          ', DEC=' + string(plug[jndx[j]].dec, format='(f11.6)')
         splog, 'Nearest mapped fibers=',strtrim(string(inearest),2)
      endfor
   endelse

   return
end
;------------------------------------------------------------------------------
