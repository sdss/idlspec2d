;+
; NAME:
;   listexptime
;
; PURPOSE:
;   List the exposure times of all (reduced) spectroscopic frames.
;
; CALLING SEQUENCE:
;   listexptime
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The exposure times of all files $BOSS_SPECTRO_REDUX/????/spFrame-b1-????????.fits
;   are logged to the file 'listexp.log'.
;   This proc just looks at the b1 exposures to save time.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   headfits()
;   splog
;
; REVISION HISTORY:
;   03-Sep-2003  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro listexptime

   spawn, '\ls -d '+getenv('BOSS_SPECTRO_REDUX') + '/????', mjdlist
   nmjd = n_elements(mjdlist)

   for imjd=0, nmjd-1 do begin
      thisname = filepath('spFrame-b1-????????.fits*', $
       root_dir=mjdlist[imjd])
      thisname = findfile(thisname, count=ct)
      if (ct GT 0) then begin
         if (NOT keyword_set(framename)) then framename = thisname $
          else framename = [framename, thisname]
      endif
   endfor

   splog, filename='listexp.log'

   nfile = n_elements(framename)
   exptime = fltarr(nfile)
   plate = lonarr(nfile)
   mjd = lonarr(nfile)
   for ifile=0, nfile-1 do begin
      print, 'Reading file ', ifile, ' of ', nfile, ': ' + framename[ifile]
      hdr = headfits(framename[ifile])
      exptime[ifile] = sxpar(hdr, 'EXPTIME')
      plate[ifile] = sxpar(hdr, 'PLATEID')
      mjd[ifile] = sxpar(hdr, 'MJD')
      splog, fileandpath(framename[ifile]), plate[ifile], $
       mjd[ifile], exptime[ifile], /noname
   endfor

   splog, /close

   return
end
;------------------------------------------------------------------------------
