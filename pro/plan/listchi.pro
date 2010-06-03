;+
; NAME:
;   listchi
;
; PURPOSE:
;   List the extraction chi of all individual spectro frames
;
; CALLING SEQUENCE:
;   listchi
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
;   The exposure times of all files $BOSS_SPECTRO_REDUX/????/spFrame-??-????????.fits
;   are logged to the file 'listchi.log'.  Also, a save-set is written
;   to the file 'listchi.ss'.
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
pro listchi

   spawn, '\ls -d '+getenv('BOSS_SPECTRO_REDUX') + '/????', mjdlist
   nmjd = n_elements(mjdlist)

   for imjd=0, nmjd-1 do begin
      thisname = filepath('spFrame-??-????????.fits*', $
       root_dir=mjdlist[imjd])
      thisname = findfile(thisname, count=ct)
      if (ct GT 0) then begin
         if (NOT keyword_set(framename)) then framename = thisname $
          else framename = [framename, thisname]
      endif
   endfor

   splog, filename='listchi.log'

   nfile = n_elements(framename)
   exptime = fltarr(nfile)
   xchi2 = fltarr(nfile)
   skychi2 = fltarr(nfile)
   xsigma = fltarr(nfile)
   wsigma = fltarr(nfile)
   framesn2 = fltarr(nfile)
   plate = lonarr(nfile)
   cartid = lonarr(nfile)
   mjd = lonarr(nfile)
   flavor = strarr(nfile)
   cameras = strarr(nfile)
   for ifile=0, nfile-1 do begin
      print, 'Reading file ', ifile, ' of ', nfile, ': ' + framename[ifile]
      hdr = headfits(framename[ifile])
      exptime[ifile] = sxpar(hdr, 'EXPTIME')
      xchi2[ifile] = sxpar(hdr, 'XCHI2')
      skychi2[ifile] = sxpar(hdr, 'SKYCHI2')
      xsigma[ifile] = sxpar(hdr, 'XSIGMA')
      wsigma[ifile] = sxpar(hdr, 'WSIGMA')
      framesn2[ifile] = sxpar(hdr, 'FRAMESN2')
      plate[ifile] = sxpar(hdr, 'PLATEID')
      mjd[ifile] = sxpar(hdr, 'MJD')
      cartid[ifile] = sxpar(hdr, 'CARTID')
      flavor[ifile] = sxpar(hdr, 'FLAVOR')
      cameras[ifile] = sxpar(hdr, 'CAMERAS')
      splog, fileandpath(framename[ifile]), plate[ifile], $
       mjd[ifile], exptime[ifile], xchi2[ifile], /noname
   endfor

   splog, /close
   save, 'listchi.ss'

   return
end
;------------------------------------------------------------------------------
