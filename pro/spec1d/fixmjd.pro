; Fix the MJD in the output files spPlate, spZbest, spZall to be consistent
; with the file name.  These are the header keywords MJD in HDU#0 of all
; files, and the MJD element of the structures in the spZbest,spZall files.
;------------------------------------------------------------------------------
pro fixmjd

   ;---------------------------------------------------------------------------
   ; Fix spPlate files

   filename = findfile('0*/spPlate*.fits', count=nfile)

   for ifile=0, nfile-1 do begin
      print, 'Reading ' + filename[ifile]

      mjdname = long( strmid( fileandpath(filename[ifile]), 13, 5) )

      hdr = headfits(filename[ifile])
      mjdhdr = sxpar(hdr, 'MJD')

      if (mjdname NE mjdhdr) then begin
         print, 'MODIFY!'
         sxaddpar, hdr, 'MJD', mjdname
         djs_modfits, filename[ifile], 0, hdr
      endif
   endfor

   ;---------------------------------------------------------------------------
   ; Fix spZbest,spZall files

   filename = findfile('0*/spZbest*.fits', count=nfile)

   for ifile=0, nfile-1 do begin
      print, 'Reading ' + filename[ifile]

      mjdname = long( strmid( fileandpath(filename[ifile]), 13, 5) )

      hdr = headfits(filename[ifile])
      mjdhdr = sxpar(hdr, 'MJD')

      if (mjdname NE mjdhdr) then begin
         print, 'MODIFY!'
         ; Modify spZbest file
         sxaddpar, hdr, 'MJD', mjdname
         djs_modfits, filename[ifile], 0, hdr
         zans = mrdfits(filename[ifile], 1, hdr1)
         zans.mjd = mjdname
         djs_modfits, filename[ifile], zans, exten_no=1

         ; Modify spZall file
         zallfile = repstr(filename[ifile], 'spZbest', 'spZall')
         hdr = headfits(zallfile)
         sxaddpar, hdr, 'MJD', mjdname
         djs_modfits, zallfile, 0, hdr
         zans = mrdfits(zallfile, 1, hdr1)
         zans.mjd = mjdname
         djs_modfits, zallfile, zans, exten_no=1
      endif
   endfor
end
;------------------------------------------------------------------------------
