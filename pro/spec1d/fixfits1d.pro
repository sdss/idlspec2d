;------------------------------------------------------------------------------
pro fixfits1d

   vers = 'v4_1_2'

   filename = findfile('0*/spZ*.fits', count=nfile)

   for ifile=0, nfile-1 do begin
      print, 'Modifying ' + filename[ifile]

      hdr = headfits(filename[ifile])

      sxaddpar, hdr, 'NAXIS', 0
      sxdelpar, hdr, 'NAXIS1'
      sxdelpar, hdr, 'NAXIS2'
      sxaddpar, hdr, 'EXTEND', 'T', after='NAXIS'
      sxdelpar, hdr, 'WCSDIM'

      sxaddpar, hdr, 'VERSUTIL', vers
      sxaddpar, hdr, 'VERSREAD', vers
      sxaddpar, hdr, 'VERS2D', vers
      sxaddpar, hdr, 'VERSCOMB', vers
      sxaddpar, hdr, 'VERS1D', vers, $
       'Version of idlspec2d for 1D reduction', after='VERSCOMB'

      modfits, filename[ifile], 0, hdr
   endfor
end

;------------------------------------------------------------------------------
pro fixfits2d

   vers = 'v4_1_2'

   filename = findfile('0*/spPlate*.fits', count=nfile)

   for ifile=0, nfile-1 do begin
      ;----------
      ; Fix the first HDU

      print, 'Modifying ' + filename[ifile]

      hdr = headfits(filename[ifile])
      nhead = n_elements(hdr)

      sxdelpar, hdr, 'WCSDIM'
      sxdelpar, hdr, 'FOCUS'

      sxaddpar, hdr, 'VERSUTIL', vers
      sxaddpar, hdr, 'VERSREAD', vers
      sxaddpar, hdr, 'VERS2D', vers
      sxaddpar, hdr, 'VERSCOMB', vers

      ; Make certain we have exactly the same number of header cards
      for iadd=0, nhead-n_elements(hdr)-1 do $
       sxaddpar, hdr, 'HISTORY', ' '

      modfits, filename[ifile], 0, hdr

      ;----------
      ; Fix the other HDU's

      for ihdu=1, 4 do begin
         nexthdr = headfits(filename[ifile], exten=ihdu)
         sxaddpar, nexthdr, 'WAT0_001', sxpar(hdr, 'WAT0_001')
         sxaddpar, nexthdr, 'WAT1_001', sxpar(hdr, 'WAT1_001')
         sxaddpar, nexthdr, 'CRVAL1', sxpar(hdr, 'CRVAL1')
         sxaddpar, nexthdr, 'CD1_1', sxpar(hdr, 'CD1_1')
         sxaddpar, nexthdr, 'CRPIX1', sxpar(hdr, 'CRPIX1')
         sxaddpar, nexthdr, 'CTYPE1', sxpar(hdr, 'CTYPE1')
         sxaddpar, nexthdr, 'DC-FLAG', sxpar(hdr, 'DC-FLAG')

         modfits, filename[ifile], 0, nexthdr, exten_no=ihdu
      endfor

   endfor
end
;------------------------------------------------------------------------------
