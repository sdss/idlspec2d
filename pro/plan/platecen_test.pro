;+
; NAME:
;   platecen_test
;
; PURPOSE:
;   Find where RA,DEC and RADEG,DECDEG disagree in the spPlate headers.
;
; CALLING SEQUENCE:
;   platecen_test
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
;   Find all the files matching '$BOSS_SPECTRO_REDUX/*/spPlate*.fits', and test
;   whether the RA,DEC positions disagree by more than 0.1 degrees from
;   RADEC,DECDEG.  This would indicate some inconsistency in the header
;   info that should be fixed.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   headfits()
;   splog
;   sxpar()
;
; REVISION HISTORY:
;   19-Dec-2001  Written by David Schlegel, Princeton (not checked in then)
;-
;------------------------------------------------------------------------------
pro platecen_test

   platefile = findfile(getenv('BOSS_SPECTRO_REDUX')+'/*/spPlate*.fits', $
    count=nfile)

   ra = fltarr(nfile)
   dec = fltarr(nfile)
   radeg = fltarr(nfile)
   decdeg = fltarr(nfile)
   for ifile=0, nfile-1 do begin
      hdr = headfits(platefile[ifile])
      ra[ifile] = sxpar(hdr, 'RA')
      dec[ifile] = sxpar(hdr, 'DEC')
      radeg[ifile] = sxpar(hdr, 'RADEG')
      decdeg[ifile] = sxpar(hdr, 'DECDEG')
   endfor

; Ignore values of RADEG that are less than 0, like -999 ?
;   i = where(radeg LT 0)
;   radeg[i] = ra[i]
;   decdeg[i] = dec[i]

   j = where(abs(ra-radeg) GT 0.1 OR abs(dec-decdeg) GT 0.1)
   for jj=0, n_elements(j)-1 do begin
      splog, fileandpath(platefile[j[jj]]), ra[j[jj]], dec[j[jj]], $
       radeg[j[jj]], decdeg[j[jj]]
   endfor

   return
end
;------------------------------------------------------------------------------
