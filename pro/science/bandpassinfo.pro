;+
; NAME:
;   bandpassinfo
; PURPOSE:
;   Return information about SDSS bandpasses u,g,r,i,z
; COMMENTS:
;   If band is passed as a string, it gets whitespace-trimmed before use.
; CALLING SEQUENCE:
;   bandpassinfo, band,index=index,name=name,wave=wave,fwhm=fwhm,zero=zero
; INPUTS:
;   band    - a name (ugriz) or index (01234), or a vector of them
; OPTIONAL KEYWORDS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
;   index   - the band's index
;   name    - name
;   wave    - central wavelength in Angstroms
;   fwhm    - width of the bandpass in Angstroms
;   zero    - the flux in Jy of a zero-magnitude source
; BUGS:
;   Returns 0's (u's) or 4's (z's) where the input is wacky.
; PROCEDURES CALLED:
; REVISION HISTORY:
;   2000-Jun-28  Written by Hogg (IAS)
;-
;------------------------------------------------------------------------------
pro bandpassinfo, band,index=index,name=name,wave=wave,fwhm=fwhm,zero=zero

; set up bandpass data
  nband= 5
  bname=  [    'u',    'g',    'r',    'i',    'z']
  bwave=  [ 3000.0, 5000.0, 7000.0, 8000.0, 9000.0]
  bfwhm=  [ 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
  bzero=  [1.0e-23,1.0e-23,1.0e-23,1.0e-23,1.0e-23]

; get size of band input
  nn= n_elements(band)

; check whether input is in terms of indices or names; make index
  if size(band,/type) EQ 7 then begin
    if nn EQ 1 then index= 0 else index= intarr(nn)
    for i=0,nn-1 do begin
      ii= where(bname EQ strtrim(band[i],2),mm)
      if mm EQ 1 then index[i]= ii[0]
    endfor
  endif else begin
    index= ((fix(band) > 0) <4)
  endelse

  name= bname[index]
  wave= bwave[index]
  fwhm= bfwhm[index]
  zero= bzero[index]
  return
end
