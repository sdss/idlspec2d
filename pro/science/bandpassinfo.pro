;+
; NAME:
;   bandpassinfo
; PURPOSE:
;   Return information about SDSS bandpasses u,g,r,i,z.
; COMMENTS:
;   This routine can be used to convert filternames into numbers and back.
;   If band is passed as a string, it gets whitespace-trimmed before use.
; CALLING SEQUENCE:
;   bandpassinfo, band,index=index,name=name,wave=wave,fwhm=fwhm,zero=zero
; INPUTS:
;   band    - a name (ugriz) or index (01234), or a vector of them
; OPTIONAL KEYWORDS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
;   index     - the band's index number (0-4)
;   name      - name ('u' through 'z')
;   wave      - central wavelength in Angstroms
;   fwhm      - width of the bandpass in Angstroms
;   zero      - the flux in Jy of a zero-magnitude source
;   blueindex - the index (0-3) of the blue bandpass of the k-correction color
;   redindex  - the index (1-4) of the red bandpass of the k-correction color
;   colorname - the name of the color ('(u-g)' through '(i-z)')
; BUGS:
;   Returns 0's (u's) or 4's (z's) where the input is wacky.
;   Computes even unnecessary things (this could be fixed with some calls to
;     keyword_set().
; PROCEDURES CALLED:
; REVISION HISTORY:
;   2000-Jun-28  Written by Hogg (IAS)
;-
;------------------------------------------------------------------------------
pro bandpassinfo, band,index=index,name=name,wave=wave,fwhm=fwhm,zero=zero, $
                  blueindex=blueindex,redindex=redindex,colorname=colorname

; set up bandpass data
  nband= 5
  bname=  [    'u',    'g',    'r',    'i',    'z']
  bwave=  [ 3543.0, 4770.0, 6231.0, 7625.0, 9134.0]
  bfwhm=  [  567.0, 1387.0, 1373.0, 1526.0,  950.0]
  bzero=  [1.0e-23,1.0e-23,1.0e-23,1.0e-23,1.0e-23]
  bblui=  [      0,      1,      2,      3,      3]
  bredi=  [      1,      2,      3,      4,      4]
  bcolorname= '('+bname[bblui]+'-'+bname[bredi]+')'

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

; output and return
  name= bname[index]
  wave= bwave[index]
  fwhm= bfwhm[index]
  zero= bzero[index]
  blueindex= bblui[index]
  redindex= bredi[index]
  colorname= bcolorname[index]
  return
end
