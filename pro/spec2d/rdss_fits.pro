;+
; NAME:
;   rdss_fits
;
; PURPOSE:
;   Read a FITS file into IDL data and header variables
;
; CALLING SEQUENCE:
;   image = rdss_fits( filename, [ hdr ] )
;
; INPUTS:
;   filename   - Scalar string containing the name of the FITS file  
;                (including extension) to be read.   If the filename has
;                a *.gz extension, it will be treated as a gzip compressed
;                file.   If it has a .Z extension, it will be treated as a
;                Unix compressed file.
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   image      - FITS data array constructed from designated record.
;                If the specified file was not found, then return -1.
;
; OPTIONAL OUTPUTS:
;   hdr        - String array containing the header from the FITS file.
;
; COMMENTS:
;   This routine will read a simple FITS image, or convert a non-standard
;   SDSS image that uses unsigned 16-bit integers.  One can pass any other
;   keywords expected by READFITS().
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   readfits()
;   sxpar()
;
; REVISION HISTORY:
;   13-May-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function rdss_fits, filename, hdr, _EXTRA=KeywordsForReadfits

   ; Need at least 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - image = rdss_fits( filename, [ hdr ] )'
      return, -1
   endif

   ; Read the image and header
   image = mrdfits(filename, 0, hdr, _EXTRA=KeywordsForReadfits)

   ; Test to see if the image is stored in the non-standard unsigned integer
   ; format of SDSS
   qsimple = strpos( (str_sep(hdr[0],'/'))[0], 'F', 8) ; GE 0 if SIMPLE=F
   result = sxpar(hdr, 'UNSIGNED', count=ct1)
   if (qsimple GE 0 AND ct1 NE 0) then begin
print,'converting from U16...'

      bitpix = sxpar(hdr, 'BITPIX')
      sxdelpar, hdr, 'UNSIGNED' ; Remove UNSIGNED keyword
      sxdelpar, hdr, ' '        ; Remove blank header cards

      ; Convert from unsigned 16-bit integers to floats
      ineg = where(image LT 0)
      image = temporary(image) + 0.0
      image[ineg] = 2.0^bitpix + image[ineg]
ineg = 0

   endif

   return, image
end
