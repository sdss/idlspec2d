;+
; NAME:
;   sdsshead
;
; PURPOSE:
;   Read in header from raw SDSS file.
;
; CALLING SEQUENCE:
;   hdr = sdsshead(infile, [ indir= ] )
;
; INPUTS:
;   infile     - Raw SDSS file name
;
; OPTIONAL KEYWORDS:
;   indir      - Input directory for INFILE
;   do_lock    - Keyword passed to SDSSPROC
;
; OUTPUTS:
;   hdr        - Processed FITS header
;
; COMMENTS:
;
; BUGS:
;
; PROCEDURES CALLED:
;   sdssproc
;
; REVISION HISTORY:
;   27-May-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

function sdsshead, infile, indir=indir, do_lock=do_lock

   sdssproc, infile, indir=indir, hdr=hdr, do_lock=do_lock

   return, hdr
end
;------------------------------------------------------------------------------
