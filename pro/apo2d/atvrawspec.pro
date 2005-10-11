;+
; NAME:
;   atvrawspec
;
; PURPOSE:
;   Display a 2D spectroscopic image with bad columns marked in red.
;
; CALLING SEQUENCE:
;   atvrawspec, filename, _EXTRA=KeywordsForATV
;
; INPUTS:
;   filename   - Either a raw SDSS file name, or a 2D pixel flat or 2D bias
;
; OPTIONAL KEYWORDS:
;   _EXTRA     - Extra keywords for ATV, such as MIN,MAX
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   Display a 2D pixel flat with a reasonable display stretch:
;   IDL> atvrawspec, 'pixflat-52069-b1.fits', min=0.9, max=1.1
;
;   Display a 2D bias with a reasonable display stretch:
;   IDL> atvrawspec, 'pixbias-52069-b1.fits', min=-5, max=10
;
; BUGS:
;
; PROCEDURES CALLED:
;   atv
;   atvplot
;   djs_filepath()
;   findopfile()
;   headfits()
;   sdssproc
;   yanny_free
;   yanny_read
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/opBC*par
;
; REVISION HISTORY:
;   26-Feb-2002  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro atvrawspec, filename, _EXTRA=KeywordsForATV

   if (n_params() LT 1) then begin
      print, 'Syntax - atvrawspec, filename'
      return
   endif

   badcolor = 'red'

   hdr = headfits(filename)
   naxis1 = sxpar(hdr, 'NAXIS1')
   naxis2 = sxpar(hdr, 'NAXIS2')
   if (naxis1 EQ 2128 AND naxis2 EQ 2069) then begin
      sdssproc, filename, image
   endif else if (naxis1 EQ 2048 AND naxis2 EQ 2048) then begin
      image = mrdfits(filename)
   endif else begin
      splog, 'Warning: Unknown file dimensions ', naxis1, naxis2
      return
   endelse

   ;-----------
   ; Display the image with ATV

   atv, image, _EXTRA=KeywordsForATV

   ;-----------
   ; Read the opBC (bad column) file

   config_dir = filepath('', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
   mjd = sxpar(hdr, 'MJD')
   bcfile = findopfile('opBC*par', mjd, config_dir, /abort_notfound)
   if (NOT keyword_set(bcfile)) then begin
      splog, 'No opBC file found for MJD=', mjd
      return
   endif
   bcfile = djs_filepath(bcfile, root_dir=config_dir)
   yanny_read, bcfile, pdata
   bc = *pdata[0]
   if (NOT keyword_set(bc)) then begin
      splog, ' Warning: Bad column file missing or empty ' + bcfile
      return
   endif
   yanny_free, pdata

   ;-----------
   ; Selct entries from the bad column file that match this CCD

   cameras = strtrim( sxpar(hdr, 'CAMERAS'), 2 )
   camnames = ['b1', 'r2', 'b2', 'r1']
   indx = where(cameras EQ camnames)
   if (indx[0] EQ -1) then begin
      splog, 'Warning: CAMERAS in image header does not match any valid value'
      return
   endif
   camcol = indx[0] + 1
   camrow = 0
   ibc = where(bc.camrow EQ camrow AND bc.camcol EQ camcol, nbc)
   bc = bc[ibc]

   ;-----------
   ; Overplot bad columns in red

   for ibc=0, nbc-1 do begin
      if (bc[ibc].dfncol EQ 1 AND bc[ibc].dfnrow EQ 1) then begin
         atvplot, bc[ibc].dfcol0, bc[ibc].dfrow0, psym=1, color=badcolor
      endif else if (bc[ibc].dfncol EQ 1 OR bc[ibc].dfnrow EQ 1) then begin
         atvplot, bc[ibc].dfcol0+[0,bc[ibc].dfncol-1], $
          bc[ibc].dfrow0+[0,bc[ibc].dfnrow-1], $
          color=badcolor
      endif else begin
         atvplot, bc[ibc].dfcol0+[0,bc[ibc].dfncol-1,bc[ibc].dfncol-1,0,0], $
          bc[ibc].dfrow0+[0,0,bc[ibc].dfnrow-1,bc[ibc].dfnrow-1,0], $
          color=badcolor
         if (bc[ibc].dfncol GE 3) then $
          for j=1, bc[ibc].dfncol-2 do $
           atvplot, bc[ibc].dfcol0+[j,j], bc[ibc].dfrow0+[0,bc[ibc].dfnrow-1], $
            color=badcolor
      endelse
   endfor

   return
end
;------------------------------------------------------------------------------
