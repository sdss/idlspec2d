;+
; NAME:
;   extract_profile
;
; PURPOSE:
;   Extract the flux with profile weighting at centroid positions.
;
; CALLING SEQUENCE:
;   fextract = extract_profile( fimage, invvar, xcen, ycen, 
;                  [ferror=ferror, fscattered=fscattered, fwidth=fwidth, 
;                   sigma=sigma, nPoly=nPoly, maxIter=maxIter,
;                   refit=refit, highrej=highrej, lowrej=lowrej, boxap=boxap])
;
; INPUTS:
;   fimage     - Image[nCol,nRow]
;   invvar     - Inverse Variance[nCol,nRow]
;   xcen       - Initial guesses for X centers[nFibers,nRow]
;   ycen       - Y positions corresponding to "xcen[nFibers]" [nRow]
;
; OPTIONAL KEYWORDS:
;   ferror     - error array (-1 if bad pixel)
;   fscattered - scattered light contribution 
;   fwidth     - final profile width used
;
;   sigma      - sigma of gaussian profile; default to 1.0
;   nPoly      - order of chebyshev scattered light background; default to 5
;   maxIter    - maximum number of profile fitting iterations; default to 5
;   refit      - order of chebyshev to fit to profile widths vs. row,
;			default is 0, which does no refitting
;   highrej    - positive sigma deviation to be rejected (default 5.0)
;   lowrej     - negative sigma deviation to be rejected (default 5.0)
;   boxap      - boxcar aperture size in pixels 
;			default is 0, which does no boxcar extraction
;
; OUTPUTS:
;   fextract   -  Extracted flux at each position specified by (xcen, ycen)
;                  [nFibers,nRow]
;		2: Scattered light estimate
;		3: Final profile width 
;		4: Boxcar extracted flux
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to extract_profile.c
;
; REVISION HISTORY:
;   24-Mar-1999  David Schlegel, Princeton.
;   24-Jun-1999  Stolen and modified by Scott Burles, Chicago.
;-
;------------------------------------------------------------------------------
function extract_profile, fimage, invvar, xcen, ycen, ferror=ferror, $
     fscattered=fscattered, fwidth=fwidth, sigma=sigma, nPoly=nPoly,$
     maxIter=maxIter, refit=refit, highrej=highrej, lowrej=lowrej, boxap=boxap,$
     ymod=ymod

   ; Need 4 parameters
   if (N_params() LT 4) then begin
      print, 'Syntax - fextract = extract_profile(fimage, invvar, xcen, ycen, [ferror=ferror,' 
      print, ' fscattered=fscattered, fwidth=fwidth, sigma=sigma, nPoly=nPoly, ymod=ymod] )'
      return, -1
   endif
   if (NOT keyword_set(sigma)) then begin
	sigma = xcen*0.0 + 1.0
   endif else if ((size(sigma))[0] LT 1) then begin
	sigma1 = sigma
	sigma = xcen*0.0 + sigma1
   endif

   if (NOT keyword_set(nPoly)) then nPoly = 5
   if (NOT keyword_set(maxIter)) then maxIter = 5
   if (NOT keyword_set(refit)) then refit = 0
   if (NOT keyword_set(highrej)) then highrej = 5.0
   if (NOT keyword_set(lowrej)) then lowrej = 5.0
   if (NOT keyword_set(boxap)) then boxap = 0.0

   if (N_elements(fimage) NE N_elements(invvar)) then $
    message, 'Number of elements in FIMAGE and INVVAR must be equal'

;
;	Check XCEN and YCEN to see if the dimensions are consistent
;	and applicable
;
;
   ndim = size(fimage,/n_dim)
   if (ndim NE 2) then $
    message, 'FIMAGE must be 2-dimensional'
   dims = size(fimage,/dim)
   nx = dims[0]
   ny = dims[1]
   nRowExtract = 1
   if ((size(xcen))[0] EQ 1) then begin
	nTrace = (size(xcen))[1]
	if ((size(ycen))[0] EQ 0) then begin
	  y1 = ycen
	  ycen = fltarr(1)
	  ycen[0] = y1
	endif else begin
	  y1 = ycen[0]
	  ycen = fltarr(1)
	  ycen[0] = y1
	endelse
;
;	Check that xcen is sorted in increasing order
;	with separations of at 3 pixels.
;
   endif else if ((size(xcen))[0] EQ 2) then begin
      nTrace = (size(xcen))[1]     
      nRowExtract = (size(xcen))[2]     
      if (N_elements(ycen) NE nRowExtract) then $
        message, 'Number of rows in XCEN does not equal nelements in YCEN'
      ycenf = reform(ycen,nRowExtract)

      check = where(xcen[0:nTrace-1,*] GE xcen[1:nTrace-2,*] - 3,count)
;      if(count GT 0) then $
;       message, 'XCEN is not sorted or not separated by greater than 3 pixels.'
   endif


   nPoly = LONG(nPoly)
   maxIter = LONG(maxIter)
   refit = LONG(refit)
		
;   if (min(xcen) LT 0 OR max(xcen) GT nx) then $
;    message, 'XCEN contains values out of range'
   if (min(ycen) LT 0 OR max(ycen) GT ny) then $
    message, 'YCEN contains values out of range'

   fextract = fltarr(nTrace,nRowExtract)
   ferror = fltarr(nTrace,nRowExtract)
   fwidth = fltarr(nTrace,nRowExtract)

   ymod = fimage*0.0
   fscattered = fimage*0.0

   soname = filepath('libspec2d.'+idlutils_so_ext(), $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='lib')

   result = call_external(soname, 'extract_profile',$
    nx, ny, float(fimage), float(invvar), float (ymod), nTrace, $
    LONG(nRowExtract), float(xcen), LONG(ycen), float(sigma), float(fextract),$
    float(ferror), float(fscattered), float(fwidth),$
    LONG(nPoly), LONG(maxIter), LONG(refit), float(highrej), float(lowrej), $
    float(boxap))
    

   return, fextract
end
;------------------------------------------------------------------------------
