;+
; NAME:
;   spx_to_sdss
;
; PURPOSE:
;   Convert an spX-extracted spectra file to an SDSS spPlate format
;
; CALLING SEQUENCE:
;   spx_to_sdss, filename, [ platefile=, outfile= ]
;
; INPUTS:
;   filename   - spX spectra file, where HDU #0 is the flux [NPIX,NFIBER],
;                HDU #1 is the inverse variance [NPIX,NFIBER],
;                and HDU #2 contains wavelengths [NPIX]
;   platefile  - Name of spPlate file from which copy HDU#2 through #5;
;                default to 'spPlate-$PLATE-$MJD' using the PLATE and MJD
;                from the FILENAME header; necessary to specify this if
;                the final spPlate file for these data is from a later MJD
;
; OPTIONAL INPUTS:
;   outfile    - Output file name; default to spPPlate-$PLATE-$MJD.fits
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The input spectra are interpolated to the log-wavelength mapping
;   of SDSS spectra using a B-spline to the flux and the variance arrays.
;
; EXAMPLES:
;   Convert an spX file and run the SDSS redshifting code:
;     IDL> spx_to_sdss,'spXvfsc-r2-00104772.fits', $
;           platefile='spPlate-3647-55178.fits', $
;           outfile='spPPlate-3647-55178.fits'
;     IDL> spreduce1d,'spPPlate-3647-55178.fits'
;
; BUGS:
;   The wavelength scale is assumed to be vacuum barycentric.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   21-Apr-2014  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro spx_to_sdss, filename, platefile=platefile1, outfile=outfile1

   if (size(filename,/tname) NE 'STRING') then begin
      splog, 'Must specify FILENAME'
      return
   endif

   ; Read the input file
   flux = mrdfits(filename, 0, hdr)
   plate = sxpar(hdr,'PLATEID')
   mjd = sxpar(hdr,'MJD')
   if (keyword_set(flux) EQ 0) then $
    message, 'Unable to read file '+filename
   invvar = mrdfits(filename, 1)
   wave = mrdfits(filename, 2)
   dims = size(flux,/dimens)
   npix = dims[0]
   nfiber = dims[1]

   ; Read other HDUs of the spPlate file from existing BOSS reductions
   readspec, plate, mjd=mjd, andmask=andmask, ormask=ormask, disp=dismap, $
    plugmap=plugmap, sky=skyimg

   ; Construt the output file name
   if (keyword_set(platefile1)) then platefile = platefile1 $
    else platefile = 'spPlate-'+plate_to_string(plate)+'-' $
     +string(mjd,format='(i5.5)')+'.fits'
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile = 'spPPlate-'+plate_to_string(plate)+'-' $
     +string(mjd,format='(i5.5)')+'.fits'

   ; Rebin to the SDSS spacing
   dloglam = 1d-4
   newloglam = wavevector(alog10(min(wave)), alog10(max(wave)), $
    binsz=dloglam)
   newwave = 10^newloglam
   nnew = n_elements(newloglam)

   ; Simple linear interpolation of bad pixels
   ; Mask all pixels with INVVAR=0 or their neighbors
   gmask = invvar NE 0 ; =1 for good
   newbadmask = fltarr(nnew,nfiber)
   for i=0, nfiber-1 do $
    newbadmask[*,i] = rebin_spectrum(float(gmask[*,i] EQ 0), wave, newwave)
   newmask = newbadmask EQ 0 ; =1 for good

   ; Resample the flux and invvar using the COMBINE1FIBER logic
   newflux = fltarr(nnew,nfiber)
   newivar = fltarr(nnew,nfiber)
   loglam = alog10(wave)
   for i=0, nfiber-1 do begin
      combine1fiber, loglam, flux[*,i], invvar[*,i], $
       newloglam=newloglam, newflux=newflux1, newivar=newivar1
      newflux[*,i] = newflux1
      newivar[*,i] = newivar1
   endfor
   newflux = float(newflux) * newmask
   newivar = float(newivar) * newmask

   sxaddpar, hdr, 'PLATEID', plate
   sxaddpar, hdr, 'MJD', mjd
   sxaddpar, hdr, 'COEFF0', newloglam[0] ; wavelength solution
   sxaddpar, hdr, 'COEFF1', dloglam ; wavelength solution
   mwrfits, newflux, outfile, hdr, /create
   mwrfits, newivar, outfile
   mwrfits, andmask, outfile
   mwrfits, ormask, outfile
   mwrfits, dispmap, outfile
   mwrfits, plugmap, outfile
   mwrfits, skyimg, outfile

   return
end
;------------------------------------------------------------------------------
