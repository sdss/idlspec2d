;+
; NAME:
;   gama_to_sdss
;
; PURPOSE:
;   Convert a 2dF GAMA extracted spectra file to an SDSS spPlate format
;
; CALLING SEQUENCE:
;   gama_to_sdss, filename, [ outfile= ]
;
; INPUTS:
;   filename   - 2dF GAMA spectra file, where HDU #0 is the flux [NPIX,NFIBER],
;                HDU #1 is the variance, and HDU #2 contains catalog info
;
; OPTIONAL INPUTS:
;   outfile    - Output file name; default to spPlate-$PLATE-$MJD.fits where
;                PLATE is derived from the GAMA field name and tile name,
;                and MJD from the data of observation.  GAMA field names
;                02, 09, 12, 15, 23 are mapped to plate numbers
;                1000, 2000, 3000, 4000, 5000 with the GAMA pointing number
;                added to that.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The 2dF spectra are interpolated to the log-wavelength mapping
;   of SDSS spectra using a B-spline to the flux and the variance arrays.
;
; EXAMPLES:
;   Convert a 2dF GAMA file and run the SDSS redshifting code:
;     IDL> gama_to_sdss,'SG09_Y3_012.fits',outfile='spPlate-2012-55246.fits'
;     IDL> spreduce1d,'spPlate-2012-55246.fits'
;
; BUGS:
;   The 2dF wavelength scale is assumed to be vacuum barycentric.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   23-Feb-2011  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro gama_to_sdss, filename, outfile=outfile1

   if (size(filename,/tname) NE 'STRING') then begin
      splog, 'Must specify FILENAME'
      return
   endif

   plate = long(strmid(fileandpath(filename),8,3))
   fieldname = strmid(fileandpath(filename),2,2)
   case fieldname of
      '02' : plate += 1000L
      '09' : plate += 2000L
      '12' : plate += 3000L
      '15' : plate += 4000L
      '23' : plate += 5000L
      else : plate += 9000L
   endcase

   dradeg = 180d0 / !dpi

   plug1 = create_struct( $
    'OBJID', lonarr(5), $
    'RA', 0d0, $
    'DEC', 0d0, $
    'OBJTYPE', '', $
    'SPECTROGRAPHID', 1L, $
    'FIBERID', 0L )

   flux = mrdfits(filename, 0, hdr)
   if (keyword_set(flux) EQ 0) then $
    message, 'Unable to read file '+filename
   var = mrdfits(filename, 1)
   cat = mrdfits(filename, 2)
   dims = size(flux,/dimens)
   npix = dims[0]
   nfiber = dims[1]
   dwave = sxpar(hdr,'CDELT1')
   wave = (dindgen(npix) - sxpar(hdr,'CRPIX1') - 1) * dwave $
    + sxpar(hdr,'CRVAL1')

   ; Construt the output file name
   mjd = sxpar(hdr,'UTMJD')
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile = 'spPlate-'+string(plate,format='(i4.4)')+'-' $
     +string(mjd,format='(i5.5)')+'.fits'

   ; Rebin to the SDSS spacing
   dloglam = 1d-4
   newloglam = wavevector(alog10(min(wave)), alog10(max(wave)), binsz=dloglam)
   newwave = 10^newloglam
   nnew = n_elements(newloglam)
   gmask = finite(flux) AND finite(var) AND var GT 0
   ibad = where(gmask EQ 0, nbad)
   if (nbad GT 0) then flux[ibad] = 0
   if (nbad GT 0) then var[ibad] = 0


; Simple linear interpolation
   newbadmask = fltarr(nnew,nfiber)
   for i=0, nfiber-1 do begin
;      newflux[*,i] = rebin_spectrum(flux[*,i], wave, newwave)
;      newvar[*,i] = rebin_spectrum(var[*,i], wave, newwave)
      newbadmask[*,i] = rebin_spectrum(float(gmask[*,i] EQ 0), wave, newwave)
   endfor
   newmask = newbadmask EQ 0
;   newivar = newmask / (newvar + (newmask EQ 0))

   newflux = fltarr(nnew,nfiber)
   newivar = fltarr(nnew,nfiber)
   loglam = alog10(wave)
   for i=0, nfiber-1 do begin
      combine1fiber, loglam, flux[*,i], 1./(var[*,i] + (var[*,i] EQ 0)), $
       newloglam=newloglam, newflux=newflux1, newivar=newivar1
      newflux[*,i] = newflux1
      newivar[*,i] = newivar1
   endfor
   newflux = float(newflux) * newmask
   newivar = float(newivar) * newmask

   ; Get the units right, making them per Ang although the spacing
   ; is no longer uniform in wavelength
   wratio = rebin((10d0^(dloglam) - 1d0) * newwave / dwave, [nnew, nfiber])
   newflux /= wratio
   newivar *= wratio^2

   andmask = lonarr(npix,nfiber)
   ormask = lonarr(npix,nfiber)
   dispmap = fltarr(npix,nfiber) + 1.
   skyimg = fltarr(npix,nfiber)
   plugmap = replicate(plug1, nfiber)
   plugmap.ra = cat.ra * dradeg
   plugmap.dec = cat.dec * dradeg
   objtype = strmid(cat.name,0,1)
   plugmap.objtype = (objtype EQ 'G' ? 'GALAXY' : '') $
    + (objtype EQ 'X' ? 'SKY' : '') $
    + (objtype EQ 'S' ? 'SPECTROPHOTO_STD' : '') $
    + (objtype EQ 'P' ? 'NA' : '') $
    + (objtype EQ 'F' ? 'NA' : '') ; fiducial/guide
   plugmap.objtype = 'GALAXY' ; default value
   plugmap.fiberid = lindgen(nfiber) + 1

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
