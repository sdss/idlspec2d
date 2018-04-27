;+
; NAME:
;   desi_to_sdss
;
; PURPOSE:
;   Convert a DESI extracted spectra file to an SDSS spPlate format
;
; CALLING SEQUENCE:
;   desi_to_sdss, filename, plate, mjd, [ outfile= ]
;
; INPUTS:
;   filename   - DESI spectra file
;   plate      - Plate label, must be valid for sdss_specobjid()
;   mjd        - MJD label, must be valid for sdss_specobjid()
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
;   The DESI spectra are interpolated to the log-wavelength mapping
;   of SDSS spectra using a B-spline to the flux and the variance arrays.
;
;   The DESI input file assumed to be of this format:
;     http://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/spectra-NSIDE/PIXGROUP/PIXNUM/spectra-NSIDE-PIXNUM.html#hdu01
;
; EXAMPLES:
;   Convert a DESI file and run the SDSS redshifting code:
;     IDL> desi_to_sdss,'weird_obj89.darksky.fits',1,56000
;     IDL> spreduce1d,'spPlate-0001-56000.fits'
;     IDL> zans=mrdfits('spZbest-0001-56000.fits',1)
;     IDL> print,zans.z,zans.z_err
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Apr-2018  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro desi_to_sdss, filename, plate, mjd, outfile=outfile1

   if (size(filename,/tname) NE 'STRING') then begin
      splog, 'Must specify FILENAME'
      return
   endif
   if (NOT keyword_set(plate)) then $
    message, 'Must set PLATE'
   if (NOT keyword_set(mjd)) then $
    message, 'Must set MJD'

   dradeg = 180d0 / !dpi

   plug1 = create_struct( $
    'OBJID', lonarr(5), $
    'RA', 0d0, $
    'DEC', 0d0, $
    'OBJTYPE', '', $
    'SPECTROGRAPHID', 1L, $
    'FIBERID', 0L )

   fibermap = mrdfits(filename, 'FIBERMAP')
   b_wave = mrdfits(filename, 'B_WAVELENGTH')
   b_flux = mrdfits(filename, 'B_FLUX')
   b_ivar = mrdfits(filename, 'B_IVAR')
   b_mask = mrdfits(filename, 'B_MASK')
   b_res = mrdfits(filename, 'B_RESOLUTION')
   r_wave = mrdfits(filename, 'R_WAVELENGTH')
   r_flux = mrdfits(filename, 'R_FLUX')
   r_ivar = mrdfits(filename, 'R_IVAR')
   r_mask = mrdfits(filename, 'R_MASK')
   r_res = mrdfits(filename, 'R_RESOLUTION')
   z_wave = mrdfits(filename, 'Z_WAVELENGTH')
   z_flux = mrdfits(filename, 'Z_FLUX')
   z_ivar = mrdfits(filename, 'Z_IVAR')
   z_mask = mrdfits(filename, 'Z_MASK')
   z_res = mrdfits(filename, 'Z_RESOLUTION')
   if (keyword_set(fibermap) EQ 0) then $
    message, 'Unable to read file '+filename

   dims = size(b_flux,/dimens)
   if (size(b_flux,/n_dimen) EQ 1) then nfiber = 1 $
    else nfiber = dims[1]

   ; Construt the output file name
   if (keyword_set(outfile1)) then outfile = outfile1 $
    else outfile = 'spPlate-'+plate_to_string(plate)+'-' $
     +string(mjd,format='(i5.5)')+'.fits'

   ; Rebin to the SDSS spacing
   dloglam = 1d-4
   newloglam = wavevector(alog10(min(b_wave)), alog10(max(z_wave)), binsz=dloglam)
   newwave = 10^newloglam
   nnew = n_elements(newloglam)

   newflux = fltarr(nnew,nfiber)
   newivar = fltarr(nnew,nfiber)

   b_npix = (size(b_flux,/dimens))[0]
   r_npix = (size(r_flux,/dimens))[0]
   z_npix = (size(z_flux,/dimens))[0]
   npix_all = (b_npix > r_npix) > z_npix
   for i=0, nfiber-1 do begin
      wave_all = dblarr(npix_all,3)
      flux_all = fltarr(npix_all,3)
      ivar_all = fltarr(npix_all,3)
      wave_all[0:b_npix-1,0] = b_wave
      if (b_npix LT npix_all) then $
       wave_all[b_npix:npix_all-1,0] = b_wave[b_npix-1] + dindgen(npix_all-b_npix) $
        * (b_wave[b_npix-1] - b_wave[b_npix-2])
      flux_all[0:b_npix-1,0] = b_flux[*,i]
      ivar_all[0:b_npix-1,0] = b_ivar[*,i]
      wave_all[0:r_npix-1,1] = r_wave
      if (r_npix LT npix_all) then $
       wave_all[r_npix:npix_all-1,1] = r_wave[r_npix-1] + dindgen(npix_all-r_npix) $
        * (r_wave[r_npix-1] - r_wave[r_npix-2])
      flux_all[0:r_npix-1,1] = r_flux[*,i]
      ivar_all[0:r_npix-1,1] = r_ivar[*,i]
      wave_all[0:z_npix-1,2] = z_wave
      if (z_npix LT npix_all) then $
       wave_all[z_npix:npix_all-1,2] = z_wave[z_npix-1] + dindgen(npix_all-z_npix) $
        * (z_wave[z_npix-1] - z_wave[z_npix-2])
      flux_all[0:z_npix-1,2] = z_flux[*,i]
      ivar_all[0:z_npix-1,2] = z_ivar[*,i]
      combine1fiber, alog10(wave_all), flux_all, ivar_all, $
       newloglam=newloglam, newflux=newflux1, newivar=newivar1
      newflux[*,i] = newflux1
      newivar[*,i] = newivar1
print,i,nfiber,string(13b),format='(i4,i4,a,$)'
   endfor
print

   andmask = lonarr(nnew,nfiber)
   ormask = lonarr(nnew,nfiber)
   dispmap = fltarr(nnew,nfiber) + 1.
   skyimg = fltarr(nnew,nfiber)
   plugmap = replicate(plug1, nfiber)
   plugmap.ra = fibermap.ra_obs
   plugmap.dec = fibermap.dec_obs
   plugmap.objtype = fibermap.objtype
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
