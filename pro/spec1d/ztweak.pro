;+
; NAME:
;   ztweak
;
; PURPOSE:
;   Tweak redshifts
;
; CALLING SEQUENCE:
;   ztweak, platefile
;
; INPUTS:
;   platefile  - Plate file from spectro-2D
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/TEMPLATEFILES
;
; PROCEDURES CALLED:
;   spec_append
;   struct_addtags()
;   sxpar()
;   veldisp
;
; REVISION HISTORY:
;   19-Jul-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro ztweak, platefile

   eigenfile = 'spEigenVstnd.fits'
   if (n_elements(eigendir) EQ 0) then $
    eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')

platefile = 'spPlate-0306-51690.fits'
djs_readcol, '/home/schlegel/idlspec2d/etc/regress1d_all.dat', $
 chicplate, junk, chicfiberid, chicz, chicclass, format='(L,L,L,F,A)'
ii=where(chicplate EQ 306)
chicz=chicz[ii]
chicclass=chicclass[ii]

   zstruct1 = create_struct( $
    'class' , '', $
    'z'     , 0.0, $
    'chi2'  , 0.0, $
    'dof'   , 0.0 )
   zstruct = replicate(zstruct1, 640)
   zstruct.class = chicclass
   zstruct.z = chicz

   ;----------
   ; Read the 2D output file

   objflux = mrdfits(platefile,0,hdr)
   npix = sxpar(hdr, 'NAXIS1')
   nobj = sxpar(hdr, 'NAXIS2')
   objivar = mrdfits(platefile,1)
   andmask = mrdfits(platefile,2)
   ormask = mrdfits(platefile,3)

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-sub residuals.

   objivar = skymask(objivar, andmask, ormask)
andmask = 0 ; Free memory
ormask = 0 ; Free memory

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')
   objloglam = objloglam0 + lindgen(npix) * objdloglam

   ;----------
   ; Read the velocity standard template file.
   ; (Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.)

   starflux = readfits(djs_filepath(eigenfile, root_dir=eigendir), shdr)
   starloglam0 = sxpar(shdr, 'COEFF0')
   stardloglam = sxpar(shdr, 'COEFF1')
   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   starloglam = starloglam0 + stardloglam * dindgen(dims[0])
   if (ndim EQ 1) then nstar = 1 $
    else nstar = dims[1]

   ;---------------------------------------------------------------------------
   ; LOOP OVER EACH OBJECT
   ;---------------------------------------------------------------------------

   for iobj=0, nobj-1 do begin
print,'OBJECT ', iobj

      ;----------
      ; Construct a list of test redshifts near the initial redshift guess

      zinit = zstruct[iobj].z
ntest = 20
      ztest = zinit + (lindgen(ntest) - ntest/2) * 20. / 3.e5 ; space at 20 km/s
      rchi2arr = fltarr(ntest)

      ;----------
      ; Loop over each possible redshift, and redshift the spectrum of
      ; the velocity standard.
      ; Create a mask where the star does not overlap the object.

      for itest=0, ntest-1 do begin

         newstarflux = fltarr(npix,nstar)
         for istar=0, nstar-1 do begin
            combine1fiber, starloglam + alog10(1 + ztest[itest]), $
             starflux[*,istar], $
             newloglam=objloglam, newflux=tmpflux, maxiter=0
            newstarflux[*,istar] = tmpflux
            if (istar EQ 0) then newstarmask = tmpflux NE 0
         endfor

         ;----------
         ; Construct template as the star flux + emission lines + poly terms

npoly = 3

         ; Convert from air wavelengths in Ang to log-10 wavelength in vacuum,
         ; redshifted to the test redshift.

         linelist = [4861.3632, 4958.911, 5006.843, 6548.05, 6562.801, $
          6583.45, 6716.44, 6730.82]
         vaclist = linelist
         airtovac, vaclist
         vaclist = alog10(vaclist * (1 + ztest[itest]))
linesig = 2.e-4 ; Set the width to sigma = 2pix = 140 km/s (FWHM = 323 km/s)
         for iline=0, n_elements(vaclist)-1 do $
          newstarflux = [ [newstarflux], $
           [gaussian(objloglam,[1.,vaclist[iline],linesig])] ]

         if (keyword_set(npoly)) then $
          newstarflux = [ [newstarflux], [poly_array(npix,npoly)] ]

         ;----------
         ; Compute the chi^2

         zans = zcompute(objflux[*,iobj], objivar[*,iobj], $
          newstarflux, newstarmask, nfind=1, pmin=0, pmax=0)
         rchi2arr[itest] = zans.chi2 / (zans.dof + (zans.dof EQ 0))

      endfor

      ;----------
      ; Determine the best reduced chi^2

      bestval = min(rchi2arr, ibest)
      if (bestval NE 0) then zstruct[iobj].z = ztest[ibest]

   endfor
stop

   return
end
;------------------------------------------------------------------------------
