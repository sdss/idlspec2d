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

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-subtraction
   ; residuals.  Grow that mask by 2 pixels in each direction.

   skymask = (andmask AND pixelmask_bits('BRIGHTSKY')) NE 0
   for iobj=0, nobj-1 do $
    skymask[*,iobj] = smooth(float(skymask[*,iobj]),5) GT 0
   ibad = where(skymask)
andmask = 0 ; Free memory
   if (ibad[0] NE -1) then objivar[ibad] = 0

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')

   ;----------
   ; Read the template files
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

   readspec, 306, 250, flux=starflux, invvar=starivar, loglam=starloglam
   starloglam0 = starloglam[0]
zstar = 0.0  ; redshift of this star - change it to zero ???

   ;----------
   ; Construct template as the star flux + emission lines + polynomial terms

npstar = n_elements(starflux)
npoly = 3

   linelist = [4861.3632, 4958.911, 5006.843, 6548.05, 6562.801, $
    6583.45, 6716.44, 6730.82]
   ; Convert from air wavelengths in Angstroms to log-10 wavelength in vacuum
   vaclist = linelist
   airtovac, vaclist
   vaclist = alog10(vaclist)
linesig = 2.e-4 ; Set the width to sigma = 2pix = 140 km/s (FWHM = 323 km/s)
   for iline=0, n_elements(vaclist)-1 do $
    starflux = [ [starflux], $
     [gaussian(starloglam,[1.,vaclist[iline],linesig])] ]

   if (keyword_set(npoly)) then $
    starflux = [ [starflux], [poly_array(npstar,npoly)] ]

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and each template

   zoffset = (objloglam0 - starloglam0) / objdloglam

   ;----------
   ; Compute the redshifts

   for iobj=0, nobj-1 do begin
print, 'Object #', iobj
      if (zstruct[iobj].class EQ 'GALAXY') then begin
         zpixobj = alog10(1 + zstruct[iobj].z) / objdloglam
         znew = computez(objflux[*,iobj], objivar[*,iobj], $
          starflux, starivar, zoffset=zoffset, $
;          starflux[*,0], starivar, zoffset=zoffset, $
;          zmin=long(zpixobj-400), zmax=long(zpixobj+400))
          zmin=long(zpixobj-20), zmax=long(zpixobj+20))
         zstruct[iobj].z = 10.^(znew.z * objdloglam) - 1.
         zstruct[iobj].chi2 = znew.chi2
         zstruct[iobj].dof = znew.dof
      endif
if (iobj EQ 20) then stop
   endfor
stop
j=where(zstruct.class EQ 'GALAXY')
splot,chicz[j],(zstruct[j].z-chicz[j])*3e5,ps=2

zdiff=(zstruct.z-chicz)*3e5
jj=where(zstruct.class EQ 'GALAXY' AND $
 zdiff LT 0 AND zdiff GT -400 AND chicz LT 0.5, ct)

get_lun, olun
openw, olun, 'eigeninputs.dat'
printf, olun, 'Plate  MJD    Fib   z        '
printf, olun, '-----  -----  ----  ---------'
for i=0, ct-1 do $
 printf, olun, 306, 51690, jj[i]+1, zstruct[jj[i]].z, $
  format='(i5,i7,i5,f12.6)'
close, olun
free_lun, olun

   return
end
;------------------------------------------------------------------------------
