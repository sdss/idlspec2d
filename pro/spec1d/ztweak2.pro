;+
; NAME:
;   ztweak2
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
pro ztweak2, platefile

platefile = 'spPlate-0321-51612-60sec.fits'

   zstruct1 = create_struct( $
    'class' , '', $
    'z'     , 0.0, $
    'chi2'  , 0.0, $
    'dof'   , 0.0 )
   zstruct = replicate(zstruct1, 640)

   ;----------
   ; Read the 2D output file

   objflux = mrdfits(platefile,0,hdr)
   npix = sxpar(hdr, 'NAXIS1')
   nobj = sxpar(hdr, 'NAXIS2')
   objivar = mrdfits(platefile,1)
   andmask = mrdfits(platefile,2)
   ormask = mrdfits(platefile,3)
   plug = mrdfits(platefile,5)

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-sub residuals.

   objivar = skymask(objivar, andmask, ormask)
andmask = 0 ; Free memory
ormask = 0 ; Free memory

zstruct.z = plug.expl / 3.e5
j = where(plug.expl NE -999 AND plug.expl NE 0)
zstruct[j].class = 'STAR'

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')

   ;----------
   ; Read the template files
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

   starflux = objflux[*,93]
   starivar = objivar[*,93]
   stardloglam = objdloglam
zstar = 33.9 / 3.e5  ; redshift of this star ???
   starloglam0 = objloglam0 - alog10(1 + zstar) ; ???

   ;----------
   ; Construct template as the star flux + emission lines + polynomial terms

npstar = n_elements(starflux)
npoly = 4

;   linelist = [4861.3632, 4958.911, 5006.843, 6548.05, 6562.801, $
;    6583.45, 6716.44, 6730.82]
;   ; Convert from air wavelengths in Angstroms to log-10 wavelength in vacuum
;   vaclist = linelist
;   airtovac, vaclist
;   vaclist = alog10(vaclist)
;linesig = 2.e-4 ; Set the width to sigma = 2pix = 140 km/s (FWHM = 323 km/s)
;   for iline=0, n_elements(vaclist)-1 do $
;    starflux = [ [starflux], $
;     [gaussian(starloglam,[1.,vaclist[iline],linesig])] ]

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
      if (zstruct[iobj].class EQ 'STAR') then begin
         zpixobj = alog10(1 + zstruct[iobj].z) / objdloglam
         znew = computez(objflux[*,iobj], objivar[*,iobj], $
          starflux, starivar NE 0, zoffset=zoffset, $
          zmin=long(zpixobj-5), zmax=long(zpixobj+5))
         zstruct[iobj].z = 10.^(znew.z * objdloglam) - 1.
         zstruct[iobj].chi2 = znew.chi2
         zstruct[iobj].dof = znew.dof
      endif
   endfor
stop
j=where(zstruct.class EQ 'STAR')
splot,plug[j].expl,3e5*zstruct[j].z-plug[j].expl,ps=3,yr=[-40,40]
j1=where(zstruct.class EQ 'STAR' AND plug.fiberid LE 320)
j2=where(zstruct.class EQ 'STAR' AND plug.fiberid GT 320)
soplot,plug[j1].expl,3e5*zstruct[j1].z-plug[j1].expl,color='red',ps=3

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
