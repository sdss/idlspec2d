;+
; NAME:
;   ztweak_gal
;
; PURPOSE:
;   Tweak galaxy redshifts
;
; CALLING SEQUENCE:
;   ztweak_gal, platefile, [subsamp= ]
;
; INPUTS:
;   platefile  - Plate file from spectro-2D
;
; OPTIONAL INPUTS:
;   subsamp    - ???
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
;  Make sure we use exactly the same object pixels for eadh test z ???
;  Convovle the emission lines just like the stars - incl. instrumental
;    response???
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
pro ztweak_gal, platefile, subsamp=subsamp

   eigenfile = 'spEigenVstnd*.fits'
   if (n_elements(eigendir) EQ 0) then $
    eigendir = concat_dir(getenv('IDLSPEC2D_DIR'), 'templates')

   if (NOT keyword_set(subsamp)) then subsamp = 10

platefile = 'spPlate-0306-51690.fits'
djs_readcol, getenv('IDLSPEC2D_DIR') + '/etc/regress1d_all.dat', $
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
   ; Find the most recent template file matching EIGENFILE

   allfiles = findfile(djs_filepath(eigenfile, root_dir=eigendir), count=ct)
   if (ct EQ 0) then $
    message, 'Unable to find EIGENFILE matching '+eigenfile
   thisfile = allfiles[ (reverse(sort(allfiles)))[0] ]
   splog, 'Selecting EIGENFILE=', thisfile

   ;----------
   ; Read the velocity standard template file.
   ; (Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.)

   starflux = readfits(thisfile, shdr)
   starloglam0 = sxpar(shdr, 'COEFF0')
   stardloglam = sxpar(shdr, 'COEFF1')
   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   npstar = dims[0]
   starloglam = starloglam0 + stardloglam * dindgen(npstar)
   if (ndim EQ 1) then nstar = 1 $
    else nstar = dims[1]

   ;---------------------------------------------------------------------------
   ; GENERATE OVER-SAMPLED EIGENSPECTRA
   ;---------------------------------------------------------------------------

   nbigpix = npstar * subsamp
   bigloglam = rebin(starloglam, nbigpix)

   bigstarflux = fltarr(nbigpix,nstar)
   for istar=0, nstar-1 do begin
      combine1fiber, starloglam, starflux[*,istar], $
       newloglam=bigloglam, newflux=tmpflux, maxiter=0
      bigstarflux[*,istar] = tmpflux
      if (istar EQ 0) then bigstarmask = tmpflux NE 0
   endfor

   ;----------
   ; Construct template as the star flux + emission lines + poly terms

npoly = 3

   ; Convert from air wavelengths in Ang to log-10 wavelength in vacuum.

   linelist = [4861.3632, 4958.911, 5006.843, 6548.05, 6562.801, $
    6583.45, 6716.44, 6730.82]
   vaclist = linelist
   airtovac, vaclist
   vaclist = alog10(vaclist)
linesig = 2.e-4 ; Set the width to sigma = 2pix = 140 km/s (FWHM = 323 km/s)
   for iline=0, n_elements(vaclist)-1 do $
    bigstarflux = [ [bigstarflux], $
     [gaussian(bigloglam,[1.,vaclist[iline],linesig])] ]

   if (keyword_set(npoly)) then $
    bigstarflux = [ [bigstarflux], [poly_array(nbigpix,npoly)] ]

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and the template.
   ; This is for the original sampling, not the over-sampled template.

   poffset = (objloglam0 - starloglam0) / objdloglam

   ;---------------------------------------------------------------------------
   ; LOOP OVER EACH OBJECT
   ;---------------------------------------------------------------------------

   for iobj=0, nobj-1 do begin
print,'OBJECT ', iobj

      ; PINIT = the galaxy pixel that should correspond to the first
      ;   template pixel; this is a fractional value

      zinit = zstruct[iobj].z
      pinit = - poffset + alog10(1. + zinit) / objdloglam

      ; P1VEC = the pixel number in the over-sampled template that corresponds
      ;   to the first galaxy pixel; this is an integer, but could be negative

ntest = 100
      p1vec = long(-pinit * subsamp + (lindgen(ntest) - ntest/2))
      zvec = 10.d^( (float(-p1vec)/subsamp + poffset) * objdloglam) - 1
; The above doesn't quite seem to be correct - always off by 1 pix???
zvec = dblarr(ntest)

      ;----------
      ; Loop over each possible redshift, and redshift the spectrum of
      ; the velocity standard.

      rchi2arr = fltarr(ntest)

      for itest=0, ntest-1 do begin

         ; P1 = First pixel in over-sampled template to use
         ; INDX = Indicies of the over-sampled template to use;
         ;   these will be separated by SUBSAMP
         ; Q1 = First pixel in object spectrum to use
         ; Q2 = Last pixel in object spectrum to use
         ; NPSMALL = Number of pixels for trimmed object and template;
         ;   if this is negative, then the two spectra do not overlap

         p1 = p1vec[itest]

         if (p1 GE 0) then begin
            q1 = 0L
         endif else begin
            q1 = (-p1+subsamp-1) / subsamp ; integer-valued
            p1 = p1 + subsamp + q1 * subsamp ; now positive-valued
         endelse

         npsmall = ((nbigpix - p1) / subsamp) < (npix - q1) ; integer-valued

         if (npsmall GT 1) then begin

            indx = p1 + subsamp * lindgen(npsmall)
zvec[itest] = 10d^(objloglam[q1] - bigloglam[indx[0]]) - 1 ; This is correct???

            ;----------
            ; Compute the chi^2

            q2 = q1 + npsmall - 1
            zans = zcompute(objflux[q1:q2,iobj], objivar[q1:q2,iobj], $
             bigstarflux[indx,*], nfind=1, pmin=0, pmax=0)
            rchi2arr[itest] = zans.chi2 / (zans.dof + (zans.dof EQ 0))
         endif

      endfor

      ;----------
      ; Determine the best reduced chi^2

      bestval = min(rchi2arr, ibest)
      if (bestval NE 0) then zstruct[iobj].z = zvec[ibest]
splot, zvec, rchi2arr

   endfor
stop

djs_readcol, getenv('IDLSPEC2D_DIR') + '/templates/eigeninput_gal.dat', $
 plate, mjd, fiberid, zold, format='(L,L,L,F)'
ii = fiberid-1
splot,zold,(zstruct[ii].z-zold)*3e5,ps=4,syms=0.4
openw, olun, 'eigeninput_gal.dat', /get_lun
for j=0, n_elements(ii)-1 do $
 printf, olun, plate[j], mjd[j], fiberid[j], zstruct[ii[j]].z, $
  format='(i5, i7, i5, f12.6)'
close, olun

vdiff = (zstruct[ii].z-zold)*3e5
vdiff = vdiff - median(vdiff)
ibad = where(abs(vdiff) GT 100)
print,fiberid[ibad]

   return
end
;------------------------------------------------------------------------------
