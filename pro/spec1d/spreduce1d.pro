;+
; NAME:
;   spreduce1d
;
; PURPOSE:
;   1-D reduction of spectra from 1 plate
;
; CALLING SEQUENCE:
;   spreduce1d, platefile
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
;   Presently, just writes a file 'z-'+PLATEFILE.
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
; INTERNAL SUPPORT ROUTINES:
;   sp1d_struct()
;
; REVISION HISTORY:
;   28-Jun-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function sp1d_struct, nobj

   result = { $
    plate               : 0L, $
    mjd                 : 0L, $
    fiberid             : 0L, $
    primtarget          : 0L, $
    sectarget           : 0L $
   }

   return, replicate(result, nobj)
end

;------------------------------------------------------------------------------
pro spreduce1d, platefile

outfile = 'z-' + platefile


djs_readcol, '/home/schlegel/idlspec2d/etc/regress1d_all.dat', $
 chicplate, junk, chicfiberid, chicz, chicclass, format='(L,L,L,F,A)'
ii=where(chicplate EQ 306)
chicz=chicz[ii]
chicclass=chicclass[ii]

   ;----------
   ; Read the 2D output file

   objflux = mrdfits(platefile,0,hdr)
   objivar = mrdfits(platefile,1)
   andmask = mrdfits(platefile,2)
;   ormask = mrdfits(platefile,3)
;   dispmap = mrdfits(platefile,4)
   plugmap = mrdfits(platefile,5)
plugmap = mrdfits(platefile,4)

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-subtraction
   ; residuals.

   ibad = where(andmask AND pixelmask_bits('BRIGHTSKY'))
andmask = 0 ; Free memory
   if (ibad[0] NE -1) then objivar[ibad] = 0

; Trim to QSO's only...
;jj = where(chicclass EQ 'QSO')
jj = where(chicclass EQ 'GALAXY')
;jj = jj[0:4] ; Trim to first few objects
chicz = chicz[jj]
chicclass = chicclass[jj]
objflux = objflux[*,jj]
objivar = objivar[*,jj]
plugmap = plugmap[jj]

   ;----------
   ; Determine the wavelength mapping for the object spectra,
   ; which are the same for all of them.

   npix = sxpar(hdr, 'NAXIS1')
   nobj = sxpar(hdr, 'NAXIS2')
nobj = n_elements(jj)
   objloglam0 = sxpar(hdr, 'COEFF0')
   objdloglam = sxpar(hdr, 'COEFF1')
   objwave = 10^(objloglam0 + objdloglam * findgen(npix))

   ;----------
   ; Look for where the S/N is unreasonably large

   for iobj=0, nobj-1 do begin
      junk = where(objflux[*,iobj] * sqrt(objivar[*,iobj]) GT 200., ct)
      if (ct GT 0) then $
       splog, 'WARNING: Fiber #', iobj+1, ' has ', ct, ' pixels with S/N > 200'
   endfor

   ;----------
   ; Read the template files
   ; Assume that the wavelength binning is the same as for the objects
   ; in log-wavelength.

templatefile = 'spTemplateQSO.dat' ; ???
zmin = 0
zmax = 6000
templatefile = ['spTemplateE.dat', 'spTemplateS.dat'] ; ???
zmin = 0
zmax = 1000
   ntemplate = n_elements(templatefile)
   for it=0, ntemplate-1 do begin
      tfile = filepath(templatefile[it], $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
      djs_readcol, tfile, twave, tflux, tivar, format='(D,F,F)'
      tloglam0 = alog10(twave[0])
      tdloglam = alog10(twave[1] / twave[0])

      ; Test if the binning is the same as for the objects.
      ; If not, then re-bin to the same pixel bin size.
      if (abs((tdloglam-objdloglam)/objdloglam) GT 0.01) then begin
         ii = where(tivar GT 0)
         tloglam0 = alog10(twave[ii[0]])
         nnewbin = alog10(max(twave[ii]) / min(twave[ii])) / objdloglam
         newloglam = tloglam0 + lindgen(nnewbin) * objdloglam
         tflux = interpol(tflux[ii], alog10(twave[ii]), newloglam)
         tivar = interpol(tivar[ii], alog10(twave[ii]), newloglam)
      endif

      if (it EQ 0) then temloglam0 = tloglam0 $
       else temloglam0 = [temloglam0, tloglam0]

      spec_append, temflux, tflux
      spec_append, temivar, tivar
   endfor

temivar = temivar * 0 + 1.0

   ;----------
   ; If the number of spectral pixels is larger in either the objects
   ; or the templates, than increase the smaller's size to agree.

   nptemp = (size(temflux,/dimens))[0]
;   if (nptemp LT npix) then begin
;      temflux = [temflux, fltarr(npix-nptemp,ntemplate)]
;      temivar = [temivar, fltarr(npix-nptemp,ntemplate)]
;   endif else if (nptemp GT npix) then begin
;      objflux = [objflux, fltarr(nptemp-npix,nobj)]
;      objivar = [objivar, fltarr(nptemp-npix,nobj)]
;   endif

   ;----------
   ; Mask out the object spectra near the 5577 sky line

   linemask = objwave GE 5573 AND objwave LE 5586 ; =1 for bad pixels
   linemask = rebin(linemask, npix, nobj)

   objflux = djs_maskinterp(objflux, linemask, iaxis=0)
   objivar = objivar * (1 - linemask)

   ;----------
   ; Compute the redshift difference between the first pixel of the object
   ; spectra and each template

   zoffset = (objloglam0 - temloglam0) / objdloglam
zoffset = zoffset[0]

   ;----------
   ; Compute the redshifts and velocity dispersions

 xvec = findgen(nptemp)/nptemp
; temflux = [ [temflux], [xvec^0], [xvec^1], [xvec^2], [xvec^3] ]
 temflux = [ [temflux], [xvec^0], [xvec^1], [xvec^2] ]
t0 = systime(1)
   daveres = computez(objflux, objivar, temflux, temivar, zoffset=zoffset, $
    zmin=zmin, zmax=zmax, nfind=5)
print,'TIME FOR COMPUTEZ = ', systime(1)-t0
davez = 10.^(objdloglam * daveres[*,0].z) - 1.
allz = 10.^(objdloglam * daveres.z) - 1.
print,[transpose(chicz),transpose(davez)]

;   vres = veldisp(objflux, objivar, temflux, temivar, $
;    zoffset=zoffset)
vres1 = vres ; backup copy

   ; Convert redshift (and error) from pixels to the conventional dimensionless
   ; value.
   vres.z = 10.^(objdloglam * vres.z) - 1.
   vres.z_err = alog(10) * objdloglam * vres.z_err * (1 + vres.z)

   result = sp1d_struct(nobj)
   result.plate = sxpar(hdr, 'PLATEID')
   result.mjd = sxpar(hdr, 'MJD')
   result.fiberid = plugmap.fiberid
   result.primtarget = plugmap.primtarget
   result.sectarget = plugmap.sectarget
   result = struct_addtags(result, vres)

   ;----------
   ; Write the output file

stop

set_plot,'ps'
device,file='qsoz-306a.ps'
djs_plot, chicz, davez,ps=4,xr=[0,3],yr=[0,3], charsize=2, $
 xtitle='Chicago-z',ytitle='z from \chi^2-min',title='QSOs on plate 306'
oplot,[0,3],[0,3]
device,/close
set_plot,'x'

set_plot,'ps'
device,file='qsoz-306b.ps'
djs_plot, chicz, ((davez-chicz)*3e5 > (-5800))<5800, $
 ps=4,yr=[-1,1]*6000,xr=[0,3], charsize=2, $
 xtitle='Chicago-z',ytitle='\Delta z from \chi^2-min',title='QSOs on plate 306'
oplot,[0,3],[0,0]
device,/close
set_plot,'x'

djs_plot, chicz, ((davez-chicz)*3e5 > (-580))<580, $
 ps=4,yr=[-1,1]*600,xr=[0,1], charsize=2
j=where( abs(davez-chicz)*3e5 LT 400 )

j=where( abs(davez-chicz)*3e5 LT 4000 )
print,mean((davez-chicz)[j]*3e5),djsig((davez-chicz)[j]*3e5)
print,n_elements(j)/float(nobj)

   mwrfits, result, outfile, /create

   return
end
;------------------------------------------------------------------------------
