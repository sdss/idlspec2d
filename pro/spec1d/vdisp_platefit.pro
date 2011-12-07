;+
; NAME:
;   vdisp_platefit
;
; PURPOSE:
;   re-run redshift-marginalized vdispfit over a plate and write out results.
;
; CALLING SEQUENCE:
;   vdisp_platefit, plate, mjd
;
; READS:
;   spPlate and spZbest files
;
; WRITES:
;   spZvfit file, which is just spZbest HDU1 with vdisp_lnl changed.
;
; WRITTEN:
;   A. Bolton, Utah, Dec 2011 (chopped out & modified from spreduce1d)
;
;------------------------------------------------------------------------------
pro vdisp_platefit, plate, mjd

pfile = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') + $
 '/' + string(plate, format='(i4.4)') + '/spPlate-' + string(plate, format='(i4.4)') + $
 '-' + string(mjd, format='(i5.5)') + '.fits'

if (not file_test(pfile)) then begin
   splog, 'spPlate file not found:'
   splog, pfile
   return
endif

zfile = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') + $
 '/' + string(plate, format='(i4.4)') + '/' + getenv('RUN1D') + $
 '/spZbest-' + string(plate, format='(i4.4)') + $
 '-' + string(mjd, format='(i5.5)') + '.fits'

ozfile = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') + $
 '/' + string(plate, format='(i4.4)') + '/' + getenv('RUN1D') + $
 '/spZvfit-' + string(plate, format='(i4.4)') + $
 '-' + string(mjd, format='(i5.5)') + '.fits'

if (not file_test(zfile)) then begin
   splog, 'spZbest file not found:'
   splog, zfile
   return
endif

objflux = mrdfits(pfile,0,hdr)
if (NOT keyword_set(hdr)) then $
   message, 'Plate file not valid: ' + pfile
npixobj = sxpar(hdr, 'NAXIS1')
nobj = sxpar(hdr, 'NAXIS2')
objivar = mrdfits(pfile,1)
andmask = mrdfits(pfile,2)
ormask = mrdfits(pfile,3)
plugmap = mrdfits(pfile,5)
skyflux = mrdfits(pfile,6)
objloglam0 = sxpar(hdr, 'COEFF0')                                       
objdloglam = sxpar(hdr, 'COEFF1')     
zhdr = headfits(zfile)
zans = mrdfits(zfile,1)
objloglam = objloglam0 + lindgen(npixobj) * objdloglam

anyandmask = transpose(andmask[0,*])
anyormask = transpose(ormask[0,*])
for ipix=1, npixobj-1 do $
   anyandmask = anyandmask OR transpose(andmask[ipix,*])
for ipix=1, npixobj-1 do $
   anyormask = anyormask OR transpose(ormask[ipix,*])
objivar = skymask(objivar, andmask, ormask)
andmask = 0                     ; Free memory
ormask = 0                      ; Free memory

;----------
; Mask out points that are unphysically negative (10-sigma negatives),
; and mask the neighboring 2 pixels in each direction.
for iobj=0L, nobj-1 do begin
   thismask = objflux[*,iobj] * sqrt(objivar[*,iobj]) LE -10.
   thismask = smooth(float(thismask),5) GT 0
   objivar[*,iobj] = objivar[*,iobj] * (1 - thismask)
endfor

;-----------
; Compute redshift-marginalized velocity-dispersion
; likelihood curves for galaxies (bolton@utah 2011aug):
wh_v = where(strmatch(zans.objtype, 'GALAXY*') and $
             strmatch(zans.class_noqso, 'GALAXY*'), n_v)
if (n_v gt 0) then begin
   vdans_new = vdispfit(objflux[*,wh_v], objivar[*,wh_v], objloglam, $
                        zobj=zans[wh_v].z_noqso, eigenfile='spEigenElodie.fits', $
                        columns=lindgen(5), npoly=5, dzpix=3, /return_chisq)
   vdans_new = reform(vdans_new)
   nchi2 = n_elements(vdans_new[0].chi2arr)
   vstruc_new = replicate({vdisp_lnl: replicate(0., nchi2)}, nobj)
   vstruc_new[wh_v].vdisp_lnl = - vdans_new.chi2arr / 2.
   newzans = struct_selecttags(zans, except_tags='vdisp_lnl')
   newzans = struct_addtags(newzans, vstruc_new)
endif

;----------
; Write the output file:
mwrfits, 0, ozfile, zhdr, /create ; Retain the original spZbest header
mwrfits, newzans, ozfile

   return
end
;------------------------------------------------------------------------------
