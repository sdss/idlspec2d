;------------------------------------------------------------------------------
; De-redshift and rebin spectra to specified wavelengths [WAVE1,WAVE2]
pro redindices_rebin, objloglam, objflux, objivar, zfit, $
 wave1, wave2, newloglam, newflux, newivar

   if (size(objflux,/n_dimen) EQ 1) then nobj = 1 $
    else nobj = (size(objflux,/dimens))[1]

   dloglam = 1.d-4 ; fix the wavelength spacing
   loglam1 = alog10(wave1) - 5 * dloglam ; pad with 5 pixels on left
   loglam2 = alog10(wave2) + 5 * dloglam ; pad with 5 pixels on right
   newloglam = loglam1 + dindgen((loglam2-loglam1) / dloglam) * dloglam
   nnew = n_elements(newloglam)
   newflux = fltarr(nnew,nobj)
   newivar = fltarr(nnew,nobj)

   ;----------
   ; Shift each spectra to z=0 and sample at the output wavelengths

   logshift = alog10(1.d + zfit)
   for iobj=0, nobj-1 do begin
      indx = where(objloglam[*,iobj] GT 0)
print, format='("Shift object ",i5," of ",i5,a1,$)', $
 iobj, nobj, string(13b)
      combine1fiber, objloglam[indx,iobj]-logshift[iobj], $
       objflux[indx,iobj], objivar[indx,iobj], $
       newloglam=newloglam, binsz=dloglam, newflux=flux1, newivar=ivar1
      newflux[*,iobj] = flux1
      newivar[*,iobj] = ivar1
   endfor
print,''

   return
end

;------------------------------------------------------------------------------
function redindices_mean, flux, ivar, mnerr=mnerr

   if (size(flux,/n_dimen) EQ 1) then nobj = 1 $
    else nobj = (size(flux,/dimens))[1]

   mn = fltarr(nobj)
   mnerr = fltarr(nobj)

   for iobj=0, nobj-1 do begin
      wtot = total(ivar[*,iobj])
      mn[iobj] = total(flux[*,iobj] * ivar[*,iobj]) / (wtot + (wtot EQ 0))
      mnerr[iobj] = 1. / sqrt(wtot + (wtot EQ 0))
   endfor

   return, mn
end

;------------------------------------------------------------------------------
function redindices_ew, objloglam, objflux, objivar, zfit, $
 lowave, midwave, hiwave

   redindices_rebin, objloglam, objflux, objivar, zfit, $
    lowave[0], hiwave[1], newloglam, newflux, newivar
   ii = where(newloglam GE alog10(lowave[0]) AND newloglam LE alog10(lowave[1]))
   mean1 = redindices_mean(newflux[ii,*], newivar[ii,*], mnerr=mnerr1)
   ii = where(newloglam GE alog10(midwave[0]) $
    AND newloglam LE alog10(midwave[1]))
   mean2 = redindices_mean(newflux[ii,*], newivar[ii,*], mnerr=mnerr2)
   ii = where(newloglam GE alog10(hiwave[0]) AND newloglam LE alog10(hiwave[1]))
   mean3 = redindices_mean(newflux[ii,*], newivar[ii,*], mnerr=mnerr3)
   continuum = (mean1 + mean3) / 2.
   ew = (mean2 - continuum) * (midwave[1] - midwave[0]) / continuum

   return, ew
end

;------------------------------------------------------------------------------
function redindices, plate, mjd=mjd

print, '--> WORKING ON PLATE ', plate
   maglimit = 18.0 ; de-reddened r-band magnitude limit

   ;----------
   ; Pre-select galaxies

   readspec, plate, mjd=mjd, zans=zans, plug=plug, tsobj=tsobj
   indx = where(strtrim(zans.class) EQ 'GALAXY' AND zans.zwarning EQ 0 $
    AND plug.mag[2] - tsobj.reddening[2,*] LT maglimit)

   ;----------
   ; Read the 2D output file

   readspec, plate, plug[indx].fiberid, mjd=mjd, $
    flux=objflux, invvar=objivar, loglam=objloglam, $
    andmask=andmask, ormask=ormask
   zans = zans[indx]
   tsobj = tsobj[indx]

   ;----------
   ; Mask around bright sky lines

   objivar = skymask(objivar, andmask, ormask)
andmask = 0 ; Free memory
ormask = 0 ; Free memory

   ;----------
   ; Measure the 4000-Ang break.

print, 'Measuring 4000-Ang break'
   lowave = [3800.,3900.]
   hiwave = [4025.,4200.]
   redindices_rebin, objloglam, objflux, objivar, zans.z, $
    lowave[0], hiwave[1], newloglam, newflux, newivar
   ii = where(newloglam GE alog10(lowave[0]) AND newloglam LE alog10(lowave[1]))
   mean1 = redindices_mean(newflux[ii,*], newivar[ii,*], mnerr=mnerr1)
   ii = where(newloglam GE alog10(hiwave[0]) AND newloglam LE alog10(hiwave[1]))
   mean2 = redindices_mean(newflux[ii,*], newivar[ii,*], mnerr=mnerr2)
   dinv4000 = mean1 / mean2

   ;----------
   ; Measure the Halpha+[NII] equiv. width

print, 'Measuring H-alpha EW'
   thiswave = 6564.
   lowave = thiswave + [-75,-35]
   midwave = thiswave + [-35,35]
   hiwave = thiswave + [35,75]
   ew_halpha = redindices_ew(objloglam, objflux, objivar, zans.z, $
    lowave, midwave, hiwave)

   ;----------
   ; Measure the Mg2 equiv. width

print, 'Measuring Mg2 EW'
   lowave = [4897,4958]
   midwave = [5156,5197]
   hiwave = [5303,5367]
   airtovac, lowave
   airtovac, midwave
   airtovac, hiwave
   ew_mg2 = redindices_ew(objloglam, objflux, objivar, zans.z, $
    lowave, midwave, hiwave)

   ;----------
   ; Now select our sample

   iselect = where(dinv4000 LT 0.6 AND ew_halpha LT 3.0)

   modelcolor = tsobj.counts_model[0:3,*] - tsobj.counts_model[1:4,*]

   redstruct = create_struct( $
    'modelcolor', fltarr(4), $
    'redsample', 0L, $
    'dinv4000', 0.0, $
    'ew_halpha', 0.0, $
    'ew_mg2', 0.0, $
    zans[0], $
    tsobj[0], $
    name='REDSTRUCT' )
   redstruct = replicate(redstruct, n_elements(tsobj))
   redstruct.modelcolor = modelcolor
   redstruct.dinv4000 = dinv4000
   redstruct.ew_halpha = ew_halpha
   redstruct.ew_mg2 = ew_mg2
   redstruct[iselect].redsample = 1
   struct_assign, zans, redstruct, /nozero
   struct_assign, tsobj, redstruct, /nozero

   return, redstruct
end
;------------------------------------------------------------------------------
