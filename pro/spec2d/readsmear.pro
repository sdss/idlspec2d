pro readsmear, thisplate, thismjd, fiberid, $
 smearflux=smearflux, smearivar=smearivar

   if (NOT keyword_set(fiberid)) then fiberid = lindgen(640) + 1

   platestr = string(thisplate, format='(i4.4)')
   mjdstr = string(thismjd, format='(i5.5)')
   platemjd = platestr + '-' + mjdstr
   thisdir = concat_dir(getenv('SPECTRO_DATA'), platestr)

   ; Read the wavelength mapping + modelled flux
   readspec, thisplate, thismjd, fiberid, loglam=loglam, synflux=synflux

   ; Figure out the smear exposure number for this plate
   thisplan = 'spPlancomb-' + platemjd +'.par'
   yanny_read, filepath(thisplan, root_dir=thisdir), pp
   spexp = *pp[0]
   yanny_free, pp
   ii = where(spexp.flavor EQ 'smear')
   if (ii[0] EQ -1) then return
   smearnames = spexp[ii[0]].name ; Select the first smear if more than one.
   smearcam = strmid(smearnames, 8, 2)

   for specid=0, 1 do begin
      findx = fiberid - 1 - 320*specid ; Index numbers in these files
      ii = where(findx GE 0 AND findx LT 320, ct)
      if (ct GT 0) then begin
         findx = findx[ii]

         if (specid EQ 0) then camnames = ['b1','r1'] $
          else camnames = ['b2','r2']

         for icam=0, 1 do begin

            ;----------
            ; Read the smear for this camera

            fmin = min(findx, max=fmax)
            range = [fmin,fmax]
            thissmear = smearnames[ where(smearcam EQ camnames[icam]) ]
            tempflux = mrdfits(filepath(thissmear, root_dir=thisdir), 0, $
             range=range)
            tempivar = mrdfits(filepath(thissmear, root_dir=thisdir), 1, $
             range=range)
            tempwset = mrdfits(filepath(thissmear, root_dir=thisdir), 3)
            traceset2xy, tempwset, xx, temploglam
            tempflux = tempflux[*,findx-fmin]
            tempivar = tempivar[*,findx-fmin]
            temploglam = temploglam[*,findx]

            ;----------
            ; Read and apply the flux-calibration vector

            fcalibprefix = 'spFluxcalib-' + platemjd
            fcalibfile = djs_filepath(fcalibprefix + $
             '-' + camnames[icam] + '.fits', root_dir=thisdir)

            calibhdr = headfits(fcalibfile)
            cwavemin = sxpar(calibhdr, 'WAVEMIN')
            cwavemax = sxpar(calibhdr, 'WAVEMAX')
            calibset = mrdfits(fcalibfile, 1)
            calibfac = bspline_valu(temploglam, calibset)

            ; Set to bad any pixels whose wavelength is outside the known
            ; flux-calibration region.
            ibad = where(temploglam LT alog10(cwavemin) $
             OR temploglam GT alog10(cwavemax))
            if (ibad[0] NE -1) then calibfac[ibad] = 0

            divideflat, tempflux, invvar=tempivar, calibfac, $
             minval=0.05*mean(calibfac)
;            temppixmask = temppixmask $
;             OR (calibfac LE 0.05*mean(calibfac)) $
;             * pixelmask_bits('BADFLUXFACTOR')

stop
         endfor ; End loop over camera b/r
      endif
   endfor ; End loop over spectrograph 1/2

stop
   return
end

; Here's an example of getting the flux-correction image for a b1 frame
; of plate 585:
; 
; filename='spFrame-b1-00009835.fits.gz'
; corrfile='spFluxcorr-00009835-1.fits'
; tempwset = mrdfits(filename, 3)
; corrset = mrdfits(corrfile, 1)
; traceset2xy, corrset, tempwave, corrimg
