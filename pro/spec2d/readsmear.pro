; Return the smear fluxes de-redshifted and rebinned to log-linear.
pro readsmear, thisplate, thismjd, fiberid, $
 smearloglam=smearloglam, smearflux=smearflux, smearivar=smearivar, $
 synflux=synflux

   if (NOT keyword_set(fiberid)) then fiberid = lindgen(640) + 1
smearloglam = 3.5d + lindgen(4500)*1.d-4 ; ???

   platestr = string(thisplate, format='(i4.4)')
   mjdstr = string(thismjd, format='(i5.5)')
   platemjd = platestr + '-' + mjdstr
   thisdir = concat_dir(getenv('SPECTRO_DATA'), platestr)

   ; Read the wavelength mapping + modelled flux
   readspec, thisplate, mjd=thismjd, fiberid, loglam=loglam, zans=zans

   ; Read the synthetic spectrum (on the new wavelength mapping)
;   synflux = multisynthspec(zans, loglam=smearloglam)

   ; Figure out the smear exposure number for this plate
   thisplan = 'spPlancomb-' + platemjd +'.par'
   yanny_read, filepath(thisplan, root_dir=thisdir), pp
   spexp = *pp[0]
   yanny_free, pp
   ii = where(spexp.flavor EQ 'smear')
   if (ii[0] EQ -1) then return
   smearnames = spexp[ii[0]].name ; Select the first smear if more than one.
   smearcam = strmid(smearnames, 8, 2)

   ; Create the arrays in which we'll accumulate the smear spectra
   nfiber = n_elements(fiberid)
   accflux = fltarr(2048,2,nfiber)
   accivar = fltarr(2048,2,nfiber)
   accloglam = fltarr(2048,2,nfiber)
   ioffset = 0

   for specid=0, 1 do begin
      findx = fiberid - 1 - 320*specid ; Index numbers in these files
      ifindx = where(findx GE 0 AND findx LT 320, nfindx)
      if (nfindx GT 0) then begin
         findx = findx[ifindx]

         if (specid EQ 0) then camnames = ['b1','r1'] $
          else camnames = ['b2','r2']

         for icam=0, 1 do begin

            ;----------
            ; Read the smear for this camera

            fmin = min(findx, max=fmax)
            range = [fmin,fmax]
            thissmear = smearnames[ (where(smearcam EQ camnames[icam]))[0] ]
            thissmear = findfile(filepath(thissmear+'*', root_dir=thisdir))
            tempflux = mrdfits(thissmear[0], 0, range=range)
            tempivar = mrdfits(thissmear[0], 1, range=range)
            tempwset = mrdfits(thissmear[0], 3)
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

; Do we want the pixel masks too???
;            temppixmask = temppixmask $
;             OR (calibfac LE 0.05*mean(calibfac)) $
;             * pixelmask_bits('BADFLUXFACTOR')

            ; Read the synthetic spectrum
            for ii=0, nfindx-1 do begin
               tempsynflux = synthspec(zans[ifindx[ii]], loglam=temploglam)
               ; Look for 5-sigma outliers from the synthetic spectrum ???
; Reject neighboring points too???
; And reject points with very low S/N ???
maxdev = 5.0
maxerr = 50.0
minsn = 0.0
               qbad = (abs(tempflux[*,ii] - tempsynflux) $
                * sqrt(tempivar[*,ii]) GT maxdev) $
                OR (abs(tempflux[*,ii]) * sqrt(tempivar[*,ii]) LT minsn) $
                OR tempivar[*,ii] LT 1./maxerr^2
               ibad = where(qbad)
; plot,10^temploglam,tempflux[*,ii]
; djs_oplot,10^temploglam,tempsynflux,color='red'
; djs_oplot,10^temploglam[ibad],0*ibad,color='green',psym=4

               if (ibad[0] NE -1) then $
                tempflux[ibad] = tempsynflux[ibad]
            endfor

            ; Accumulate these flux values into our arrays
            accflux[*,icam,ioffset:ioffset+nfindx-1] = tempflux
            accivar[*,icam,ioffset:ioffset+nfindx-1] = tempivar
            accloglam[*,icam,ioffset:ioffset+nfindx-1] = temploglam
         endfor ; End loop over camera b/r
         ioffset = ioffset + nfindx
      endif
   endfor ; End loop over spectrograph 1/2

   ;----------
   ; Now rebin each spectrum to log-linear and de-redshifted

   npix = n_elements(smearloglam)
   smearflux = fltarr(npix,nfiber)
   smearivar = fltarr(npix,nfiber)
   dloglam = smearloglam[1] - smearloglam[0]

   for ifiber=0, nfiber-1 do begin
print,'OBJECT ',ifiber
      combine1fiber, (accloglam[*,*,ifiber])[*] - alog10(1+zans[ifiber].z), $
       (accflux[*,*,ifiber])[*], (accivar[*,*,ifiber])[*], $
       newloglam=smearloglam, binsz=dloglam, newflux=flux1, newivar=ivar1
      smearflux[*,ifiber] = flux1
      smearivar[*,ifiber] = ivar1
   endfor


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
