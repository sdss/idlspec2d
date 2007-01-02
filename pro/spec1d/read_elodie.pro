;------------------------------------------------------------------------------
; Read one Elodie echelle spectrum and put on the SDSS pixelization.
; The spectra are returned very crudely normalized.
; 
;------------------------------------------------------------------------------
function read_elodie, filename, loglam=newloglam, hdr=hdr, minloglamclip=minloglamclip

   if (n_params() LT 1) then begin
      print, 'Syntax: flux = read_elodie(filename, [loglam=, hdr=, minloglamclip= ])'
      return, 0
   endif

   if (not keyword_set(dospline)) then dospline=0

   ; Version 2 ELODIE spectra have 13501 0.2-Ang pixels from 4100 - 6800 Ang
   ; Version 3 ELODIE spectra have 14001 0.2-Ang pixels from 4000 - 6800 Ang
   ; There doesn't seem to a a version card, so use NAXIS1
   ;
   objflux = mrdfits(filename, 0, hdr)
   if (NOT keyword_set(objflux)) then $
    message, 'File not found: ' + filename
   objwave = sxpar(hdr,'CRVAL1') + dindgen(sxpar(hdr,'NAXIS1')) $
    * sxpar(hdr,'CDELT1')

   dloglam = 1.d-4
   nsubsamp = 8
   ; Keep perfect backward compatibility
   if (sxpar(hdr,'NAXIS1') EQ 13501) then begin
       loglam0 = 3.6120
       nnew = 2220L
       dospline=0
   end else begin
       loglam0 = alog10(sxpar(hdr,'CRVAL1'))
       loglam1 = alog10(sxpar(hdr,'CRVAL1') + $
                        sxpar(hdr,'NAXIS1') * sxpar(hdr,'CDELT1'))
       nnew = long((loglam1 - loglam0) / dloglam)
       dospline=1               ; There are many masked pixels in Elodie3 spectra; enough
                                ; to seriously nick (and effectively
                                ; shift) real features.
                                ; See 5890 in 00001, say.
   end


   sigma = 1.0 ; Gaussian sigma in pixels on the final sampled vector


   ; Convert to heliocentric wavelengths in vacuum
   airtovac, objwave
; It looks like the spectrum is already de-redshifted, so don't do the below
;   cspeed = 2.99792458d5
;   vhelio = sxpar(hdr,'VR')
;   objwave = objwave * (1.d0 + vhelio / cspeed)

   ; Resample to a vector starting at bigloglam0,
   ; spaced every dloglam/nsubsamp
;   bigloglam0 = loglam0 - 0.5 * (nsubsamp - 1) * dloglam / nsubsamp
;   ibigpix = (alog10(objwave) - bigloglam0) / (dloglam/nsubsamp)
;   newflux = fltarr(nnew*nsubsamp)
;   igood = where(finite(objflux))
;   populate_image, newflux, ibigpix[igood], weights=objflux[igood], $
;    assign='cic'

   igood = where(finite(objflux))
   bigloglam0 = loglam0 - 0.5 * (nsubsamp - 1) * dloglam / nsubsamp
   bigloglam = bigloglam0 + lindgen(nnew*nsubsamp) * (dloglam/nsubsamp)
   newflux = interpol(objflux[igood], alog10(objwave[igood]), bigloglam, spline=dospline)

   ; Convolve with a gaussian
   nkpix = long(10*sigma*nsubsamp)
   kern = exp( -0.5 * (findgen(nkpix*2+1) - nkpix)^2 / (sigma*nsubsamp)^2 )
   kern = kern / total(kern)
   newflux = convol(newflux, kern, /center)

   ; Now rebin to the final sampling
   newflux = rebin(newflux, nnew) * nsubsamp
   newloglam = loglam0 + lindgen(nnew) * dloglam

   ; Set any pixels outside of the wavelength range to zero.
   ; Insist that we be at least 3*sigma from the enpoint measurements.
   minloglam = alog10( min(objwave[igood]) ) + 3 * sigma * dloglam
   maxloglam = alog10( max(objwave[igood]) ) - 3 * sigma * dloglam

                                ; Elodie 2 spectra go partway down into
                                ; H-delta. Let the caller optionally
                                ; avoid that.
   if keyword_set(minloglamclip) then $
      minloglam = minloglam > minloglamclip

   newflux = newflux * (newloglam GT minloglam) * (newloglam LT maxloglam)

                                ; In Elodie 2 the spectra are
                                ; normalized to the flux between 5550-5600.
                                ; Elodie 3 spectra are not. There does
                                ; not appear to be a card describing
                                ; this. BUNIT might work but I don't
                                ; trust it to be constant enough.
   if (sxpar(hdr,'NAXIS1') NE 13501) then begin
       if (sxpar(hdr, 'BUNIT') NE 'erg/cm^2/s/0.1nm') then $
         splog, 'unexpected value for BUNIT: ', sxpar(hdr, 'BUNIT')

                                ; Approximate what Elodie 2 claims is done.
       mflux_i = where(newloglam GE alog10(5550) and newloglam LE alog10(5600))
       newflux /= mean(newflux[mflux_i])
   end

   return, newflux
end
;------------------------------------------------------------------------------
