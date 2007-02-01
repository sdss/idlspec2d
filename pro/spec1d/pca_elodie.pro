;------------------------------------------------------------------------------
; Generate the eigen-vectors from the Elodie eschelle data,
; for use in fitting velocity dispersions.
;------------------------------------------------------------------------------
pro pca_elodie

   outfile = 'spEigenElodie.fits'

   ;----------
   ; Read all the Elodie spectra

   allfiles = elodie_filelist(minwave=minwave)
   nfile = n_elements(allfiles)
   
   t0 = systime(1)
   for ifile=0, nfile-1 do begin
      splog, 'Reading file ', ifile+1, ' of ', nfile
      thisflux = read_elodie(allfiles[ifile], loglam=loglam, minloglamclip=alog10(minwave))
      if (NOT keyword_set(allflux)) then allflux = thisflux $
       else allflux = [[allflux],[thisflux]]
   endfor
   npix = n_elements(loglam)
   splog, 'Time to read all files = ', systime(1) - t0

   ;----------
   ; Trim wavelengths to those covered by the majority of the objects

   fracgpix = total(allflux NE 0, 2) / nfile
   igood = where(fracgpix GT 0.95)
   i1 = igood[0]
   i2 = (reverse(igood))[0]
   loglam = loglam[i1:i2]
   allflux = allflux[i1:i2,*]

   ;----------
   ; Interpolate over bad data, of which there is very little

   allflux = djs_maskinterp(allflux, allflux EQ 0, /const, iaxis=0)

   ;----------
   ; Compute something like the equiv. width for H-alpha, so that we
   ; can reject active stars

   lwav = 6564.6
   i1 = where(loglam GT alog10(lwav-5) AND loglam LT alog10(lwav+5), num1)
   i2 = where(loglam GT alog10(lwav-10) AND loglam LT alog10(lwav+10), num2)

   sum1 = total(allflux[i1,*],1) / num1
   sum2 = total(allflux[i2,*],1) / num2
   fline = sum1 / (sum2 + (sum2 LE 0))

   ;----------
   ; Compute the wavelength coverage for each object

   fraccov = total(allflux NE 0,1) / npix

   ;----------
   ; Select objects for PCA

   istar = where(fline LT 1.1 AND fraccov GT 0.95, nstar)
   splog, 'Selecting ', nstar, ' objects of ', nfile
   t0 = systime(1)
   pres = pcomp(transpose(allflux[*,istar]), eigenval=eigenval, /double)
   splog, 'Time to compute PCA = ', systime(1) - t0
   pres = transpose(pres)

   ;----------
   ; Keep only the most significant eigenvectors

   ikeep = where(eigenval GT 0.01)
   eigenval = eigenval[ikeep]
   pres = float( pres[*,ikeep] )

   ;----------
   ; Write output file

   dloglam = loglam[1] - loglam[0]
   sxaddpar, hdr, 'OBJECT', 'GALAXY'
   sxaddpar, hdr, 'COEFF0', loglam[0]
   sxaddpar, hdr, 'COEFF1', dloglam
   for i=0, n_elements(eigenval)-1 do $
    sxaddpar, hdr, 'EIGEN'+strtrim(string(i),1), eigenval[i]

   mwrfits, pres, outfile, hdr, /create

   return
end
;------------------------------------------------------------------------------
