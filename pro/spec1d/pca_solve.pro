;------------------------------------------------------------------------------
function pca_solve, objflux, objivar, objloglam, zfit, $
 wavemin=wavemin, wavemax=wavemax, $
 niter=niter, nkeep=nkeep, newloglam=newloglam, eigenval=eigenval

   if (NOT keyword_set(niter)) then niter = 10
   if (NOT keyword_set(nkeep)) then nkeep = 3

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then nobj = 1 $
    else nobj = dims[1]
   objdloglam = objloglam[1] - objloglam[0]

   splog, 'Building PCA from ', nobj, ' object spectra'

   ;----------
   ; The redshift of each object in pixels would be LOGSHIFT/OBJDLOGLAM

   logshift = alog10(1.d + zfit)

   ;----------
   ; Determine the new wavelength mapping

   logmin = min(objloglam) - max(logshift)
   logmax = max(objloglam) - min(logshift)
   if (keyword_set(wavemin)) then logmin = logmin > alog10(wavemin)
   if (keyword_set(wavemax)) then logmax = logmax < alog10(wavemax)
   newloglam = wavevector(logmin, logmax, binsz=objdloglam)
   nnew = n_elements(newloglam)
   newflux = fltarr(nnew,nobj)
   newivar = fltarr(nnew,nobj)

   ;----------
   ; Shift each spectra to z=0

   for iobj=0, nobj-1 do begin
      indx = where(objloglam[*,iobj] GT 0)
print,'OBJECT ',iobj
      combine1fiber, objloglam[indx,iobj]-logshift[iobj], $
       objflux[indx,iobj], objivar[indx,iobj], $
       newloglam=newloglam, binsz=objdloglam, newflux=flux1, newivar=ivar1
      newflux[*,iobj] = flux1
      newivar[*,iobj] = ivar1
   endfor

   ;----------
   ; Construct the synthetic weight vector, to be used when replacing
   ; the low-S/N object pixels with the reconstructions.

   synwvec = fltarr(nnew) + 1 ; Set to 1 if no data for this wavelength
   for ipix=0, nnew-1 do begin
      indx = where(newivar[ipix,*] NE 0)
      if (indx[0] NE -1) then $
       synwvec[ipix] = mean(newivar[ipix,indx])
   endfor

   ;----------
   ; Compute a mean spectrum, and use this to replace masked pixels.
   ; Use only the NUSE spectra with flux levels at least 5% of the median
   ; flux level.  For wavelengths with no unmasked data in any spectrum,
   ; just average all the spectra for lack of anything better to do.

   normflux = total(newflux,1) / nnew
   iuse = where(normflux GT 0.05 * median(normflux), nuse)
   synflux = fltarr(nnew)
   for ipix=0, nnew-1 do begin
      ibad = where(newivar[ipix,iuse] EQ 0, nbad)
      if (nbad LT nuse) then begin
         synflux[ipix] = total( newflux[ipix,iuse] * newivar[ipix,iuse]) $
          / total(newivar[ipix,iuse] * normflux[iuse])
      endif else begin
         synflux[ipix] = total( newflux[ipix,iuse] / normflux[iuse]) / nuse
      endelse
   endfor

   for iobj=0, nobj-1 do begin
      ibad = where(newivar[*,iobj] EQ 0)
      if (ibad[0] NE -1) then $
       newflux[ibad,iobj] = synflux[ibad] * normflux[iobj]
   endfor

   ;----------
   ; If there is only 1 object spectrum, then all we can do is return it
   ; (after it has been re-binned).

   if (nobj EQ 1) then begin
      eigenval = 1.0
      return, newflux
   endif

   ;----------
   ; Iteratively do the PCA solution

   filtflux = newflux

   t0=systime(1)
   for iiter=0, niter-1 do begin
      pres = pcomp(transpose(filtflux), coeff=coeff, eigenval=eigenval,/double)

      sqivar = sqrt(newivar)
      bvec = filtflux * sqivar
      mmatrix = pres[0:nkeep-1,*]
      for i=0, nkeep-1 do $
       mmatrix[i,*] = mmatrix[i,*] * sqivar
      mmatrixt = transpose(mmatrix)

      for iobj=0, nobj-1 do begin
         junk = computechi2(newflux[*,iobj], sqrt(newivar[*,iobj]), $
          transpose(pres[0:nkeep-1,*]), acoeff=theta)
         synflux = theta # pres[0:nkeep-1,*]
         filtflux[*,iobj] = (newivar[*,iobj] * newflux[*,iobj] + $
          synwvec * synflux) / (newivar[*,iobj] + synwvec)
;splot,filtflux[*,iobj]
;soplot,synflux,color='red'
      endfor

writefits, 'test-'+strtrim(string(iiter),1)+'.fits', $
 float(transpose(pres[0:nkeep-1,*]))
      splog, 'Elapsed time for iteration #', iiter, ' = ', systime(1)-t0
   endfor

   eigenval = eigenval[0:nkeep-1]
   return, transpose(pres[0:nkeep-1,*])
end
;------------------------------------------------------------------------------
