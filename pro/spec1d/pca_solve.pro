;------------------------------------------------------------------------------
function pca_solve, objflux, objivar, objloglam, zfit, $
 niter=niter, nkeep=nkeep, newloglam=newloglam, eigenval=eigenval

   if (NOT keyword_set(niter)) then niter = 10
   if (NOT keyword_set(nkeep)) then nkeep = 3

   dims = size(objflux, /dimens)
   npix = dims[0]
   nobj = dims[1]
   objdloglam = objloglam[1] - objloglam[0]

   splog, 'Building PCA from ', nobj, ' object spectra'

   ;----------
   ; The redshift of each object in pixels would be LOGSHIFT/OBJDLOGLAM

   logshift = alog10(1.d + zfit)

   ;----------
   ; Determine the new wavelength mapping

   newloglam = wavevector(min(objloglam) - max(logshift), $
    max(objloglam) - min(logshift), binsz=objdloglam)
;newloglam=newloglam[1000:4000] ; ???
   nnew = n_elements(newloglam)
   newflux = fltarr(nnew,nobj)
   newivar = fltarr(nnew,nobj)

   ;----------
   ; Shift each spectra to z=0

   for iobj=0, nobj-1 do begin
print,'OBJECT ',iobj
      combine1fiber, objloglam[*,iobj]-logshift[iobj], $
       objflux[*,iobj], objivar[*,iobj], $
       newloglam=newloglam, binsz=objdloglam, newflux=flux1, newivar=ivar1
      newflux[*,iobj] = flux1
      newivar[*,iobj] = ivar1
   endfor

   ;----------
   ; Construct the synthetic weight vector

   synwvec = fltarr(nnew) + 1 ; Set to 1 if no data for this wavelength
   for ipix=0, nnew-1 do begin
      indx = where(newivar[ipix,*] NE 0)
      if (indx[0] NE -1) then $
       synwvec[ipix] = mean(newivar[ipix,indx])
   endfor

   ;----------
   ; Compute a mean spectrum, and use this to replace masked pixels.

   normflux = total(newflux,1) / nnew
   synflux = fltarr(nnew)
   for ipix=0, nnew-1 do begin
      ibad = where(newivar[ipix,*] EQ 0, nbad)
      if (nbad LT nobj) then begin
         synflux[ipix] = total( newflux[ipix,*] * newivar[ipix,*] / normflux) $
          / total(newivar[ipix,*])
         if (nbad GT 0) then $
          newflux[ipix,ibad] = synflux[ipix] * normflux[ibad]
      endif else begin
         synflux[ipix] = total( newflux[ipix,*] / normflux) / nobj
      endelse
   endfor

   ;----------

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

mwrfits, float(transpose(pres[0:nkeep-1,*])), $
 'test-'+strtrim(string(iiter),1)+'.fits'
      splog, 'Elapsed time for iteration #', iiter, ' = ', systime(1)-t0
   endfor

   eigenval = eigenval[0:nkeep-1]
   return, transpose(pres[0:nkeep-1,*])
end
;------------------------------------------------------------------------------
