;------------------------------------------------------------------------------
pro zeigenq, platefile

plate = 306
djs_readcol, '/home/schlegel/idlspec2d/etc/regress1d_all.dat', $
 chicplate, junk, chicfiberid, zfit, chicclass, format='(L,L,L,F,A)'
ii=where(chicplate EQ plate AND chicclass EQ 'QSO')
zfit=zfit[ii]
fiber=chicfiberid[ii]

   ;----------
   ; Read the 2D output file

   readspec, plate, fiber, mjd=mjd, flux=objflux, invvar=objivar, $
    andmask=andmask, plugmap=plugmap, loglam=objloglam

   dims = size(objflux, /dimens)
   npix = dims[0]
   nobj = dims[1]
   objdloglam = objloglam[1] - objloglam[0]

objivar = objivar < 10.0 ; ??? Set limit for maximum inverse-variance

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-subtraction
   ; residuals.  Grow that mask by 2 pixels in each direction.

   skymask = (andmask AND pixelmask_bits('BRIGHTSKY')) NE 0
   for iobj=0, nobj-1 do $
    skymask[*,iobj] = smooth(float(skymask[*,iobj]),5) GT 0
   ibad = where(skymask)
andmask = 0 ; Free memory
   if (ibad[0] NE -1) then objivar[ibad] = 0

   ;----------
   ; The redshift of each object in pixels would be LOGSHIFT/OBJDLOGLAM

   logshift = alog10(1.d + zfit)
;chicpix = alog10(1.d + zfit) / objdloglam

   ;----------
   ; Determine the new wavelength mapping

   newloglam = wavevector(min(objloglam)-max(logshift), $
    max(objloglam)-min(logshift), binsz=objdloglam)
;newloglam=newloglam[1000:4000] ; ???
   nnew = n_elements(newloglam)
   newflux = fltarr(nnew,nobj)
   newivar = fltarr(nnew,nobj)

   ;----------
   ; Shift to z=0

   for iobj=0, nobj-1 do begin
print,'OBJECT ',iobj
      combine1fiber, objloglam[*,iobj]-logshift[iobj], $
       objflux[*,iobj], objivar[*,iobj], $
       newloglam=newloglam, binsz=objdloglam, newflux=flux1, newivar=ivar1
      newflux[*,iobj] = flux1
      newivar[*,iobj] = ivar1
   endfor

; Interpolate over bright emission lines
;bmask = (10^newloglam GT 6535 AND 10^newloglam LT 6600) $
;     OR (10^newloglam GT 6700 AND 10^newloglam LT 6745)
;for iobj=0, nobj-1 do $
; newflux[*,iobj] = djs_maskinterp(newflux[*,iobj], bmask, /const)

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
stop

   ;----------

niter = 10
nkeep = 4 ; Number of modes to keep in replacement
   filtflux = newflux

   for iiter=0, niter-1 do begin
t0=systime(1)
      pres = pcomp(transpose(filtflux), coeff=coeff, eigenval=eigenval,/double)
print,systime(1)-t0

      sqivar = sqrt(newivar)
      bvec = filtflux * sqivar
      mmatrix = pres[0:nkeep-1,*]
      for i=0, nkeep-1 do $
       mmatrix[i,*] = mmatrix[i,*] * sqivar
      mmatrixt = transpose(mmatrix)

      for iobj=0, nobj-1 do begin
         zz = computez(newflux[*,iobj], newivar[*,iobj], $
          transpose(pres[0:nkeep-1,*]), 0*transpose(pres[0:nkeep-1,*])+1, $
          zoffset=0, zmin=0, zmax=0, nfind=1)
         synflux = zz.theta # pres[0:nkeep-1,*]
         filtflux[*,iobj] = (newivar[*,iobj] * newflux[*,iobj] + $
          synwvec * synflux) / (newivar[*,iobj] + synwvec)
;splot,filtflux[*,iobj]
;soplot,synflux,color='red'
      endfor
   endfor

save,file='eigen.ss'
stop
   return
end
;------------------------------------------------------------------------------
