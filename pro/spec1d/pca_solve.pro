;+
; NAME:
;   pca_solve
;
; PURPOSE:
;   Iteratively find PCA solution for noisy or gappy spectra.
;
; CALLING SEQUENCE:
;   res = pca_solve( objflux, objivar, objloglam, [ zfit, $
;    wavemin=, wavemax=, newloglam=, $
;    maxiter=, niter=, nkeep=, nreturn=, eigenval=, acoeff=, outmask=, $
;    usemask=, _EXTRA= ] )
;
; INPUTS:
;   objflux        - Object fluxes [NPIX,NSPEC]
;   objivar        - Object inverse variances [NPIX,NSPEC]
;
; OPTIONAL INPUTS:
;   objloglam      - Object wavelengths in log10(Angstroms)
;                    [NPIX] if the same wavelength mapping for all spectra,
;                    or [NPIX,NSPEC] if the wavelength mappings are different.
;   zfit           - Redshifts of each input spectrum [NSPEC]; if set, then
;                    each input spectrum is de-redshifted to z=0.
;   wavemin        - Minimum wavelength to use in PCA solution, in Angstroms;
;                    default to the minimum (de-redshifted) input wavelength.
;   wavemax        - Maximum wavelength to use in PCA solution, in Angstroms
;                    default to the minimum (de-redshifted) input wavelength.
;   newloglam      - PCA wavelength sampling in log-10(Angstroms) [NNEWPIX]
;   maxiter        - Number of rejection iterations; default to 0 (no rejection)
;   niter          - Number of PCA iterations; default to 10.
;   nkeep          - Number of PCA components to keep in each iteration
;                    and use in replacing noisy or missing data; default to 3.
;   nreturn        - Number of PCA components to return; default to the same as
;                    NKEEP.
;   _EXTRA         - Keywords for DJS_REJECT().
;
; OUTPUTS:
;   res            - PCA spectra in rest-frame [NNEWPIX,NKEEP]
;
; OPTIONAL OUTPUTS:
;   newloglam      - PCA wavelength sampling in log-10(Angstroms) [NNEWPIX]
;   newflux        - Rebinned OBJFLUX on the wavelength-mapping NEWLOGLAM.
;   newivar        - Rebinned OBJIVAR on the wavelength-mapping NEWLOGLAM.
;   eigenval       - Eigenvalue for each output eigenspectra [NRETURN]
;   acoeff         - PCA coefficients [NRETURN,NOBJ]
;   outmask        - Output mask from DJS_REJECT() [NNEWPIX,NOBJ]
;   usemask        - Number of unmasked spectra used for each pixel, so these
;                    are integers in the range 0 to NSPEC [NNEWPIX]; this is
;                    equivalent to TOTAL(OUTMASK,2).
;
; COMMENTS:
;   The best-fit eigenspectra for each of the input spectra can be determined
;   for object number IOBJ by ACOEFF[*,IOBJ] # RES.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   combine1fiber
;   computechi2()
;   djs_mean()
;   djs_reject()
;   splog
;   wavevector()
;
; REVISION HISTORY:
;   10-Oct-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function pca_solve, objflux, objivar, objloglam, zfit, $
 wavemin=wavemin, wavemax=wavemax, newloglam=newloglam, $
 newflux=newflux, newivar=newivar,  maxiter=maxiter, $
 niter=niter, nkeep=nkeep, nreturn=nreturn, eigenval=eigenval, acoeff=acoeff, $
 outmask=outmask, usemask=usemask, _EXTRA=KeywordsForReject

   if (n_elements(maxiter) EQ 0) then maxiter = 0
   if (NOT keyword_set(niter)) then niter = 10
   if (NOT keyword_set(nkeep)) then nkeep = 3
   if (NOT keyword_set(nreturn)) then nreturn = nkeep

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then nobj = 1 $
    else nobj = dims[1]

   splog, 'Building PCA from ', nobj, ' object spectra'

   ;----------
   ; The redshift of each object in pixels would be LOGSHIFT/OBJDLOGLAM

   if (keyword_set(zfit)) then $
    logshift = alog10(1.d + zfit) $
   else $
    logshift = fltarr(nobj)

   ;----------
   ; Determine the new wavelength mapping

if (keyword_set(objloglam)) then begin ; ???

   if (NOT keyword_set(newloglam)) then begin
      igood = where(objloglam NE 0)
      objdloglam = abs(objloglam[1] - objloglam[0])
      logmin = min(objloglam[igood]) - max(logshift)
      logmax = max(objloglam[igood]) - min(logshift)
      if (keyword_set(wavemin)) then logmin = logmin > alog10(wavemin)
      if (keyword_set(wavemax)) then logmax = logmax < alog10(wavemax)
      newloglam = wavevector(logmin, logmax, binsz=objdloglam)
   endif else begin
      objdloglam = abs(newloglam[1] - newloglam[0])
   endelse
   nnew = n_elements(newloglam)
   newflux = fltarr(nnew,nobj)
   newivar = fltarr(nnew,nobj)

   ndim = size(objloglam, /n_dimen)
   if (ndim EQ 1) then begin
      qwavevec = 0B
   endif else begin
      qwavevec = 1B
      if ((size(objloglam, /dimens))[1] NE nobj) then $
       message, 'Wrong number of dimensions for OBJLOGLAM'
   endelse

   ;----------
   ; Shift each spectra to z=0 and sample at the output wavelengths

   for iobj=0, nobj-1 do begin
      indx = where(objloglam[*,iobj*qwavevec] GT 0)
print,'OBJECT ',iobj
      combine1fiber, objloglam[indx,iobj*qwavevec]-logshift[iobj], $
       objflux[indx,iobj], objivar[indx,iobj], $
       newloglam=newloglam, binsz=objdloglam, newflux=flux1, newivar=ivar1
      newflux[*,iobj] = flux1
      newivar[*,iobj] = ivar1
   endfor

endif else begin
   newflux = objflux
   newivar = objivar
   nnew = (size(objflux,/dimens))[0]
endelse

   ;----------
   ; Construct the synthetic weight vector, to be used when replacing
   ; the low-S/N object pixels with the reconstructions.

   synwvec = fltarr(nnew) + 1 ; Set to 1 if no data for this wavelength
   for ipix=0, nnew-1 do begin
      indx = where(newivar[ipix,*] NE 0)
      if (indx[0] NE -1) then $
       synwvec[ipix] = djs_mean(newivar[ipix,indx])
   endfor

   ;----------
   ; Compute a mean spectrum, and use this to replace masked pixels.
   ; Use only the NUSE spectra with flux levels at least 5% of the median
   ; flux level.  For wavelengths with no unmasked data in any spectrum,
   ; just average all the spectra for lack of anything better to do.

;   normflux = total(newflux,1) / nnew
;   iuse = where(normflux GT 0.05 * median(normflux), nuse)
;   synflux = fltarr(nnew)
;   usemask = lonarr(nnew)
;   for ipix=0, nnew-1 do begin
;      ibad = where(newivar[ipix,iuse] EQ 0, nbad)
;      usemask[ipix] = nuse - nbad
;      if (nbad LT nuse) then begin
;         synflux[ipix] = total( newflux[ipix,iuse] * newivar[ipix,iuse]) $
;          / total(newivar[ipix,iuse] * normflux[iuse])
;      endif else begin
;         synflux[ipix] = total( newflux[ipix,iuse] / normflux[iuse]) / nuse
;      endelse
;   endfor
;
;   for iobj=0, nobj-1 do begin
;      ibad = where(newivar[*,iobj] EQ 0)
;      if (ibad[0] NE -1) then $
;       newflux[ibad,iobj] = synflux[ibad] * normflux[iobj]
;   endfor

; Construct the USEMASK from the output mask (OUTMASK) instead of from NEWIVAR.
;   if (nobj EQ 1) then $
;    usemask = newivar NE 0 $
;   else $
;    usemask = total(newivar NE 0, 2)

   ;----------
   ; If there is only 1 object spectrum, then all we can do is return it
   ; (after it has been re-binned).

   if (nobj EQ 1) then begin
      if (arg_present(eigenval)) then eigenval = 1.0
      if (arg_present(acoeff)) then acoeff = 1.0
      if (arg_present(outmask)) then outmask = lonarr(nnew) + 1
      if (arg_present(usemask)) then usemask = lonarr(nnew) + 1
      return, newflux
   endif

   ;----------
   ; Rejection iteration loop

   qdone = 0
   iiter = 0
   ; Begin with all points good (unless the inverse variance is zero).
;   outmask = make_array(dimension=size(newflux,/dimens), /byte) + 1B
   outmask = 0
   inmask = newivar NE 0

   while ((qdone EQ 0) AND (iiter LE maxiter)) do begin

      qdone = djs_reject(newflux, ymodel, inmask=inmask, outmask=outmask, $
       invvar=newivar, _EXTRA=KeywordsForReject)

      ;----------
      ; Iteratively do the PCA solution

      filtflux = newflux
      acoeff = fltarr(nkeep,nobj)

      t0=systime(1)
      for ipiter=0, niter-1 do begin
         eigenval = 1 ; Set so that the PCOMP() routine returns this.
         coeff = 1 ; Set so that the PCOMP() routine returns this.

         totflux = fltarr(nobj)
         for iobj=0, nobj-1 do $
          totflux[iobj] = total(abs(filtflux[*,iobj] - filtflux[0,iobj]))
         igoodobj = where(totflux GT 0, ngoodobj)
         if (ngoodobj EQ nobj) then begin
            pres = pcomp(transpose(filtflux), eigenval=eigenval, /double)
         endif else begin
            tmp_pres = pcomp(transpose(filtflux[*,igoodobj]), $
             eigenval=tmp_eigenval, /double)
            pres = dblarr(nobj,nnew)
            pres[igoodobj,*] = tmp_pres
            eigenval = dblarr(1,nobj)
            eigenval[0,igoodobj] = tmp_eigenval
         endelse

         maskivar = newivar * outmask
         sqivar = sqrt(maskivar)
         bvec = filtflux * sqivar
         mmatrix = pres[0:nkeep-1,*]
         for i=0, nkeep-1 do $
          mmatrix[i,*] = mmatrix[i,*] * sqivar
         mmatrixt = transpose(mmatrix)

         for iobj=0, nobj-1 do begin
            junk = computechi2(newflux[*,iobj], sqivar[*,iobj], $
             transpose(pres[0:nkeep-1,*]), acoeff=theta)
            synflux = theta # pres[0:nkeep-1,*]
            filtflux[*,iobj] = (maskivar[*,iobj] * newflux[*,iobj] + $
             synwvec * synflux) / (maskivar[*,iobj] + synwvec)
            acoeff[*,iobj] = theta
;splot,filtflux[*,iobj]
;soplot,synflux,color='red'
         endfor

;writefits, 'test-'+strtrim(string(ipiter),1)+'.fits', $
; float(transpose(pres[0:nkeep-1,*]))
         splog, 'Elapsed time for iteration #', ipiter, ' = ', systime(1)-t0
      endfor ; End PCA iterations

      ;----------
      ; Now set YMODEL for rejecting points

      ymodel = acoeff ## transpose(pres[0:nkeep-1,*])

      iiter = iiter + 1
   endwhile ; End rejection iterations

   if (arg_present(usemask)) then begin
      if (nobj EQ 1) then $
       usemask = outmask $
      else $
       usemask = total(outmask, 2)
   endif

   eigenval = eigenval[0:nreturn-1]
   return, transpose(pres[0:nreturn-1,*])
end
;------------------------------------------------------------------------------
