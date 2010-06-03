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
;    usemask=, /quiet, _EXTRA= ] )
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
;   quiet          - Minimal output to splog
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
;   07-Jul-2003  added quiet keyword - C. Tremonti
;
; VERSION:
;   $Id$
;
;-
FUNCTION pca_solve, objflux, objivar, objloglam, zfit, $
    wavemin=wavemin, wavemax=wavemax, newloglam=newloglam, $
    newflux=newflux, newivar=newivar,  maxiter=maxiter, $
    niter=niter, nkeep=nkeep, nreturn=nreturn, eigenval=eigenval, $
    acoeff=acoeff, outmask=outmask, usemask=usemask, quiet=quiet, $
    _EXTRA=KeywordsForReject
    ;
    ;
    ;
    IF (N_ELEMENTS(maxiter) EQ 0) THEN maxiter = 0
    IF ~KEYWORD_SET(niter) THEN niter = 10
    IF ~KEYWORD_SET(nkeep) THEN nkeep = 3
    IF ~KEYWORD_SET(nreturn) THEN nreturn = nkeep
    ndim = SIZE(objflux, /N_DIMEN)
    dims = SIZE(objflux, /DIMENS)
    npix = dims[0]
    IF (ndim EQ 1) THEN nobj = 1 ELSE nobj = dims[1]
    splog, 'Building PCA from ', nobj, ' object spectra'
    ;
    ; The redshift of each object in pixels would be LOGSHIFT/OBJDLOGLAM
    ;
    IF KEYWORD_SET(zfit) THEN logshift = ALOG10(1.0D + zfit) $
    ELSE logshift = DBLARR(nobj)
    ;
    ; Determine the new wavelength mapping
    ;
    IF KEYWORD_SET(objloglam) THEN BEGIN ; ???
        IF KEYWORD_SET(newloglam) THEN $
            objdloglam = ABS(newloglam[1] - newloglam[0]) $
        ELSE BEGIN
            igood = WHERE(objloglam NE 0)
            objdloglam = ABS(objloglam[1] - objloglam[0])
            logmin = MIN(objloglam[igood]) - MAX(logshift)
            logmax = MAX(objloglam[igood]) - MIN(logshift)
            IF KEYWORD_SET(wavemin) THEN logmin = logmin > ALOG10(wavemin)
            IF KEYWORD_SET(wavemax) THEN logmax = logmax < ALOG10(wavemax)
            newloglam = wavevector(logmin, logmax, binsz=objdloglam)
        ENDELSE
        nnew = N_ELEMENTS(newloglam)
        newflux = DBLARR(nnew,nobj)
        newivar = DBLARR(nnew,nobj)
        ndim = SIZE(objloglam, /N_DIMEN)
        IF (ndim EQ 1) THEN qwavevec = 0B $
        ELSE BEGIN
            qwavevec = 1B
            IF ((SIZE(objloglam, /DIMENS))[1] NE nobj) THEN $
                MESSAGE, 'Wrong number of dimensions for OBJLOGLAM'
        ENDELSE
        ;
        ; Shift each spectra to z=0 and sample at the output wavelengths
        ;
        FOR iobj=0L, nobj-1L DO BEGIN
            indx = WHERE(objloglam[*,iobj*qwavevec] GT 0)
            splog,'OBJECT ',iobj
            combine1fiber, objloglam[indx,iobj*qwavevec]-logshift[iobj], $
                objflux[indx,iobj], objivar[indx,iobj], $
                newloglam=newloglam, binsz=objdloglam, newflux=flux1, $
                newivar=ivar1 ;, /verbose
            newflux[*,iobj] = flux1
            newivar[*,iobj] = ivar1
        ENDFOR
    ENDIF ELSE BEGIN
        newflux = objflux
        newivar = objivar
        nnew = (SIZE(objflux,/DIMENS))[0]
    ENDELSE
    ;
    ; Construct the synthetic weight vector, to be used when replacing
    ; the low-S/N object pixels with the reconstructions.
    ;
    synwvec = DBLARR(nnew) + 1.0D ; Set to 1 if no data for this wavelength
    FOR ipix=0L, nnew-1L DO BEGIN
        indx = WHERE(newivar[ipix,*] NE 0)
        IF (indx[0] NE -1) THEN $
            synwvec[ipix] = djs_mean(newivar[ipix,indx])
    ENDFOR
    ;
    ; Compute a mean spectrum, and use this to replace masked pixels.
    ; Use only the NUSE spectra with flux levels at least 5% of the median
    ; flux level.  For wavelengths with no unmasked data in any spectrum,
    ; just average all the spectra for lack of anything better to do.
    ;
    ; normflux = TOTAL(newflux,1) / nnew
    ; iuse = WHERE(normflux GT 0.05 * MEDIAN(normflux), nuse)
    ; synflux = DBLARR(nnew)
    ; usemask = LONARR(nnew)
    ; FOR ipix=0L, nnew-1L DO BEGIN
    ;     ibad = WHERE(newivar[ipix,iuse] EQ 0, nbad)
    ;     usemask[ipix] = nuse - nbad
    ;     IF (nbad LT nuse) THEN BEGIN
    ;         synflux[ipix] = TOTAL( newflux[ipix,iuse] * newivar[ipix,iuse]) $
    ;             / TOTAL(newivar[ipix,iuse] * normflux[iuse])
    ;     ENDIF ELSE BEGIN
    ;         synflux[ipix] = TOTAL( newflux[ipix,iuse] / normflux[iuse]) / nuse
    ;     ENDELSE
    ; ENDFOR
    ; FOR iobj=0L, nobj-1L DO BEGIN
    ;     ibad = WHERE(newivar[*,iobj] EQ 0)
    ;     IF (ibad[0] NE -1) THEN $
    ;         newflux[ibad,iobj] = synflux[ibad] * normflux[iobj]
    ; ENDFOR
    ;
    ; Construct the USEMASK from the output mask (OUTMASK) instead of
    ; from NEWIVAR.
    ;
    ; IF (nobj EQ 1) THEN usemask = newivar NE 0 $
    ; ELSE usemask = TOTAL(newivar NE 0, 2)
    ;
    ; If there is only 1 object spectrum, then all we can do is return it
    ; (after it has been re-binned).
    ;
    IF (nobj EQ 1) THEN BEGIN
        IF ARG_PRESENT(eigenval) THEN eigenval = 1.0D
        IF ARG_PRESENT(acoeff) THEN acoeff = 1.0D
        IF ARG_PRESENT(outmask) THEN outmask = LONARR(nnew) + 1L
        IF ARG_PRESENT(usemask) THEN usemask = LONARR(nnew) + 1L
        RETURN, newflux
    ENDIF
    ;
    ; Rejection iteration loop
    ;
    qdone = 0
    iiter = 0
    ;
    ; Begin with all points good (unless the inverse variance is zero).
    ;
    ; outmask = MAKE_ARRAY(DIMENSION=SIZE(newflux,/DIMENS), /BYTE) + 1B
    outmask = 0
    inmask = newivar NE 0
    WHILE ((qdone EQ 0) AND (iiter LE maxiter)) DO BEGIN
        qdone = djs_reject(newflux, ymodel, inmask=inmask, outmask=outmask, $
            invvar=newivar, _EXTRA=KeywordsForReject)
        ;
        ; Iteratively do the PCA solution
        ;
        filtflux = newflux
        acoeff = DBLARR(nkeep,nobj)
        t0=SYSTIME(1)
        FOR ipiter=0L, niter-1L DO BEGIN
            eigenval = 1 ; Set so that the PCOMP() routine returns this.
            coeff = 1 ; Set so that the PCOMP() routine returns this.
            totflux = DBLARR(nobj)
            FOR iobj=0L, nobj-1L DO $
                totflux[iobj] = TOTAL(ABS(filtflux[*,iobj] - filtflux[0,iobj]))
            igoodobj = WHERE(totflux GT 0, ngoodobj)
            IF (ngoodobj EQ nobj) THEN BEGIN
                pres = PCOMP(TRANSPOSE(filtflux), EIGENVAL=eigenval, /DOUBLE)
            ENDIF ELSE BEGIN
                tmp_pres = PCOMP(TRANSPOSE(filtflux[*,igoodobj]), $
                    EIGENVAL=tmp_eigenval, /DOUBLE)
                pres = DBLARR(nobj,nnew)
                pres[igoodobj,*] = tmp_pres
                eigenval = DBLARR(1,nobj)
                eigenval[0,igoodobj] = tmp_eigenval
            ENDELSE
            maskivar = newivar * outmask
            sqivar = SQRT(maskivar)
            ; bvec = filtflux * sqivar
            ; mmatrix = pres[0:nkeep-1,*]
            ; FOR i=0L, nkeep-1L DO $
            ;     mmatrix[i,*] = mmatrix[i,*] * sqivar
            ; mmatrixt = TRANSPOSE(mmatrix)
            FOR iobj=0L, nobj-1L DO BEGIN
                chi2 = computechi2(newflux[*,iobj], sqivar[*,iobj], $
                    TRANSPOSE(pres[0:nkeep-1,*]), acoeff=theta)
                synflux = theta # pres[0:nkeep-1,*]
                filtflux[*,iobj] = (maskivar[*,iobj] * newflux[*,iobj] + $
                    synwvec * synflux) / (maskivar[*,iobj] + synwvec)
                acoeff[*,iobj] = theta
                ; splot,filtflux[*,iobj]
                ; soplot,synflux,color='red'
            ENDFOR
            ; writefits, 'test-'+STRTRIM(STRING(ipiter),1)+'.fits', $
            ;     FLOAT(TRANSPOSE(pres[0:nkeep-1,*]))
            IF ~KEYWORD_SET(quiet) THEN $
                splog, 'Elapsed time for iteration #', ipiter, ' = ', $
                    SYSTIME(1)-t0
        ENDFOR ; End PCA iterations
        ;
        ; Now set YMODEL for rejecting points
        ;
        ymodel = acoeff ## TRANSPOSE(pres[0:nkeep-1,*])
        iiter = iiter + 1
    ENDWHILE ; End rejection iterations
    IF ARG_PRESENT(usemask) THEN BEGIN
        IF (nobj EQ 1) THEN usemask = outmask $
        ELSE usemask = TOTAL(outmask, 2)
    ENDIF
    eigenval = eigenval[0:nreturn-1]
    RETURN, TRANSPOSE(pres[0:nreturn-1,*])
END
