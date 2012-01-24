;+
; NAME:
;   pca_star
;
; PURPOSE:
;   Wrapper on pca_solve.pro
;
; CALLING SEQUENCE:
;   pca_star, [inputfile=inputfile], $
;            [binsz=binsz], [niter=niter], [savefile=savefile], [/flux]
;
; INPUTS:
;   None
;
; OPTIONAL INPUTS:
;   inputfile - File containing plate, mjd, fiber, redshift.  If not
;               specified, it will read eigeninput_star.par from
;               $IDLSPEC2D_DIR/templates
;   flux      - Plot input spectra.
;   niter     - Number of iterations to pass to pca_solve, default 10.
;   savefile  - Save the input spectra to a named file for debugging.
;
; OUTPUTS:
;   None, but creates the files spEigenStar-MJD.fits & spEigenStar-MJD.ps
;
; DATA FILES:
;   $IDLSPEC2D_DIR/templates/eigeninput_star.dat
;
; PROCEDURES CALLED:
;   get_juldate
;   dfpsplot
;   yanny_read
;   yanny_free
;   readspec
;   skymask()
;   wavevector()
;   pca_solve()
;   computechi2()
;   djs_median()
;   djs_plot
;   djs_oplot
;   djs_xyouts
;   sxaddpar
;   mwrfits
;   dfpsclose
;
; REVISION HISTORY:
;   Written a long time ago.
;   Updated for BOSS inputs by B. A. Weaver, NYU
;
; VERSION:
;   $Id$
;
;------------------------------------------------------------------------------
PRO pca_star, inputfile=inputfile, flux=flux, niter=niter, savefile=savefile
    ;
    ; Define initial parameters
    ;
    IF ~KEYWORD_SET(inputfile) THEN BEGIN
        inputfile = FILEPATH('eigeninput_star.par', $
            ROOT_DIR=GETENV('IDLSPEC2D_DIR'), SUBDIRECTORY='templates')
    ENDIF
    wavemin = 0
    wavemax = 0
    snmax = 100
    IF ~KEYWORD_SET(niter) THEN niter = 10
    cspeed = 2.99792458d5
    ;
    ; Name the output files
    ;

    get_juldate, jd
    mjdstr = STRING(LONG(jd-2400000.5), FORMAT='(I5)')
    outfile = 'spEigenStar-' + mjdstr + '.fits'
    plotfile = 'spEigenStar-' + mjdstr + '.ps'
    dfpsplot, plotfile, /color, /landscape
    colorvec = ['default', 'red', 'green', 'blue', 'magenta', 'cyan']
    ;
    ; Read the input spectra
    ;
    yanny_read, inputfile, pdat
    slist = *pdat[0]
    yanny_free, pdat
    readspec, slist.plate, slist.fiberid, mjd=slist.mjd, $
        flux=objflux, invvar=objivar, $
        andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=objloglam, /align
    ;
    ; Insist that all of the requested spectra exist
    ;
    imissing = WHERE(plugmap.fiberid EQ 0, nmissing)
    IF (nmissing GT 0) THEN BEGIN
        FOR i=0, nmissing-1 DO $
            PRINT, 'Missing plate=', plate[imissing[i]], $
                ' mjd=', mjd[imissing[i]], $
                ' fiber=', fiber[imissing[i]]
        MESSAGE, STRING(nmissing) + ' missing object(s)'
    ENDIF
    ;
    ; Do not fit where the spectrum may be dominated by sky-sub residuals.
    ;
    objivar = skymask(objivar, andmask, ormask)
    andmask = 0 ; Free memory
    ormask = 0 ; Free memory
    nobj = (SIZE(objflux, /dimens))[1]
    objdloglam = objloglam[1] - objloglam[0]
    IF (KEYWORD_SET(snmax)) THEN BEGIN
        ifix = WHERE(objflux^2 * objivar GT snmax^2)
        IF (ifix[0] NE -1) THEN objivar[ifix] = (snmax/objflux[ifix])^2
    ENDIF
    ;
    ; Find the list of unique star types
    ;
    isort = SORT(slist.class)
    classlist = slist[isort[UNIQ(slist[isort].class)]].class
    ;
    ; LOOP OVER EACH STAR TYPE
    ;
    FOR iclass=0, N_ELEMENTS(classlist)-1 DO BEGIN
        ;
        ; Find the subclasses for this stellar type
        ;
        indx = WHERE(slist.class EQ classlist[iclass], nindx)
        thesesubclass = slist[indx].subclass
        isort = SORT(thesesubclass)
        subclasslist = thesesubclass[isort[UNIQ(thesesubclass[isort])]]
        nsubclass = N_ELEMENTS(subclasslist)
        ;
        ; Solve for 2 eigencomponents if we have specified subclasses
        ; for this stellar type.
        ;
        IF (nsubclass EQ 1) THEN nkeep = 1 ELSE nkeep = 2
        newloglam = objloglam
        pcaflux = pca_solve(objflux[*,indx], objivar[*,indx], objloglam, $
            slist[indx].cz/cspeed, wavemin=wavemin, wavemax=wavemax, $
            niter=niter, nkeep=nkeep, newloglam=newloglam, $
            eigenval=eigenval, acoeff=acoeff, usemask=usemask)
        ;
        ; Interpolate over bad flux values in the middle of a spectrum,
        ; and set fluxes to zero at the blue+red ends of the spectrum.
        ;
        ; minuse = 1 ; ?
        minuse = FLOOR((nindx+1) / 3.)
        qbad = usemask LT minuse
        ;
        ; Interpolate over all bad pixels
        ;
        FOR j=0, nkeep-1 DO $
            pcaflux[*,j] = djs_maskinterp(pcaflux[*,j], qbad EQ 1, /const)
        ;
        ; Set bad pixels at the very start or end of the spectrum to zero instead
        ;
        npix = N_ELEMENTS(qbad)
        IF (qbad[0]) THEN $
            pcaflux[0:(where(qbad EQ 0))[0]-1,*] = 0
        IF (qbad[npix-1]) THEN $
            pcaflux[(REVERSE(WHERE(qbad EQ 0)))[0]+1:npix-1,*] = 0
        ;
        ; The following would plot the 0th object and overplot the best-fit PCA
        ;
        ; ii=0
        ; splot,10^newloglam,objflux[*,indx[ii]]
        ; junk=pcaflux[*,0] * (acoeff[0,ii])[0] + pcaflux[*,1] * (acoeff[1,ii])[0]
        ; soplot,10^newloglam,junk,color='red'
        ;
        ; Re-normalize the first eigenspectrum to a mean of 1
        ;
        norm = MEAN(pcaflux[*,0])
        pcaflux = pcaflux / norm
        acoeff = acoeff * norm
        ;
        ; Now loop through each stellar subclass and reconstruct
        ; an eigenspectrum for that subclass
        ;
        thesesubclassnum = LONARR(N_ELEMENTS(thesesubclass))
        FOR isub=0, nsubclass-1 DO BEGIN
            ii = WHERE(thesesubclass EQ subclasslist[isub])
            thesesubclassnum[ii] = isub
            IF (nkeep EQ 1) THEN BEGIN
                thisflux = pcaflux
            ENDIF ELSE BEGIN
                aratio = acoeff[1,ii] / acoeff[0,ii]
                thisratio = MEDIAN(aratio, /EVEN)
                thisflux = pcaflux[*,0] + thisratio * pcaflux[*,1]
            ENDELSE
            ;
            ; The output wavelength mapping is the same for everything,
            ; so we can simply stack the PCA spectra.
            ;
            IF KEYWORD_SET(fullflux) THEN fullflux = [[fullflux], [thisflux]] $
                ELSE fullflux = thisflux
            IF KEYWORD_SET(namearr) THEN namearr = [namearr, subclasslist[isub]] $
                ELSE namearr = subclasslist[isub]
            plotflux = thisflux / MAX(thisflux) ; Re-scale for plotting
            IF (isub EQ 0) THEN $
                djs_plot, 10^newloglam, plotflux, color=colorvec[0], $
                xtitle='Wavelength [\AA]', ytitle='Flux [arbitrary units]', $
                title='STAR '+classlist[iclass]+': Eigenspectra Reconstructions' $
            ELSE $
                djs_oplot, 10^newloglam, plotflux, $
                color=colorvec[isub MOD N_ELEMENTS(colorvec)]
            nnew = N_ELEMENTS(newloglam)
            XYOUTS, 10^newloglam[nnew-1], plotflux[nnew-1], $
                subclasslist[isub], align=-0.5, $
                color=djs_icolor(colorvec[isub MOD N_ELEMENTS(colorvec)])
        ENDFOR
        IF (nkeep GT 1) THEN BEGIN
            allratio = TRANSPOSE(acoeff[1,*] / acoeff[0,*])
            isort = SORT(thesesubclassnum)
            djs_plot, thesesubclassnum[isort], allratio[isort], ps=-4, $
                xrange=[-1,nsubclass], xstyle=1, xtickname=subclasslist, $
                xtickv=LINDGEN(nsubclass), xticks=nsubclass-1, $
                xtitle='Subclass', ytitle='Eigenvalue Ratio (a_1/a_0)', $
                title='STAR '+classlist[iclass]+': Eigenvalue Ratios'
            FOR j=0, N_ELEMENTS(indx)-1 DO $
                XYOUTS, thesesubclassnum[isort[j]], allratio[isort[j]], $
                    align=0.0, orient=45, $
                    STRING(slist[indx[isort[j]]].plate, slist[indx[isort[j]]].fiberid, $
                    FORMAT='(I4,"-",I3)')
        ENDIF
    ENDFOR
    ;
    ; Construct header for output file
    ;
    sxaddpar, hdr, 'OBJECT', 'STAR'
    sxaddpar, hdr, 'COEFF0', newloglam[0]
    sxaddpar, hdr, 'COEFF1', objdloglam
    sxaddpar, hdr, 'IDLUTILS', idlutils_version(), 'Version of idlutils'
    sxaddpar, hdr, 'SPEC2D', idlspec2d_version(), 'Version of idlspec2d'
    sxaddpar, hdr, 'RUN2D', GETENV('RUN2D'), 'Version of 2d reduction'
    sxaddpar, hdr, 'RUN1D', GETENV('RUN1D'), 'Version of 1d reduction'
    ;
    ; Add a space to the name below, so that 'F' appears as a string and
    ; not as a logical.
    ;
    FOR i=0, N_ELEMENTS(namearr)-1 do $
        sxaddpar, hdr, 'NAME'+STRTRIM(STRING(i),2), namearr[i]+' '
    ;
    ; Write output file
    ;
    mwrfits, FLOAT(fullflux), outfile, hdr, /create
    ;
    ; Create a table of inputs
    ;
    sxaddpar, hdr2, 'FILENAME', inputfile
    inputs0 = {inputs, plate:0L, mjd:0L, fiber:0L, redshift:0.0D}
    inputs = REPLICATE(inputs0,N_ELEMENTS(slist))
    inputs.plate = slist.plate
    inputs.mjd = slist.mjd
    inputs.fiber = slist.fiberid
    inputs.redshift = slist.cz/cspeed
    mwrfits, inputs, outfile, hdr2
    dfpsclose
    RETURN
END
;------------------------------------------------------------------------------
