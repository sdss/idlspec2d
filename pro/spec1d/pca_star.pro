;------------------------------------------------------------------------------
pro spappend, newloglam, pcaflux, fullloglam, fullflux

   if (NOT keyword_set(fullloglam)) then begin
      fullloglam = newloglam
      fullflux = pcaflux
      return
   endif

   npix1 = n_elements(newloglam)
   npix2 = (size(fullloglam,/dimens))[0]

   if (newloglam[0] EQ fullloglam[0]) then begin
      npix = min(npix1,npix2)
      fullloglam = fullloglam[0:npix-1]
      fullflux = [ [fullflux[0:npix-1,*]], [pcaflux[0:npix-1]] ]
   endif else if (newloglam[0] GT fullloglam[0]) then begin
      ; Assume NEWLOGLAM[0] = FULLLOGLAM[PSHIFT], and trim the first PSHIFT
      ; elements from FULLLOGLAM and FULLFLUX.
      junk = min(abs(fullloglam - newloglam[0]), pshift)
      npix = min(npix1,npix2-pshift)
      fullloglam = fullloglam[pshift:pshift+npix-1]
      fullflux = [ [fullflux[pshift:pshift+npix-1,*]], [pcaflux[0:npix-1]] ]
   endif else begin
      ; Assume NEWLOGLAM[PSHIFT] = FULLLOGLAM[0], and trim the first PSHIFT
      ; elements from NEWLOGLAM and PCAFLUX.
      junk = min(abs(newloglam - fullloglam[0]), pshift)
      npix = min(npix1-pshift,npix2)
      fullloglam = fullloglam[0:npix-1]
      fullflux = [ [fullflux[0:npix-1,*]], [pcaflux[pshift:pshift+npix-1]] ]
   endelse

   return
end

;------------------------------------------------------------------------------
pro pca_star

   wavemin = 0
   wavemax = 0
   snmax = 100
   niter = 10
   nkeep = 1

   get_juldate, jd
   mjdstr = string(long(jd-2400000L), format='(i5)')
   outfile = 'spEigenStar-' + mjdstr + '.fits'

   ;----------
   ; Read the input spectra

   eigenfile = filepath('eigeninput_star.dat', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   djs_readcol, eigenfile, skip=2, plate, mjd, fiber, zfit, subclass, $
    format='(L,L,L,D,A)'

   readspec, plate, fiber, mjd=mjd, flux=objflux, invvar=objivar, $
    andmask=andmask, plugmap=plugmap, loglam=objloglam

   nobj = (size(objflux, /dimens))[1]
   objdloglam = objloglam[1] - objloglam[0]

   if (keyword_set(snmax)) then begin
      ifix = where(objflux^2 * objivar GT snmax^2)
      if (ifix[0] NE -1) then objivar[ifix] = (snmax/objflux[ifix])^2
   endif

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
   ; Find the list of unique star types

   isort = sort(subclass)
   classlist = subclass[isort[uniq(subclass[isort])]]

   ;----------
   ; LOOP OVER EACH STAR TYPE

   for iclass=0, n_elements(classlist)-1 do begin

      indx = where(subclass EQ classlist[iclass])

      pcaflux = pca_solve(objflux[*,indx], objivar[*,indx], objloglam[*,indx], $
       zfit[indx], wavemin=wavemin, wavemax=wavemax, $
       niter=niter, nkeep=nkeep, newloglam=newloglam, eigenval=eigenval)

      spappend, newloglam, pcaflux, fullloglam, fullflux

   endfor

   ;----------
   ; Construct header for output file

   sxaddpar, hdr, 'OBJECT', 'STAR'
   sxaddpar, hdr, 'COEFF0', fullloglam[0]
   sxaddpar, hdr, 'COEFF1', objdloglam
   ; Add a space to the name below, so that 'F' appears as a string and
   ; not as a logical.
   for iclass=0, n_elements(classlist)-1 do $
    sxaddpar, hdr, 'NAME'+strtrim(string(iclass),2), classlist[iclass]+' '

   ;----------
   ; Write output file

   mwrfits, float(fullflux), outfile, hdr, /create

   return
end
;------------------------------------------------------------------------------
