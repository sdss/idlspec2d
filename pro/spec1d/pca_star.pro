; Generate both an output FITS file and a PostScript plot.
;------------------------------------------------------------------------------
pro pca_star, filename

   wavemin = 0
   wavemax = 0
   snmax = 100
   niter = 10
   cspeed = 2.99792458d5

   get_juldate, jd
   mjdstr = string(long(jd-2400000L), format='(i5)')
   outfile = 'spEigenStar-' + mjdstr + '.fits'
   plotfile = 'spEigenStar-' + mjdstr + '.ps'

   dfpsplot, plotfile, /color, /landscape
   colorvec = ['default', 'red', 'green', 'blue']

   ;----------
   ; Read the input spectra

   if (NOT keyword_set(filename)) then $
    filename = filepath('eigeninput_star.par', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   thisfile = (findfile(filename))[0]
   if (NOT keyword_set(thisfile)) then begin
      print, 'File not found: ' + filename
      return
   endif
   yanny_read, thisfile, pdat
   slist = *pdat[0]
   yanny_free, pdat

   readspec, slist.plate, slist.fiberid, mjd=slist.mjd, $
    flux=objflux, invvar=objivar, $
    andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=objloglam, /align

   ;----------
   ; Insist that all of the requested spectra exist

   imissing = where(plugmap.fiberid EQ 0, nmissing)
   if (nmissing GT 0) then begin
      for i=0, nmissing-1 do $
       print, 'Missing plate=', slist[imissing[i]].plate, $
        ' mjd=', slist[imissing[i]].mjd, $
        ' fiber=', slist[imissing[i]].fiberid
      message, string(nmissing) + ' missing object(s)'
   endif

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-sub residuals.

   objivar = skymask(objivar, andmask, ormask)
andmask = 0 ; Free memory
ormask = 0 ; Free memory

   nobj = (size(objflux, /dimens))[1]
   objdloglam = objloglam[1] - objloglam[0]

   if (keyword_set(snmax)) then begin
      ifix = where(objflux^2 * objivar GT snmax^2)
      if (ifix[0] NE -1) then objivar[ifix] = (snmax/objflux[ifix])^2
   endif

   ;----------
   ; Find the list of unique star types

   isort = sort(slist.class)
   classlist = slist[isort[uniq(slist[isort].class)]].class

   ;----------
   ; LOOP OVER EACH STAR TYPE

   for iclass=0, n_elements(classlist)-1 do begin

      ;----------
      ; Find the subclasses for this stellar type

      indx = where(slist.class EQ classlist[iclass], nindx)
      thesesubclass = slist[indx].subclass
      isort = sort(thesesubclass)
      subclasslist = thesesubclass[isort[uniq(thesesubclass[isort])]]
      nsubclass = n_elements(subclasslist)

      ;----------
      ; Solve for 2 eigencomponents if we have specified subclasses
      ; for this stellar type.

      if (nsubclass EQ 1) then nkeep = 1 $
       else nkeep = 2
      newloglam = objloglam
      pcaflux = pca_solve(objflux[*,indx], objivar[*,indx], objloglam, $
       slist[indx].cz/cspeed, wavemin=wavemin, wavemax=wavemax, $
       niter=niter, nkeep=nkeep, newloglam=newloglam, $
       eigenval=eigenval, acoeff=acoeff, usemask=usemask)

      ;----------
      ; Interpolate over bad flux values in the middle of a spectrum,
      ; and set fluxes to zero at the blue+red ends of the spectrum.

;      minuse = 1 ; ?
      minuse = floor((nindx+1) / 3.)
      qbad = usemask LT minuse

      ; Interpolate over all bad pixels
      for j=0, nkeep-1 do $
       pcaflux[*,j] = djs_maskinterp(pcaflux[*,j], qbad EQ 1, /const)

      ; Set bad pixels at the very start or end of the spectrum to zero instead
      npix = n_elements(qbad)
      if (qbad[0]) then $
       pcaflux[0:(where(qbad EQ 0))[0]-1,*] = 0
      if (qbad[npix-1]) then $
       pcaflux[(reverse(where(qbad EQ 0)))[0]+1:npix-1,*] = 0

; The following would plot the 0th object and overplot the best-fit PCA
;ii=0
;splot,10^newloglam,objflux[*,indx[ii]]
;junk=pcaflux[*,0] * (acoeff[0,ii])[0] + pcaflux[*,1] * (acoeff[1,ii])[0]
;soplot,10^newloglam,junk,color='red'

      ;----------
      ; Re-normalize the first eigenspectrum to a mean of 1
      norm = mean(pcaflux[*,0])
      pcaflux = pcaflux / norm
      acoeff = acoeff * norm

      ;----------
      ; Now loop through each stellar subclass and reconstruct
      ; an eigenspectrum for that subclass

      thesesubclassnum = lonarr(n_elements(thesesubclass))
      for isub=0, nsubclass-1 do begin
         ii = where(thesesubclass EQ subclasslist[isub])
         thesesubclassnum[ii] = isub
         if (nkeep EQ 1) then begin
            thisflux = pcaflux
         endif else begin
            aratio = acoeff[1,ii] / acoeff[0,ii]
            thisratio = median(aratio, /even)
            thisflux = pcaflux[*,0] + thisratio * pcaflux[*,1]
         endelse

         ; The output wavelength mapping is the same for everything,
         ; so we can simply stack the PCA spectra.
         if (NOT keyword_set(fullflux)) then fullflux = thisflux $
          else fullflux = [[fullflux], [thisflux]]

         if (NOT keyword_set(namearr)) then namearr = subclasslist[isub] $
          else namearr = [namearr, subclasslist[isub]]

         plotflux = thisflux / max(thisflux) ; Re-scale for plotting
         if (isub EQ 0) then $
          djs_plot, 10^newloglam, plotflux, color=colorvec[0], $
           xtitle='Wavelength [Ang]', ytitle='Flux [arbitrary units]', $
           title='STAR '+classlist[iclass]+': Eigenspectra Reconstructions' $
         else $
          djs_oplot, 10^newloglam, plotflux, $
           color=colorvec[isub MOD n_elements(colorvec)]
         nnew = n_elements(newloglam)
         xyouts, 10^newloglam[nnew-1], plotflux[nnew-1], $
          subclasslist[isub], align=-0.5, $
          color=djs_icolor(colorvec[isub MOD n_elements(colorvec)])
      endfor

      if (nkeep GT 1) then begin
         allratio = transpose(acoeff[1,*] / acoeff[0,*])
         isort = sort(thesesubclassnum)
         djs_plot, thesesubclassnum[isort], allratio[isort], ps=-4, $
          xrange=[-1,nsubclass], xstyle=1, xtickname=subclasslist, $
          xtickv=lindgen(nsubclass), xticks=nsubclass-1, $
          xtitle='Subclass', ytitle='Eigenvalue Ratio (a_1/a_0)', $
          title='STAR '+classlist[iclass]+': Eigenvalue Ratios'
         for j=0, n_elements(indx)-1 do $
          xyouts, thesesubclassnum[isort[j]], allratio[isort[j]], $
           align=0.0, orient=45, $
           string(slist[indx[isort[j]]].plate, slist[indx[isort[j]]].fiberid, $
           format='(i4,"-",i3)')
      endif

   endfor

   ;----------
   ; Construct header for output file

   sxaddpar, hdr, 'OBJECT', 'STAR'
   sxaddpar, hdr, 'COEFF0', newloglam[0]
   sxaddpar, hdr, 'COEFF1', objdloglam
   ; Add a space to the name below, so that 'F' appears as a string and
   ; not as a logical.
   for i=0, n_elements(namearr)-1 do $
    sxaddpar, hdr, 'NAME'+strtrim(string(i),2), namearr[i]+' '

   ;----------
   ; Write output file

   mwrfits, float(fullflux), outfile, hdr, /create

   dfpsclose

   return
end
;------------------------------------------------------------------------------
