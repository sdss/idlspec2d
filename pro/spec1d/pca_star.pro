; Generate both an output FITS file and a PostScript plot.
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

   get_juldate, jd
   mjdstr = string(long(jd-2400000L), format='(i5)')
   outfile = 'spEigenStar-' + mjdstr + '.fits'
   plotfile = 'spEigenStar-' + mjdstr + '.ps'

   dfpsplot, plotfile, /color
   colorvec = ['default', 'red', 'green', 'blue']

   ;----------
   ; Read the input spectra

   filename = filepath('eigeninput_star.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   yanny_read, filename, pdat
   slist = *pdat[0]
   yanny_free, pdat

   readspec, slist.plate, slist.fiberid, mjd=slist.mjd, $
    flux=objflux, invvar=objivar, $
    andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=objloglam

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

      indx = where(slist.class EQ classlist[iclass])
      thesesubclass = slist[indx].subclass
      isort = sort(thesesubclass)
      subclasslist = thesesubclass[isort[uniq(thesesubclass[isort])]]
      nsubclass = n_elements(subclasslist)

      ;----------
      ; Solve for 2 eigencomponents if we have specified subclasses
      ; for this stellar type.

      if (nsubclass EQ 1) then nkeep = 1 $
       else nkeep = 2
      pcaflux = pca_solve(objflux[*,indx], objivar[*,indx], objloglam[*,indx], $
       slist[indx].z, wavemin=wavemin, wavemax=wavemax, $
       niter=niter, nkeep=nkeep, newloglam=newloglam, $
       eigenval=eigenval, acoeff=acoeff)

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
         spappend, newloglam, thisflux, fullloglam, fullflux
         if (NOT keyword_set(namearr)) then namearr = subclasslist[isub] $
          else namearr = [namearr, subclasslist[isub]]

         if (isub EQ 0) then $
          djs_plot, 10^newloglam, thisflux, color=colorvec[0], $
           xtitle='Wavelength [Ang]', ytitle='Flux [arbitrary units]', $
           title='STAR '+classlist[iclass]+': Eigenspectra Reconstructions' $
         else $
          djs_oplot, 10^newloglam, thisflux, $
           color=colorvec[isub MOD n_elements(colorvec)]
         nnew = n_elements(newloglam)
         xyouts, 10^newloglam[nnew-1], thisflux[nnew-1], $
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
   sxaddpar, hdr, 'COEFF0', fullloglam[0]
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
