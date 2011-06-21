; Generate both an output FITS file and a PostScript plot.
;
; Loosely based on pca_star by Schlegel (????),
; Substantially rewritten by Bolton (utah2011june) to inherit Indo-US archetype spectra,
; to do PCA of BOSS spectra for weird types not represented in Indo-US,
; and to simplify some of the pleasantries.
;------------------------------------------------------------------------------
pro pca_star_boss, filename

; Read in the BOSS-sampled indo-US archetype spectra
; and define the target baseline from them:
ius_fname = filepath('spArchIndoUS.fits', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')

influx = mrdfits(ius_fname,0,ius_hdr)
coeff0 = sxpar(ius_hdr, 'COEFF0')
coeff1 = sxpar(ius_hdr, 'COEFF1')
npix = sxpar(ius_hdr, 'NAXIS1')
newloglam = findgen(npix) * coeff1 + coeff0
n_ius = sxpar(ius_hdr, 'NAXIS2')

   snmax = 100
   niter = 10
   cspeed = 2.99792458d5

  minfrac = 0.75 ; minimum fraction of spectra contributing to a pixel
  pbuff = 25 ; pixel buffer for averaging for padding value

   get_juldate, jd
   mjdstr = string(long(jd-2400000L), format='(i5)')
   outfile = 'spEigenStar-' + mjdstr + '.fits'
   plotfile = 'spEigenStar-' + mjdstr + '.ps'

    ;----------
   ; Read the input spectra

   if (NOT keyword_set(filename)) then $
    filename = filepath('eigeninput_star_boss.dat', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   thisfile = (findfile(filename))[0]
   if (NOT keyword_set(thisfile)) then begin
      print, 'File not found: ' + filename
      return
   endif
   readcol, thisfile, plate, mjd, fiberid, zstar, subclass, $
    format='L,L,L,F,A', comment='#'

   readspec, plate, fiberid, mjd=mjd, plugmap=plugmap, $
    flux=flux, invvar=invvar, loglam=loglam, andmask=andmask, ormask=ormask, /align

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

   invvar = skymask(invvar, andmask, ormask)
   andmask = 0 ; Free memory
   ormask = 0 ; Free memory

   nobj = (size(flux, /dimens))[1]
   objdloglam = loglam[1] - loglam[0]

   if (keyword_set(snmax)) then begin
      ifix = where(flux^2 * invvar GT snmax^2)
      if (ifix[0] NE -1) then invvar[ifix] = (snmax/flux[ifix])^2
   endif

   ;----------
   ; Find the list of unique star types

   classlist = subclass[sort(subclass)]
   classlist = classlist[uniq(classlist)]

; How many classes are there?
nclass = n_elements(classlist)

; Array to hold the new templates:
tflux = fltarr(npix, nclass)

   ;----------
   ; LOOP OVER EACH STAR TYPE

   for iclass=0, n_elements(classlist)-1 do begin

      ;----------
      ; Find the examples of this stellar type

      indx = where(subclass EQ classlist[iclass], nindx)

; Do the PCA:
      pcaflux = pca_solve(flux[*,indx], invvar[*,indx], loglam, $
                          zstar[indx], niter=niter, usemask=usemask, $
                          newloglam=newloglam, nkeep=1, nreturn=1)

; Pad the cruddy bits by hand:
      minuse = round(minfrac * nindx)
      wh_use = where(usemask ge minuse)
      umax = max(wh_use)
      umin = min(wh_use)
      hi_mean = mean(pcaflux[umax-pbuff+1:umax])
      lo_mean = mean(pcaflux[umin:umin+pbuff-1])
      pcaflux[0:umin-1] = lo_mean
      pcaflux[umax+1:npix-1] = hi_mean

; Assign to output array:

      tflux[*,iclass] = pcaflux

   endfor

   ;----------
   ; Construct header for output file

hdr = ius_hdr
sxdelpar, hdr, 'NAXIS1'
sxdelpar, hdr, 'NAXIS2'
for i = 0L, nclass-1 do sxaddpar, hdr, 'NAME'+strtrim(string(i)+n_ius,2), classlist[i], $
 'BOSS PCA template'

; Concatenate output arrays:
outflux = [[influx], [tflux]]
n_out = (size(outflux))[2]

; Normalize to have average flux of unity in each spectrum:
for i = 0L, n_out-1 do outflux[*,i] = outflux[*,i] / mean(outflux[*,i])

   mwrfits, float(outflux), outfile, hdr, /create

   return
end
;------------------------------------------------------------------------------
