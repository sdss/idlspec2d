;+
; NAME:
;   pca_lrgtest
;
; PURPOSE:
;   Build PCA templates for LRGs within a specified redshift range.
;
; CALLING SEQUENCE:
;   pca_lrgtest, [ platenums, nkeep=, zrange= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   platenums  - Plate number(s) from which to select the LRGs; if not set,
;                then select all 'main' survey plates with QSURVEY=1
;                in the plate list file (which are required to be unique
;                tiles with good quality observations)
;   nkeep      - Number of PCA templates to keep (and to use for
;                noisy data replacement); default to 2.
;   zrange     - 2-element array with redshift range for fitting; if not set,
;                then call this routine iteratively with redshift ranges
;                starting at [0,0.05] and extending to [0.45,0.50], spaced
;                every 0.05 in redshift.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The output are written to FITS files named "spLRG_xxx_yyy.fits" where
;   "xxx" is the starting redshift times 100, and "yyy" is the ending
;   redshift times ten.
;
; EXAMPLES:
;   Create a set of PCA templates for LRGs from plates 400 through 409:
;     IDL> pca_lrgtest, 400+lindgen(10)
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;   djs_median
;   mwrfits
;   pca_solve()
;   platelist
;   readspec
;   splog
;   skymask()
;   sxaddhist
;   sxaddpar
;   wavevector()
;
; REVISION HISTORY:
;   23-Mar-2001  Written by David Schlegel, Princeton.
;   17-Sep-2003  Modified for N. Padmanabhan
;-
;------------------------------------------------------------------------------
pro pca_lrgtest, platenums, nkeep=nkeep, zrange=zrange

   if (NOT keyword_set(zrange)) then begin
      for z1=0.0, 0.46, 0.05 do begin
         pca_lrgtest, platenums, nkeep=nkeep, zrange=[z1,z1+0.05]
      endfor
      return
   endif

   ;----------
   ; Set defaults

   wavemin = 2500.
   wavemax = 9200.
   snmax = 100
   niter = 10
   if (NOT keyword_set(nkeep)) then nkeep = 2
   minuse = 5
   outfile = string(long(zrange[0]*100), long(zrange[1]*100), $
    format='("spLRG_", i3.3, "_", i3.3, ".fits")')

   ;----------
   ; Get the list of plates if not specified

   if (NOT keyword_set(platenums)) then begin
      platelist, plist=plist
      platenums = plist[ where(plist.qsurvey) ].plate
   endif

   ;----------
   ; Read the input spectra.
   ; Select good spectra of galaxies that were actually targetted
   ; as LRGs in the main survey.

   readspec, platenums, zans=zans, plug=plug
   indx = where(((plug.primtarget AND 2LL^5) NE 0 $
           OR (plug.primtarget AND 2LL^26) NE 0) $
    AND strtrim(zans.class) EQ 'GALAXY' $
    AND zans.zwarning EQ 0 $
    AND zans.z GE zrange[0] AND zans.z LE zrange[1])
   zans = zans[indx]
   splog, 'Number of objects = ', n_elements(zans)

   readspec, zans.plate, zans.fiberid, mjd=zans.mjd, $
    flux=objflux, invvar=objivar, $
    andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=objloglam

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
   ; Set the new wavelength mapping here...

   newloglam = wavevector(alog10(wavemin), alog10(wavemax), binsz=objdloglam)

   ;----------
   ; Do PCA solution

   pcaflux = pca_solve(objflux, objivar, objloglam, zans.z, $
    niter=niter, nkeep=nkeep, newloglam=newloglam, eigenval=eigenval, $
    usemask=usemask)
   pcaflux = float(pcaflux)

   ;----------
   ; Fill in bad data with a running median of good data

   qgood = usemask GE minuse
   igood = where(qgood, ngood)
   ibad = where(qgood EQ 0, nbad)
   medflux = 0 * pcaflux
   if (nbad GT 0) then begin
      for i=0, nkeep-1 do begin
         medflux[igood,i] = $
          djs_median(pcaflux[igood,i], width=51, boundary='nearest')
         medflux[*,i] = djs_maskinterp(medflux[*,i], qgood EQ 0, /const)
      endfor
      pcaflux[ibad,*] = medflux[ibad,*]
   endif

   ;----------
   ; Write output file

   sxaddpar, hdr, 'OBJECT', $
    'LRG z=' + string(zrange[0]) + ' ' + string(zrange[1])
   sxaddpar, hdr, 'COEFF0', newloglam[0]
   sxaddpar, hdr, 'COEFF1', objdloglam
   sxaddpar, hdr, 'NPLATE', n_elements(platenums), ' Number of plates used'
   sxaddpar, hdr, 'NGALAXY', n_elements(zans), ' Number of galaxies in fit'
   for i=0, n_elements(eigenval)-1 do $
    sxaddpar, hdr, 'EIGEN'+strtrim(string(i),1), eigenval[i]
   for i=0, n_elements(platenums)-1 do $
    sxaddhist, 'Plate ' + strtrim(platenums[i],2) + ' with ' $
     + strtrim(long(total(zans.plate EQ platenums[i])),2) + ' galaxies', hdr

   mwrfits, pcaflux, outfile, hdr, /create

   return
end
;------------------------------------------------------------------------------
