
;------------------------------------------------------------------------------
; Find isolated fibers and return a boxcar extraction + optimal extraction
; Use this only on sparse-fiber flats, such as in MJD=51441
; Ex:
;   extract2, 'sdR-01-00001409.fit', fbox1, flux1, indir='51441'
;   extract2, 'sdR-01-00001415.fit', fbox2, flux2, indir='51441'
;   extract2, 'sdR-01-00001417.fit', fbox3, flux3, indir='51441'

pro extract2, flatname, fbox, flux, indir=indir, xcen=xcen, ycen=ycen

   sdssproc, flatname, flatimg, flativar, indir=indir

   xcen = trace320crude(flatimg, yset=ycen, maxdev=0.15)
   xy2traceset, ycen, xcen, tset, ncoeff=5, maxdev=0.1
   traceset2xy, tset, ycen, xsol

   fbox = extract_boxcar(flatimg, xsol)

   ; Select isolated fibers that are lit
   illum = total(fbox,1) / 2048. GT 1.e4 ; Illuminated or not
   isol = where(illum AND shift(illum,-1) EQ 0 AND shift(illum,1) EQ 0)

   ; Re-extract with wide boxcar
   fbox = extract_boxcar(flatimg, xsol, radius=5.0)

   extract_image, flatimg, flativar, xsol, 1.0, flux, fluxivar, $
    proftype=1, wfixed=[1], highrej=20, lowrej=25, nPoly=6, relative=1

   ; Trim to isolated fibers if there are any
   if (isol[0] EQ -1) then isol = indgen(320)
   fbox = fbox[*,isol]
   flux = flux[*,isol]
   xcen = xsol[*,isol]
   ycen = ycen[*,isol]

   return
end

