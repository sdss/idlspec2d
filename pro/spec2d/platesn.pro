;+
; NAME:
;   platesn
;
; PURPOSE:
;   Generate S/N plots for an entire plate.
;
; CALLING SEQUENCE:
;   platesn, finalflux, finalivar, finalandmask, finalplugmap, loglam, $
;    [ hdr=, plotfile=, snvec=, synthmag= ]
;
; INPUTS:
;   finalflux      - 
;   finalivar      - 
;   finalandmask   - 
;   finalplugmap   - 
;   loglam         - 
;
; REQUIRED KEYWORDS:
;   fcalibprefix   - Prefix for flux-calibration files.
;
; OPTIONAL KEYWORDS:
;   hdr            - FITS header; if specified, then keywords are added.
;                    The PLATE and MJD for the plot title are from this header.
;   plotfile       - If set, then write PostScript plot to this file.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   hdr            - [Modified]
;   snvec          - S/N vector for g,r,i bands
;   synthmag       - Synthetic magnitudes from convolution with fiducial
;                    filter curves 
;
; COMMENTS:
;   Median fluxes are used in each band-pass to generate the synthesized
;   magnitudes.
;   
;   Header keywords ROFFSET & RSIGMA store the offset and standard
;   deviation of the r-band difference in the spectro and photo mags.  
;   GIOFF and GISIGMA store the offset and scatter in the (g-i) color.
;   
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_mean()
;   djs_median()
;   filter_thru()
;   pixelmask_bits()
;   plotsn
;   sdss_flagname()
;   splog
;   sxaddpar
;   sxpar()
;
; REVISION HISTORY:
;   06-Oct-2000  Written by S. Burles & D. Schlegel
;   20-Oct-2002  C. Tremonti added header keywords to access spectrophotometry
;-
;------------------------------------------------------------------------------
pro platesn, finalflux, finalivar, finalandmask, finalplugmap, loglam, $
 hdr=hdr, plotfile=plotfile, snvec=snvec, synthmag=synthmag, filtsz=filtsz

   if NOT keyword_set(filtsz) then filtsz=25

   common com_maskbits, maskbits

   dims = size(finalflux, /dimens)
   npix = dims[0]
   nfiber = dims[1]
   
   gwave = where(loglam GT alog10(4000) AND loglam LT alog10(5500))
   rwave = where(loglam GT alog10(5600) AND loglam LT alog10(6900))
   iwave = where(loglam GT alog10(6910) AND loglam LT alog10(8500))

   snimg = finalflux * sqrt(finalivar)
   snvec = fltarr(3, nfiber)

;  Do the same S/N calculation as in apo2d/quickextract.pro
;  Horribly bulky, but what can we do???

   for ifib=0, nfiber-1 do begin
     sntemp = 0.0
     ig = where(finalivar[gwave,ifib] GT 0, nwave)
     if (nwave GT filtsz) then $
       sntemp = djs_median(snimg[gwave[ig],ifib], $
                           width=filtsz, boundary='reflect')
     sng = djs_mean(sntemp)

     sntemp = 0.0
     ig = where(finalivar[rwave,ifib] GT 0, nwave)
     if (nwave GT filtsz) then $
       sntemp = djs_median(snimg[rwave[ig],ifib], $
                           width=filtsz, boundary='reflect')
     snr = djs_mean(sntemp)

     sntemp = 0.0
     ig = where(finalivar[iwave,ifib] GT 0, nwave)
     if (nwave GT filtsz) then $
       sntemp = djs_median(snimg[iwave[ig],ifib], $
                           width=filtsz, boundary='reflect')
     sni = djs_mean(sntemp)

     snvec[*,ifib] = [sng, snr, sni]
   endfor

   ;--------------------------------------------------------------------
   ; Spectra are already in 10^-17 flambda
   ; That's why we add 2.5*17 to the magnitude
  
   waveimg = 10^(loglam) 
   flambda2fnu = (waveimg*waveimg / 2.99792e18) # replicate(1,nfiber)

   filter = transpose(filter_thru(finalflux*flambda2fnu, waveimg=waveimg, $
    mask=(finalivar LE 0)))

   synthmag = fltarr(3,nfiber)
   posfilter = where(filter[1:3,*] GT 0)
   if posfilter[0] NE -1 then $
    synthmag[posfilter] = $
     -2.5 * alog10((filter[1:3,*])[posfilter]) - 48.6 + 2.5*17.0

   ;----------
   ; Make S/N plot

   plottitle = 'PLATE=' + strtrim(string(sxpar(hdr,'PLATEID')),2) $
    + '  MJD=' + strtrim(string(sxpar(hdr,'MJD')),2)
   plotsn, snvec, finalplugmap, plotfile=plotfile, plottitle=plottitle, $
    synthmag=synthmag, snplate=snplate, roffset = roffset, rsigma = rsigma, $
    gioffset = gioffset, gisigma = gisigma

   ;----------
   ; Print roll call of bad fibers and bad pixels.
   ; Assume that all mask bits under bit #16 are for the entire fiber.

   bitlabel = sdss_flagname('SPPIXMASK', 2UL^32-1, /silent)
   bitnum = where(bitlabel NE '', nlabel)
   bitlabel = bitlabel[bitnum]

   ncam = 4
   camnames = ['b1','r1','b2','r2']
   camwave1 = [3500, 6000, 3500, 6000]
   camwave2 = [6000, 9500, 6000, 9500]
   specidimg = bytarr(npix,nfiber) + 1B
   specidimg[*,nfiber/2:nfiber-1] = 2

   rollcall = fltarr(ncam)

   splog, ' '
   splog, camnames, format='(24x,4a7)'
   splog, format='(24x,4(" ------"))'
   for ilabel=0, nlabel-1 do begin
      for icam=0, ncam-1 do begin
         specid = fix( strmid(camnames[icam],1) )
         thiscam = strmid(camnames[icam],0,1)
         fib1 = (specid-1) * (nfiber/2)
         fib2 = fib1 + (nfiber/2) - 1

         if (bitnum[ilabel] LT 16) then begin
            ; CASE: Fiber mask
            ; Look at only the g-band or i-band wavelengths, so that
            ; we ignore any overlap in wavelength between cameras.
            if (thiscam EQ 'b') then indx = gwave $
             else indx = iwave
            qmask = (finalandmask[indx,fib1:fib2] AND 2L^bitnum[ilabel]) NE 0
            rollcall[icam] = total( total(qmask,1) NE 0 )
         endif else begin
            ; CASE: Pixel mask
            ; Include overlap in wavelength between cameras, calling
            ; any wavelengths < 6000 Ang the blue camera, and
            ; any wavelengths > 6000 Ang the red camera.
            indx = where(loglam GT alog10(camwave1[icam]) $
             AND loglam LT alog10(camwave2[icam]), nindx)
            qmask = (finalandmask[indx,fib1:fib2] AND 2L^bitnum[ilabel]) NE 0
            rollcall[icam] = 100.0 * total(qmask NE 0) $
             / (nindx * (fib2-fib1+1))
         endelse
      endfor
      if (bitnum[ilabel] LT 16) then begin
         splog, 'N(fiber) '+bitlabel[ilabel]+'               ', $
          long(rollcall), format='(a24,4i7)'
      endif else begin
         splog, '%(pixel) '+bitlabel[ilabel]+'               ', $
          rollcall, format='(a24,4f7.2)'
      endelse
   endfor

   ;----------
   ; Add the keywords SPEC1_G,SPEC1_R,... to the header.

   if (keyword_set(hdr)) then begin
      bands = ['G','R','I']
      snmag = [20.2, 20.25, 19.9]

      for ispec=1, 2 do begin
         for bb=0, n_elements(bands)-1 do begin
            key1 = 'SPEC'+ strtrim(ispec,2)+'_'+bands[bb]
            comment = string(format='(a,i2,a,f5.2)', $
             '(S/N)^2 for spec ', ispec, ' at mag ', snmag[bb])
            sxaddpar, hdr, key1, snplate[ispec-1,bb], comment, before='LOWREJ'
         endfor
      endfor
 
      ;-----------------
      ; Add keywords describing the agreement of the syntmags and the
      ; photo mags in the plugmap.  

      sxaddpar, hdr, 'ROFFSET1', roffset[0], $
                'Mean r-band mag difference (spectro mag - photo mag)'
      sxaddpar, hdr, 'RSIGMA1', rsigma[0], $
             'Stddev of r-band mag difference (spectro mag - photo mag)'
      sxaddpar, hdr, 'GIOFF1', gioffset[0], $
                'Mean (g - i) color difference (spectro mag - photo mag)'
      sxaddpar, hdr, 'GISIGMA1', gisigma[0], $
                'Stddev of (g - i) color difference (spectro mag - photo mag)'
      sxaddpar, hdr, 'ROFFSET2', roffset[1], $
                'Mean r-band mag difference (spectro mag - photo mag)'
      sxaddpar, hdr, 'RSIGMA2', rsigma[1], $
             'Stddev of r-band mag difference (spectro mag - photo mag)'
      sxaddpar, hdr, 'GIOFF2', gioffset[1], $
                'Mean (g - i) color difference (spectro mag - photo mag)'
      sxaddpar, hdr, 'GISIGMA2', gisigma[1], $
                'Stddev of (g - i) color difference (spectro mag - photo mag)'

   endif

   return
end
;------------------------------------------------------------------------------
