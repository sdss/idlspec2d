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

   ;----------
   ; Do the same S/N calculation as in apo2d/quickextract.pro

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

   ;----------
   ; Spectra are already in 10^-17 flambda
   ; That's why we add 2.5*17 to the magnitude
  
   waveimg = 10^(loglam) 
   flambda2fnu = (waveimg*waveimg / 2.99792e18) # replicate(1,nfiber)

   synflux = transpose(filter_thru(finalflux*flambda2fnu, waveimg=waveimg, $
    mask=(finalivar LE 0)))

   synthmag = fltarr(5,nfiber)
   igood = where(synflux GT 0, ngood)
   if (ngood GT 0) then $
    synthmag[igood] = -2.5 * alog10(synflux[igood]) - 48.6 + 2.5*17.0

   ;----------
   ; Make S/N plot

   if (keyword_set(hdr)) then $
    plottitle = 'PLATE=' + strtrim(string(sxpar(hdr,'PLATEID')),2) $
     + '  MJD=' + strtrim(string(sxpar(hdr,'MJD')),2)
   plotsn, snvec, finalplugmap, plotfile=plotfile, plottitle=plottitle, $
    synthmag=synthmag, snplate=snplate

   ;----------
   ; Add the keywords SPEC1_G,SPEC1_R,... to the header.
; Add median offset and RMS keywords to header ???
; Move all these sxaddpar statements out of this routine?

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
   endif

   return
end
;------------------------------------------------------------------------------
