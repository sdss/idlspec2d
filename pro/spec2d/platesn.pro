;+
; NAME:
;   spcoadd_frames
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
;   plotfile       - If set, then write PostScript plot to this file.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   hdr            - [Modified]
;   snvec          - S/N vector for g,r,i bands
;   synthmag       - Synthetic magnitudes from convolution with fiducial
;                      Filter curves 
;
; COMMENTS:
;   Median fluxes are used in each band-pass to generate the synthesized
;   magnitudes.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_median()
;   filter_thru()
;   pixelmask_bits()
;   plotsn
;   splog
;   sxaddpar
;
; REVISION HISTORY:
;   06-Oct-2000  Written by S. Burles & D. Schlegel
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
;  Horribly bulky, but what can we do?

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
    mask=(finalivar LE 0), /norm))

   synthmag = fltarr(3,nfiber)
   posfilter = where(filter[1:3,*] GT 0)
   if posfilter[0] NE -1 then $
    synthmag[posfilter] = $
     -2.5 * alog10((filter[1:3,*])[posfilter]) - 48.6 + 2.5*17.0

   ;----------
   ; Make S/N plot

   plotsn, snvec, finalplugmap, plotfile=plotfile, plottitle=plottitle, $
    synthmag=synthmag, snplate=snplate

   ;----------
   ; Bad Fiber roll call
  
   rollcall = lonarr(3,8)

   for i=0, 7 do $
    rollcall[*,i] = $
     [ total(total((finalandmask[gwave,*] AND 2L^i) GT 0, 1) GT 0), $
       total(total((finalandmask[rwave,*] AND 2L^i) GT 0, 1) GT 0), $
       total(total((finalandmask[iwave,*] AND 2L^i) GT 0, 1) GT 0)  ]
   
   splog, "Fiber Problem",  "  # in g'",  "  # in r'", "  # in i'" , $
          format='(3x,a20, 3(a9))'

   ; The following sets up the names of the pixelmask bits...
   junk = pixelmask_bits('NODATA')

   for i=0, 7 do begin
      label = maskbits[where(maskbits.bit EQ i)].label
      info = string(format='(a20,3(i9.3))',label,rollcall[*,i])
      splog, label, rollcall[*,i], format='(a20,3(i9))'
   endfor

   ;----------
   ; Modify the header

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
