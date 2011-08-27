;+
; NAME:
;   platesn
;
; PURPOSE:
;   Generate S/N plots for an entire plate.
;
; CALLING SEQUENCE:
;   platesn, [ objflux, objivar, andmask, plugmap, loglam, $
;    filtsz=, hdr=, platefile=, plotfile=, snvec=, synthmag= ]
;
; INPUTS:
;
; OPTIONAL KEYWORDS:
;   objflux    - 
;   objivar    - 
;   andmask    - 
;   plugmap    - 
;   loglam     - 
;   filtsz     - Filter size for median-filtering the S/N values within
;                a spectrum; default to 25 pix.
;   hdr        - FITS header; if specified, then keywords are added.
;                The PLATE and MJD for the plot title are from this header.
;   platefile  - If set, then read OBJFLUX and all the other inputs from
;                this spPlate file instead of using those inputs;
;                also, generate PostScript plots using the PLATESN procedure.
;   plotfile   - If set, then write PostScript plot to this file.
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
;   djs_iterstat
;   djs_mean()
;   djs_median()
;   filter_thru()
;   pixelmask_bits()
;   plotsn
;   sdss_flagval()
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
pro platesn, objflux, objivar, andmask, plugmap, loglam, $
 hdr=hdr, platefile=platefile, plotfile=plotfile, $
 snvec=snvec, synthmag=synthmag, filtsz=filtsz

   common com_maskbits, maskbits

   if NOT keyword_set(filtsz) then filtsz=25

   if (keyword_set(platefile)) then begin
      objflux = mrdfits(platefile,0,hdr)
      objivar = mrdfits(platefile,1)
      andmask = mrdfits(platefile,2)
      plugmap = mrdfits(platefile,5)
      loglam = sxpar(hdr, 'COEFF0') $
       + dindgen(sxpar(hdr, 'NAXIS1')) * sxpar(hdr, 'COEFF1')
   endif

   dims = size(objflux, /dimens)
   npix = dims[0]
   nfiber = dims[1]
   
   gwave = where(loglam GT alog10(4000) AND loglam LT alog10(5500))
   rwave = where(loglam GT alog10(5600) AND loglam LT alog10(6900))
   iwave = where(loglam GT alog10(6910) AND loglam LT alog10(8500))

   snimg = objflux * sqrt(objivar)
   snvec = fltarr(3, nfiber)

   ;----------
   ; Do the same S/N calculation as in apo2d/quickextract.pro

   for ifib=0, nfiber-1 do begin
      sntemp = 0.0
      ig = where(objivar[gwave,ifib] GT 0, nwave)
      if (nwave GT filtsz) then $
       sntemp = djs_median(snimg[gwave[ig],ifib], $
        width=filtsz, boundary='reflect')
      sng = djs_mean(sntemp)

      sntemp = 0.0
      ig = where(objivar[rwave,ifib] GT 0, nwave)
      if (nwave GT filtsz) then $
       sntemp = djs_median(snimg[rwave[ig],ifib], $
        width=filtsz, boundary='reflect')
      snr = djs_mean(sntemp)

      sntemp = 0.0
      ig = where(objivar[iwave,ifib] GT 0, nwave)
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

   synflux = transpose(filter_thru(objflux*flambda2fnu, waveimg=waveimg, $
    mask=(objivar LE 0)))

   synthmag = fltarr(5,nfiber)
   igood = where(synflux GT 0, ngood)
   if (ngood GT 0) then $
    synthmag[igood] = -2.5 * alog10(synflux[igood]) - 48.6 + 2.5*17.0

   ;----------
   ; Make S/N plot

   if (keyword_set(hdr)) then $
    plottitle = 'PLATE=' + strtrim(string(sxpar(hdr,'PLATEID')),2) $
     + '  MJD=' + strtrim(string(sxpar(hdr,'MJD')),2)
   filter = ['g','r','i']
   plotsn, snvec, plugmap, plotfile=plotfile, plottitle=plottitle, $
    sncode='spcombine', filter=filter, synthmag=synthmag, $
    snplate=snplate, dered_snplate=dered_snplate, specsnlimit=specsnlimit, $
    redden=sxpar(hdr,'REDDEN*')

   ;----------
   ; Add header keywords if HDR is passed.

   if (keyword_set(hdr)) then begin
      ;----------
      ; Add the keywords SPEC1_G,SPEC1_R,... with the (S/N)^2 values
      ; per spectrograph.
      ; Also add keywords SN2EXT1G, SN2EXT1R, ... with extinction
      ; corrected (S/N)^2 values

      bands = ['G','R','I']

      for ispec=1, 2 do begin
         for bb=0, n_elements(bands)-1 do begin
            ; Standard (S/N)^2
            key1 = 'SPEC'+ strtrim(ispec,2)+'_'+strupcase(filter[bb])
            comment = string(format='(a,i2,a,f5.2)', $
             ' (S/N)^2 for spec ', ispec, ' at mag ', specsnlimit[bb].snmag)
            sxaddpar, hdr, key1, snplate[ispec-1,bb], comment, before='LOWREJ'
            
            ; Extinction corrected (S/N)^2
            key1 = 'SN2EXT'+strtrim(ispec,2)+strupcase(filter[bb])
            comment = ' Extinction corrected (S/N)^2'
            sxaddpar, hdr, key1, dered_snplate[ispec-1,bb], $
             comment, before='LOWREJ'
         endfor
      endfor

      ;----------
      ; Add the keywords GOFFSTD,ROFFSTD,... and GRMSSTD,RRMSSTD,...
      ; with the spectro-photometric offsets (OFF) + dispersions (RMS) in the
      ; three bands G,R,I for standard stars (STD) and "MAIN" galaxies (GAL).
      ; This is computed for objects on both spectrographs.

      fname = ['u','g','r','i','z']
      mlimit = 22.5 ; Mag limit in gri bands for both photo + spectro-photo mags
      if (tag_exist(plugmap,'CALIBFLUX')) then begin
         minflux = 0.1
         thismag = (22.5 - 2.5*alog10(plugmap.calibflux > minflux)) $
          * (plugmap.calibflux GT minflux)
         if (total(plugmap.calibflux_ivar GT 0) GT 0) then $
          thismag = thismag * (plugmap.calibflux_ivar GT 0)
      endif
      qgood = thismag[1,*] GT 0 AND thismag[1,*] LT mlimit $
       AND thismag[2,*] GT 0 AND thismag[2,*] LT mlimit $
       AND thismag[3,*] GT 0 AND thismag[3,*] LT mlimit $
       AND synthmag[1,*] GT 0 AND synthmag[1,*] LT mlimit $
       AND synthmag[2,*] GT 0 AND synthmag[2,*] LT mlimit $
       AND synthmag[3,*] GT 0 AND synthmag[3,*] LT mlimit
      qstd = strtrim(plugmap.objtype,2) EQ 'SPECTROPHOTO_STD' $
       OR strtrim(plugmap.objtype,2) EQ 'REDDEN_STD'
      qgal = (plugmap.primtarget AND sdss_flagval('TARGET','GALAXY')) NE 0
      istd = where(qstd AND qgood, nstd)
      igal = where(qgal AND qgood, ngal)

      for itype=0, 1 do begin
         if (itype EQ 0) then begin
            indx = istd
            tstring = 'STD'
            tname = 'std stars'
         endif else begin
            indx = igal
            tstring = 'GAL'
            tname = 'main galaxies'
         endelse
         if (indx[0] NE -1) then begin
            sxaddpar, hdr, 'N'+tstring, n_elements(indx), $
             ' Number of (good) '+tname, before='LOWREJ'
            for iband=1, 3 do begin
               djs_iterstat, synthmag[iband,indx] - thismag[iband,indx], $
                median=med1, sigma=sig1
               key1 = strupcase(fname[iband]) + 'OFF' + tstring
               comment = ' Spectrophoto offset for ' + tname + ' in ' $
                + strupcase(fname[iband])+'-band'
               sxaddpar, hdr, key1, med1, comment, before='LOWREJ'
               key1 = strupcase(fname[iband]) + 'RMS' + tstring
               comment = ' Spectrophoto RMS for ' + tname + ' in ' $
                + strupcase(fname[iband])+'-band'
               sxaddpar, hdr, key1, sig1, comment, before='LOWREJ'
            endfor
            for icolor=0, 1 do begin
               djs_iterstat, synthmag[1+icolor,indx] - thismag[1+icolor,indx] $
                - (synthmag[2+icolor,indx] - thismag[2+icolor,indx]), $
                median=med1, sigma=sig1
               colorname = strupcase(fname[icolor+1]+fname[icolor+2])
               key1 = strupcase(colorname) + 'OFF' + tstring
               comment = ' Spectrophoto offset for ' + tname + ' in (' $
                + colorname + ')'
               sxaddpar, hdr, key1, med1, comment, before='LOWREJ'
               key1 = strupcase(colorname) + 'RMS' + tstring
               comment = ' Spectrophoto RMS for ' + tname + ' in (' $
                + colorname + ')'
               sxaddpar, hdr, key1, sig1, comment, before='LOWREJ'
            endfor
         endif
      endfor
   endif

   return
end
;------------------------------------------------------------------------------
