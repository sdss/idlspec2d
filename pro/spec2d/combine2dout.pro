;+
; NAME:
;   combine2dout
;
; PURPOSE:
;   Combine several reduced frames of the same objects
;
; CALLING SEQUENCE:
;   combine2dout, filenames, outputroot, spectrographid, $
;    binsz, zeropoint, nord=nord, wavemin=wavemin, $
;    bkptbin=bkptbin, everyn=everyn, display=display, window=window, $
;    individual=individual
;
; INPUTS:
;   filenames      - Name(s) of files written by SPREDUCE
;   outputroot     - Root name for all output files
;   spectrographid - Spectrograph ID (1 or 2) for use in computing fiber
;                    number for output file name; fibers from spectro-1
;                    are numbered from 1, those from spectro-2 from 321.
;
; REQUIRED KEYWORDS:
;
; OPTIONAL KEYWORDS:
;
;   display        - Show QA plot for each combined spectrum
;   individual     - Append individual spectra in HDU's after combined spectra
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;	This routine also outputs original 2048 spectra with mask pixels
;       replaced with their b-spline values.  Helpful for FFT
;
;
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_oplot
;   djs_plot
;   djs_maskinterp()
;   mrdfits()
;   writefits
;
; INTERNAL SUPPORT PROCEDURES:
;   makelabel()
;
; REVISION HISTORY:
;   ??-Sep-1999  Written by S. Burles
;   02-Jan-2000  Modified by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------

function makelabel, hdr

   plate = strtrim(string(sxpar(hdr, 'PLATEID')),2)
   camera = strtrim(sxpar(hdr, 'CAMERAS'),2)
   mjd = strtrim(string(sxpar(hdr, 'MJD')),2)
   seqid =  strtrim(string(sxpar(hdr, 'SEQID')),2)
   expos =  strtrim(string(sxpar(hdr, 'EXPOSURE')),2)

   return, plate+'-'+camera+'-'+mjd+'-'+seqid+'-'+expos
end

;------------------------------------------------------------------------------

pro combine2dout, filenames, outputroot, spectrographid, $
 binsz, zeropoint, nord=nord, wavemin=wavemin, $
 bkptbin=bkptbin, everyn=everyn, display=display, window=window, $
 individual=individual, maxsep=maxsep

   ; Initial binning was 69 km/s per pixel (with a sigma of 1.0 pixel)
   ; 69.02977415 km/s is log lambda 10^-4

   if (NOT keyword_set(binsz)) then binsz = 1.0d-4 $
   else binsz = double(binsz)

   if (NOT keyword_set(zeropoint)) then zeropoint = 3.5d
   if (NOT keyword_set(nord)) then nord = 3
   if (NOT keyword_set(bkptbin)) then bkptbin = binsz*1.2

   if (NOT keyword_set(maxsep)) then maxsep = 2.0*binsz

;------------------------------------------------------------------------
;  We should sort here to get blue together and then reds afterward
;

   filenames = filenames[sort(filenames)]

   nfiles = N_elements(filenames)
   if (nfiles EQ 0) then return

   redfiles = 0
   bluefiles = 0

   ;------------------------------------------------------------
   ; Somewhere up here we should add a flag combination if gives a
   ; positive intersection with fibermask... To reject single fibers.


   ; Read first file

   flux     = mrdfits(filenames[0], 0, hdr)
   fluxivar = mrdfits(filenames[0], 1)
   plugmap  = mrdfits(filenames[0], 2)
   wset     = mrdfits(filenames[0], 3)
   traceset2xy, wset, pixnorm, wave

   dims = size(flux, /dimens)
   npix  = dims[0]

   nfiber = dims[1]
   specnum = bytarr(npix)
   isred = (strpos(sxpar(hdr, 'CAMERAS'),'r') EQ 0)
   bluered = bytarr(npix) + isred

   ;-------------------------------------------------------------
   ; Try to read in pixelmask

   tt = mrdfits(filenames[0], 4)
   if (size(tt,/n_elem) EQ 1) then pixelmask = lonarr(npix,nfiber) $
   else pixelmask = tt

   ;-------------------------------------------------------------
   ; Try to read in fibermask

   tt = mrdfits(filenames[0], 5)
   if (size(tt,/n_elem) EQ 1) then fibermask = transpose(bytarr(nfiber)) $
   else fibermask = transpose(tt)


   if (isred) then redfiles = redfiles + 1 
   if (strpos(sxpar(hdr, 'CAMERAS'),'b') EQ 0) then bluefiles = bluefiles + 1 

   label =  makelabel(hdr)
   exptime = sxpar(hdr, 'EXPTIME')

   for i=1, nfiles-1 do begin
      tempflux = mrdfits(filenames[i], 0, hdr)
      tempivar = mrdfits(filenames[i], 1)
      tempplug = mrdfits(filenames[i], 2)
      tempwset = mrdfits(filenames[i], 3)

      traceset2xy, tempwset, pixnorm, tempwave

      npixhere = (size(tempflux, /dimens))[0]
      if (npixhere NE npix) then print, 'trouble npix NE npixhere'

      nfiberhere = (size(tempflux, /dimens))[1]

      if (nfiberhere NE nfiber) then begin
         splog, 'ABORT: Different files have different number of fibers??'
         return
      endif


      ;-------------------------------------------------------------
      ; Try to read in pixelmask

      tt = mrdfits(filenames[i], 4)
      if (size(tt,/n_elem) EQ 1) then pixelmask = [pixelmask, lonarr(npix,nfiber)] $
      else pixelmask = [pixelmask, tt]

      ;-------------------------------------------------------------
      ; Try to read in fibermask

      tt = mrdfits(filenames[i], 5)
      if (size(tt,/n_elem) EQ 1) then $
       fibermask = [fibermask, transpose(bytarr(nfiber))] $
      else $
       fibermask = [fibermask, transpose(tt)]

      if (keyword_set(window)) then $
       tempivar[0:window] = tempivar[0:window] * findgen(window+1) / window

      flux = [flux, tempflux]
      fluxivar = [fluxivar, tempivar]
      wave = [wave, tempwave]
      specnum = [specnum, bytarr(npix) + i]

      isred = (strpos(sxpar(hdr, 'CAMERAS'),'r') EQ 0)
      bluered = [bluered, bytarr(npix) + isred] 

      if (isred) then redfiles = redfiles + 1 
      if (strpos(sxpar(hdr, 'CAMERAS'),'b') EQ 0) then $
                 bluefiles = bluefiles + 1 
        
      label = [label, makelabel(hdr)]
      exptime = exptime + sxpar(hdr, 'EXPTIME')
   endfor


   splog, 'Found '+string(redfiles)+' Red files'
   splog, 'Found '+string(bluefiles)+' Blue files'
   if (redfiles LT 2 OR bluefiles LT 2) then begin
      splog, 'For the time being, I expect at least 2 of each red and blue to work'
      return
   endif

;set_plot,'x'
;ifib=237
;colorv=['default','red','green','blue','magenta','yellow','cyan','default']
;lam = 10^wave
;splot, lam[where(specnum EQ 0),ifib], $
; flux[where(specnum EQ 0),ifib], color=colorv[0]
;for jnum=1, nfiles-1 do $
; soplot, lam[where(specnum EQ jnum),ifib], $
;  flux[where(specnum EQ jnum),ifib], color=colorv[jnum]
;stop

   totalpix = (size(flux, /dimens))[0]

   redpix = where(bluered, numred)
   bluepix = where(bluered EQ 0, numblue)
   if (numblue GT 0 AND numred GT 0) then $
    exptime = exptime * 0.5

   ; Fix up new header, any one should do to start with

   ncoeff = sxpar(hdr, 'NWORDER')
   for i=2, ncoeff-1 do sxdelpar, hdr, 'COEFF'+strtrim(string(i),2)

   sxaddpar, hdr, 'PIXMIN', 0.000, 'Place holding'
   sxaddpar, hdr, 'PIXMAX', float(npix - 1), 'Place holding'
   sxaddpar, hdr, 'CREATORS', 'Burles & Schlegel (1999) IDLspec', after='SDSS'

   ; Now get rid of exposure, and add list of exposures

   sxdelpar, hdr, 'EXPOSURE'
   sxdelpar, hdr, 'SEQID'
   sxaddpar, hdr, 'NEXP', nfiles, $
    'Number of exposures in this file', after='TELESCOP'
   for i=0,nfiles-1 do $
    sxaddpar, hdr, 'EXPID'+strtrim(string(i),2), label[i], $
     'ID string for exposure '+strtrim(string(i),2), before='EXPTIME'

   sxaddpar, hdr, 'EXPTIME', exptime, 'total exposure time (seconds)'
   sxaddpar, hdr, 'COMBINE2', systime(), $
    'COMBINE2DSCOTT finished', after='EXPTIME'

   scale = fltarr(nfiber)
   blueflux = fltarr(nfiber)
   redflux = fltarr(nfiber)

   nonzero = where(fluxivar GT 0.0)
   minfullwave = min(wave[nonzero])
   maxfullwave = max(wave[nonzero])

   ; Get max and min from good pixels

   if (NOT keyword_set(wavemin)) then begin
      spotmin = fix((minfullwave - zeropoint)/binsz) + 1
      spotmax = fix((maxfullwave - zeropoint)/binsz)
      wavemin = spotmin * binsz + zeropoint
      wavemax = spotmax * binsz + zeropoint
   endif else begin
      spotmin = 0
      if (NOT keyword_set(wavemax)) then begin
        spotmax = fix((maxfullwave - wavemin)/binsz)
        wavemax = spotmax * binsz + wavemin
      endif else spotmax = fix((wavemax - wavemin)/binsz)
   endelse

   npix = spotmax - spotmin + 1
   newwave = dindgen(npix)*binsz + wavemin

   for i=0, nfiber-1 do begin
 
      scale[i] = 1.0

      splog, plugmap[i].fiberid, ' ', plugmap[i].objtype, plugmap[i].mag, $
       format = '(i4.3, a, a, f6.2, f6.2, f6.2, f6.2, f6.2)'
      fullwave = wave[*,i]
      fullspec = flux[*,i]
      fullivar = fluxivar[*,i]
      fullcombmask = bytarr(n_elements(fullspec))
      fullpixelmask = pixelmask[*,i]
      fullfibermask = fibermask[*,i]

      print, 'FULL fibermask ', fullfibermask

      outputfile = outputroot+'-' $
       +string(format='(i3.3,a)',i+1+(spectrographid-1)*320)+'.fits'

      bad = 0
      bestguess = fltarr(npix)
      bestivar = bestguess*0.0
      besterr = bestivar
      nonzero = where(fullivar GT 0.0,ngood)

      if (ngood EQ 0 OR strtrim(plugmap[i].objtype,2) EQ 'NA') then begin

         splog, 'No good points OR no plugmap entry'
         outputpixelmask = lonarr(npix) - 1L
         bad = 1

      endif else begin

;	 Now let's break sorted wavelengths into groups where
;	 pixel separations are larger than maxsep
;

         goodsort = nonzero[sort(fullwave[nonzero])]
	 wavesort = fullwave[goodsort]

	 padwave = [min(wavesort) - 2.0*maxsep, wavesort, max(wavesort) + 2.0*maxsep]

         startgroup = where(padwave[1:ngood] - padwave[0:ngood-1] GT maxsep,nstart)
         endgroup = where(padwave[2:ngood+1] - padwave[1:ngood] GT maxsep,nend)

         if (nstart NE nend) then $
           message, 'ABORT: grouping tricks did not work!'


         for igrp=0,nstart - 1 do begin
         
           ss = goodsort[startgroup[igrp]: endgroup[igrp]]    
;           bkptmin = wavesort[startgroup[igrp]]
;           bkptmax = wavesort[endgroup[igrp]]
;           nbkpt = fix((bkptmax - bkptmin)/bkptbin) + 1
           bkpt = 0

           fullbkpt = slatec_splinefit(fullwave[ss], fullspec[ss], coeff, $
             nord=nord, eachgroup=1, maxiter=20, upper=3.0, lower=3.0, $
             bkspace=bkptbin, bkpt=bkpt, invvar=fullivar[ss], mask=mask, /silent)

           inside = where(newwave GE min(bkpt) $
                        AND newwave LE max(bkpt), numinside)

           if (total(abs(coeff)) EQ 0.0 OR numinside EQ 0) then begin
             if (numinside EQ 0) then $
               splog,'WARNING: No wavelengths inside breakpoints'
             if (total(abs(coeff)) EQ 0.0) then $
               splog,'WARNING: All B-spline coefficients have been set to zero!'
           endif else begin         

             fwave = float(newwave[inside])
             bestguess[inside] = slatec_bvalu(fwave,fullbkpt,coeff)

             splog, 'Masked ', fix(total(1-mask)), ' of', $
               n_elements(mask), ' pixels'

             replace = where(mask EQ 0)

           ;-----------------------------------------------------------------
           ;  Here replace original flux values of masked pixels with b-spline
           ;  evaluations

             if (replace[0] NE -1) then begin 
               fullspec[ss[replace]] = $
                 slatec_bvalu(fullwave[ss[replace]],fullbkpt,coeff)

               fullpixelmask[ss[replace]] = fullpixelmask[ss[replace]] OR $
                                      pixelmask_bits('COMBINEREJ')
             endif

           endelse
           fullcombmask[ss] = mask

         endfor


         ;------------------------------------------------------------------------
         ; Also keep 2 spots for each file in case 2 pixels are touching newwave

         andmask = lonarr(npix) - 1L
         ormask = lonarr(npix)

         medflux = fltarr(nfiles)
         for j=0, nfiles-1 do begin
            these = where(specnum EQ j)
            if (these[0] NE -1) then begin
	       medflux[j] = djs_median(fullspec[these])

               inbetween = where(newwave GE min(fullwave[these]) AND $
                                 newwave LE max(fullwave[these]))
               if (inbetween[0] NE -1) then begin

                  ; Let's conserve inverse variance

                  totalbefore = total(fullivar[these] * fullcombmask[these])
                  result = interpol(fullivar[these] * fullcombmask[these], $
                   fullwave[these], newwave[inbetween])

                  if (total(result) GT 0.0) then begin
                    conservevariance = totalbefore / total(result)
                    bestivar[inbetween] = bestivar[inbetween] + $
                      result * conservevariance
                  endif

                  lowside = fix((fullwave[these]-wavemin)/binsz)
                  highside = lowside + 1
                  andmask[lowside]  = andmask[lowside] AND fullpixelmask[these]
                  andmask[highside] = andmask[highside] AND fullpixelmask[these]
                  ormask[lowside]   = ormask[lowside] OR fullpixelmask[these]
                  ormask[highside]  = ormask[highside] OR fullpixelmask[these]

                  ; We can also try some bitmask tricks in here
               endif

            endif
         endfor
	 splog, 'Medians:', string(format='(f7.2)',medflux)

      ;-------------------------------------------------------------------------
      ;  Replace -1's in andmask
      ;
      andmask = andmask * (andmask NE -1)

      outputpixelmask = ormask OR ishft(andmask,16)
      
         nonzero = where(bestivar GT 0.0)
         if (nonzero[0] NE -1) then $
          besterr[nonzero] = 1.0/sqrt(bestivar[nonzero])

         if (keyword_set(display)) then begin
            djs_plot, 10^fullwave, fullspec, ps=3, xr=[3700,4800], yr=[-2,10], $
             title=string(format='(i4,x,a,5(f7.3))', $
             plugmap[i].fiberid, plugmap[i].objtype, plugmap[i].mag), $
             xtitle='\lambda [A]', ytitle='Flux (10^-17 cgs)'
  
            djs_oplot, 10^newwave, bestguess,color='red',ps=10
         endif

      endelse



      newhdr = hdr

      sxaddpar, newhdr, 'OBJID', string(format='(5(i))', plugmap[i].objid)
      sxaddpar, newhdr, 'MAG', string(format='(5(f8.3))', plugmap[i].mag)
      sxaddpar, newhdr, 'RAOBJ', plugmap[i].ra, 'RA (deg) of object'
      sxaddpar, newhdr, 'DECOBJ', plugmap[i].dec, 'DEC (deg) of object'
      sxaddpar, newhdr, 'OBJTYPE', plugmap[i].objtype

      ; Need to be more clever when multiple plugmap's are used ???
      ; This is true throughout this entire routine!

      sxaddpar, newhdr, 'XFOCAL', plugmap[i].xfocal
      sxaddpar, newhdr, 'YFOCAL', plugmap[i].yfocal
      sxaddpar, newhdr, 'SPECID', plugmap[i].spectrographId

      sxaddpar, newhdr, 'PRIMTARG', plugmap[i].primtarget
      sxaddpar, newhdr, 'SECTARGE', plugmap[i].sectarget
      sxaddpar, newhdr, 'FIBERID', plugmap[i].fiberId

      sxaddpar, newhdr, 'NWORDER', 2, 'Linear-log10 coefficients'
      sxaddpar, newhdr, 'WFITTYPE', 'LOG-LINEAR', 'Linear-log10 dispersion'
      sxaddpar, newhdr, 'COEFF0', wavemin, $
       'Center wavelength (log10) of first pixel'
      ;  !!!  make sure binsz  is double for precision  !!!
      sxaddpar, newhdr, 'COEFF1', binsz, 'Log10 dispersion per pixel'
      sxaddpar, newhdr, 'REDSCAL', scale[i], $
       'Red scaling to match blue overlap', AFTER='EXPTIME'

      sxaddpar, newhdr, 'NAXIS1', n_elements(bestguess)
      sxaddpar, newhdr, 'NAXIS2', 2
      sxaddpar, newhdr, 'WAT0_001', 'system=linear'
      sxaddpar, newhdr, 'WAT1_001', $
       'wtype=linear label=Wavelength units=Angstroms'
      sxaddpar, newhdr, 'CRVAL1', wavemin, 'Iraf zero point'
      ;  !!!  make sure binsz  is double for precision  !!!
      sxaddpar, newhdr, 'CD1_1', binsz, 'Iraf dispersion'
      sxaddpar, newhdr, 'CRPIX1', 1, 'Iraf starting pixel'
      sxaddpar, newhdr, 'CTYPE1', 'LINEAR'
      sxaddpar, newhdr, 'WCSDIM', 2
      sxaddpar, newhdr, 'DC-FLAG', 1, 'Log-linear flag'

      ;  Let's also add in parameters for fibermask

;      for jj=0,n_elements(fullfibermask) - 1 do $
;        sxaddpar, newhdr, string(format='(a,i3.3)','FMASK',jj), long(fullfibermask[jj]), $
;                       'FiberMask bits,  see idlspec2d/pro/fibermask_bits.pro'


      ; This is place holding for a fiber with no counts at all!

      if (bad EQ 1) then begin
         bestguess = fltarr(npix)
         besterr = fltarr(npix)
      endif

      ; Replace NaN's in combined spectra; this should really never happen

      inff = where(finite(bestguess) EQ 0 OR finite(besterr) EQ 0)
      if (inff[0] NE -1) then begin
         splog, 'WARNING: NaNs in combined spectra ', N_elements(inff)
         bestguess[inff] = 0.0
         besterr[inff] = 0.0
      endif

      ; Interpolate over masked pixels, just for aesthetic purposes

      bestguess = djs_maskinterp(bestguess, besterr EQ 0, /const)

      ; 1st HDU is flux and error
      mwrfits, [[bestguess],[besterr]], outputfile, newhdr, /create
      ; 2nd HDU is pixelmask
      mwrfits, outputpixelmask, outputfile

      ;----------------------------------------------------------------------
      ;  For individual files,
      ;      let's append wave, flux and err for each spectrum
      ;

      if (keyword_set(individual)) then begin
        individualfile = individual+'-' $
           +string(format='(i3.3,a)',i+1+(spectrographid-1)*320)+'.fits'

        mwrfits, reform(fullwave,npixhere,nfiles), individualfile, $
                  newhdr, /create
        mwrfits, reform(fullspec,npixhere,nfiles), individualfile
        mwrfits, reform(fullivar,npixhere,nfiles), individualfile
        mwrfits, reform(fullpixelmask,npixhere,nfiles), individualfile
      endif            

      if (keyword_set(display)) then begin
         plot, 10^newwave, bestguess, /xstyle, yr=[-3,10]
         djs_oplot, 10^newwave, besterr, color='red'
      endif

   endfor

   return
end
;------------------------------------------------------------------------------
