;+
; NAME:
;   combine2dout
;
; PURPOSE:
;   Combine several reduced frames of the same objects
;
; CALLING SEQUENCE:
;   combine2dout, filenames, outputname, spectrographid, $
;    binsz, zeropoint, [ nord=, wavemin=, $
;    bkptbin=, /everyn, window=, maxsep= ]
;
; INPUTS:
;   filenames      - Name(s) of files written by SPREDUCE
;   outputname     - Name for all output file
;   spectrographid - Spectrograph ID (1 or 2) for use in computing fiber
;                    number for output file name; fibers from spectro-1
;                    are numbered from 1, those from spectro-2 from 321.
;
; REQUIRED KEYWORDS:
;
; OPTIONAL KEYWORDS:
;   binsz          - ???
;   zeropoint      - ???
;   nord           - ???
;   wavemin        - ???
;   bkptbin        - ???
;   everyn         - ???
;   window         - ???
;   maxsep         - ???
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This routine also outputs original 2048 spectra with mask pixels
;   replaced with their b-spline values.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
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
pro combine1fiber, specnum, nfinalpix, $
 fullwave, fullspec, fullivar, fullpixelmask, fullfibermask, $
 newwave, bestguess, besterr, outputpixelmask, $
 nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep

   if (NOT keyword_set(nord)) then nord = 3
   if (NOT keyword_set(maxsep)) then maxsep = 2.0 * binsz
   if (NOT keyword_set(bkptbin)) then bkptbin = 1.2 * binsz

   fullcombmask = bytarr(n_elements(fullspec))

   print, 'FULL fibermask ', fullfibermask

   bestguess = fltarr(nfinalpix)
   bestivar = bestguess*0.0
   besterr = bestivar

   nonzero = where(fullivar GT 0.0, ngood)

   if (ngood EQ 0) then begin

      splog, 'No good points'
      outputpixelmask = - 1L
      return

   endif else begin

      ; Now let's break sorted wavelengths into groups where
      ; pixel separations are larger than maxsep

      goodsort = nonzero[sort(fullwave[nonzero])]
      wavesort = fullwave[goodsort]

      padwave = [min(wavesort) - 2.0*maxsep, wavesort, $
       max(wavesort) + 2.0*maxsep]

      startgroup = where(padwave[1:ngood] - padwave[0:ngood-1] GT maxsep,nstart)
      endgroup = where(padwave[2:ngood+1] - padwave[1:ngood] GT maxsep,nend)

      if (nstart NE nend) then $
       message, 'ABORT: grouping tricks did not work!'

      for igrp=0, nstart-1 do begin
         
         ss = goodsort[startgroup[igrp]: endgroup[igrp]]    
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

            bestguess[inside] = slatec_bvalu(newwave[inside],fullbkpt,coeff)
fwave = float(newwave[inside])
print,10^[min(fwave),max(fwave)]

            splog, 'Masked ', fix(total(1-mask)), ' of', $
             n_elements(mask), ' pixels'

            replace = where(mask EQ 0)

            ;-----------------------------------------------------------------
            ;  Here replace original flux values of masked pixels with b-spline
            ;  evaluations ??? Not needed if INDIVIDUAL keyword gone???

            if (replace[0] NE -1) then begin 
               fullspec[ss[replace]] = $
                slatec_bvalu(fullwave[ss[replace]],fullbkpt,coeff)

               fullpixelmask[ss[replace]] = fullpixelmask[ss[replace]] OR $
                pixelmask_bits('COMBINEREJ')
            endif

         endelse
         fullcombmask[ss] = mask

      endfor

      ;---------------------------------------------------------------------
      ; Also keep 2 spots for each file in case 2 pixels are touching newwave

      andmask = lonarr(nfinalpix) - 1L
      ormask = lonarr(nfinalpix)

      for j=0, max(specnum) do begin
         these = where(specnum EQ j)
         if (these[0] NE -1) then begin

            inbetween = where(newwave GE min(fullwave[these]) AND $
                              newwave LE max(fullwave[these]))
            if (inbetween[0] NE -1) then begin

               ; Conserve inverse variance by doing a linear interpolation
               ; on that quantity.

               result = interpol(fullivar[these] * fullcombmask[these], $
                fullwave[these], newwave[inbetween])

               ; Grow the fullcombmask below to reject any new sampling
               ; containing even a partial masked pixel.

               smask = interpol(float(fullcombmask[these]), $
                fullwave[these], newwave[inbetween])
               ibad = where(smask LT 1.0)
               if (ibad[0] NE -1) then result[ibad] = 0

               bestivar[inbetween] = bestivar[inbetween] + result

               lowside = fix((fullwave[these]-newwave[0])/binsz)
               highside = lowside + 1
               andmask[lowside]  = andmask[lowside] AND fullpixelmask[these]
               andmask[highside] = andmask[highside] AND fullpixelmask[these]
               ormask[lowside]   = ormask[lowside] OR fullpixelmask[these]
               ormask[highside]  = ormask[highside] OR fullpixelmask[these]
            endif

         endif
      endfor
      splog, 'Medians:', string(format='(f7.2)', $
       djs_median(fullspec[these]))

      ;----------
      ;  Replace -1's in andmask

      andmask = andmask * (andmask NE -1)

      outputpixelmask = ormask OR ishft(andmask,16)
      
      nonzero = where(bestivar GT 0.0)
      if (nonzero[0] NE -1) then $
       besterr[nonzero] = 1.0 / sqrt(bestivar[nonzero])

   endelse

   ; Replace NaN's in combined spectra; this should really never happen

   inff = where(finite(bestguess) EQ 0 OR finite(besterr) EQ 0)
   if (inff[0] NE -1) then begin
      splog, 'WARNING: NaNs in combined spectra ', N_elements(inff)
      bestguess[inff] = 0.0
      besterr[inff] = 0.0
   endif

   ; Interpolate over masked pixels, just for aesthetic purposes

   bestguess = djs_maskinterp(bestguess, besterr EQ 0, /const)

   return
end
;------------------------------------------------------------------------------

pro combine2dout, filenames, outputname, spectrographid, $
 binsz, zeropoint, nord=nord, wavemin=wavemin, $
 bkptbin=bkptbin, everyn=everyn, window=window, maxsep=maxsep

   ; Initial binning was 69 km/s per pixel (with a sigma of 1.0 pixel)
   ; 69.02977415 km/s is log lambda 10^-4

   if (NOT keyword_set(binsz)) then binsz = 1.0d-4 $
    else binsz = double(binsz)

   if (NOT keyword_set(zeropoint)) then zeropoint = 3.5D

   ;----------
   ; Sort filenames such that this list contains first the blue then the red

   filenames = filenames[sort(filenames)]

   nfiles = N_elements(filenames)
   if (nfiles EQ 0) then return

   redfiles = 0
   bluefiles = 0

   ;----------
   ; Somewhere up here we should add a flag combination if gives a
   ; positive intersection with fibermask... To reject single fibers. ???

   ;---------------------------------------------------------------------------

   exptime = 0

   ;---------------------------------------------------------------------------
   ; Loop through each 2D output and read in the data

   for ifile=0, nfiles-1 do begin
      tempflux = mrdfits(filenames[ifile], 0, hdr)
      tempivar = mrdfits(filenames[ifile], 1)
      tempplug = mrdfits(filenames[ifile], 2)
      tempwset = mrdfits(filenames[ifile], 3)

      traceset2xy, tempwset, pixnorm, tempwave

      dims = size(tempflux, /dimens)
      npix = dims[0]
      nfiber = dims[1]

      if (keyword_set(nfibersave)) then begin
         if (nfiber NE nfibersave) then begin
            splog, 'ABORT: Different files have different number of fibers!'
            return
         endif
      endif
      nfibersave = nfiber

      ;----------
      ; Read in pixelmask

      tt = mrdfits(filenames[ifile], 4)
      if (keyword_set(pixelmask)) then pixelmask = [pixelmask, tt] $
       else pixelmask = tt

      ;----------
      ; Read in fibermask

      tt = mrdfits(filenames[ifile], 5)
      if (keyword_set(fibermask)) then fibermask = [fibermask, transpose(tt)] $
       else fibermask = transpose(tt)

      if (keyword_set(window)) then $
       tempivar[0:window] = tempivar[0:window] * findgen(window+1) / window

      if (keyword_set(flux)) then flux = [flux, tempflux] $
       else flux = tempflux

      if (keyword_set(fluxivar)) then fluxivar = [fluxivar, tempivar] $
       else fluxivar = tempivar

      if (keyword_set(wave)) then wave = [wave, tempwave] $
       else wave = tempwave

      if (keyword_set(specnum)) then specnum = [specnum, bytarr(npix) + ifile] $
       else specnum = bytarr(npix) + ifile

      isred = (strpos(sxpar(hdr, 'CAMERAS'),'r') EQ 0)
      if (keyword_set(bluered)) then bluered = [bluered, bytarr(npix) + isred] $
       else bluered = bytarr(npix) + isred

      if (isred) then redfiles = redfiles + 1 
      if (strpos(sxpar(hdr, 'CAMERAS'),'b') EQ 0) then $
                 bluefiles = bluefiles + 1 

      if (NOT keyword_set(label)) then $
       label = makelabel(hdr) $
      else $
       label = [label, makelabel(hdr)]

      exptime = exptime + sxpar(hdr, 'EXPTIME')

      ;----------
      ; Read in the plug map (need to only do this once)

      if (ifile EQ 0) then begin
         plugmap  = mrdfits(filenames[ifile], 2)
      endif

   endfor

   splog, 'Found '+string(redfiles)+' red files'
   splog, 'Found '+string(bluefiles)+' blue files'
   if (redfiles LT 2 OR bluefiles LT 2) then begin
      splog, 'ABORT: For the time being, I expect at least 2 of each red and blue to combine'
      return
   endif

;set_plot,'x'
;ifib=193
;colorv=['default','red','green','blue','magenta','yellow','cyan','default']
;lam = 10^wave
;splot, lam[where(specnum EQ 0),ifib], $
; flux[where(specnum EQ 0),ifib], color=colorv[0]
;for jnum=1, nfiles-1 do $
; soplot, lam[where(specnum EQ jnum),ifib], $
;  flux[where(specnum EQ jnum),ifib], color=colorv[jnum]
;
;splot, lam[where(specnum EQ 0),ifib], $
; pixelmask[where(specnum EQ 0),ifib], color=colorv[0]
;for jnum=1, nfiles-1 do $
; soplot, lam[where(specnum EQ jnum),ifib], $
;  pixelmask[where(specnum EQ jnum),ifib], color=colorv[jnum]
;stop ; ???

   totalpix = (size(flux, /dimens))[0]

   redpix = where(bluered, numred)
   bluepix = where(bluered EQ 0, numblue)
   if (numblue GT 0 AND numred GT 0) then $
    exptime = exptime * 0.5

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

   nfinalpix = spotmax - spotmin + 1
   newwave = dindgen(nfinalpix) * binsz + wavemin

nfiber = 320 ; ???
   finalflux = fltarr(nfinalpix, nfiber)
   finalerr = fltarr(nfinalpix, nfiber)
   finalpixelmask = fltarr(nfinalpix, nfiber)

   ;---------------------------------------------------------------------------
   ; Combine each fiber, one at a time

   for ifiber=0, nfiber-1 do begin
;   for ifiber=0, 10 do begin ; ???
      splog, plugmap[ifiber].fiberid, ' ', plugmap[ifiber].objtype, $
       plugmap[ifiber].mag, $
       format = '(i4.3, a, a, f6.2, f6.2, f6.2, f6.2, f6.2)'

      if (strtrim(plugmap[ifiber].objtype,2) NE 'NA') then begin
         combine1fiber, specnum, nfinalpix, $
          wave[*,ifiber], flux[*,ifiber], fluxivar[*,ifiber], $
          pixelmask[*,ifiber], fibermask[*,ifiber], $
          newwave, bestguess, besterr, outputpixelmask, $
          nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep

         finalflux[*,ifiber] = bestguess
         finalerr[*,ifiber] = besterr
         finalpixelmask[*,ifiber] = outputpixelmask
      endif else begin
         splog, 'No plugmap entry'
         finalpixelmask[*,ifiber] = -1L
      endelse
   endfor

   ;---------------------------------------------------------------------------
   ; Create the output header

   sxaddpar, hdr, 'CREATORS', 'Burles & Schlegel (1999) IDLspec', after='SDSS'

   ncoeff = sxpar(hdr, 'NWORDER')
   for i=2, ncoeff-1 do sxdelpar, hdr, 'COEFF'+strtrim(string(i),2)

   ; Now get rid of exposure, and add list of exposures

   sxdelpar, hdr, 'EXPOSURE'
   sxdelpar, hdr, 'SEQID'
   sxaddpar, hdr, 'NEXP', nfiles, $
    'Number of exposures in this file', after='TELESCOP'
   for i=0,nfiles-1 do $
    sxaddpar, hdr, 'EXPID'+strtrim(string(i),2), label[i], $
     'ID string for exposure '+strtrim(string(i),2), before='EXPTIME'

   sxaddpar, hdr, 'EXPTIME', exptime, 'Total exposure time (seconds)'
   sxaddpar, hdr, 'COMBINE2', systime(), $
    'COMBINE2DOUT finished', after='EXPTIME'

   sxaddpar, hdr, 'NWORDER', 2, 'Linear-log10 coefficients'
   sxaddpar, hdr, 'WFITTYPE', 'LOG-LINEAR', 'Linear-log10 dispersion'
   sxaddpar, hdr, 'COEFF0', wavemin, $
    'Central wavelength (log10) of first pixel'
   sxaddpar, hdr, 'COEFF1', binsz, 'Log10 dispersion per pixel'

   sxaddpar, hdr, 'NAXIS1', n_elements(bestguess)
   sxaddpar, hdr, 'NAXIS2', 2
   sxaddpar, hdr, 'WAT0_001', 'system=linear'
   sxaddpar, hdr, 'WAT1_001', $
    'wtype=linear label=Wavelength units=Angstroms'
   sxaddpar, hdr, 'CRVAL1', wavemin, 'Iraf zero point'
   sxaddpar, hdr, 'CD1_1', binsz, 'Iraf dispersion'
   sxaddpar, hdr, 'CRPIX1', 1, 'Iraf starting pixel'
   sxaddpar, hdr, 'CTYPE1', 'LINEAR'
   sxaddpar, hdr, 'WCSDIM', 2
   sxaddpar, hdr, 'DC-FLAG', 1, 'Log-linear flag'

   ;---------------------------------------------------------------------------
   ; 1st HDU is flux
   mwrfits, finalflux, outputname, hdr, /create

   ; 2nd HDU is error
   mwrfits, finalerr, outputname

   ; 3rd HDU is pixelmask
   mwrfits, finalpixelmask, outputname

   return
end
;------------------------------------------------------------------------------
