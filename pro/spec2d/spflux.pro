; INPUTS:
;   sciname        - Name of science and smear reduced images, for only
;                    one of the spectrographs.
;   fcalibprefix   - Prefix for flux-calibration files.

;------------------------------------------------------------------------------
pro spflux, sciname, fcalibprefix, adderr=adderr

   if (NOT keyword_set(adderr)) then adderr = 0.03

   ;----------
   ; Select the exposure number to use as the fiducial spectro-photo exposure.
   ; If there is a FLAVOR='smear' image, use that.

   nscience = n_elements(sciname)
   exposure = lonarr(nscience)
   flavor = strarr(nscience)
   cameras = strarr(nscience)
   camcolor = strarr(nscience)
   camspecid = strarr(nscience)
   exptime = fltarr(nscience)

   for i=0, nscience-1 do begin
      hdr = headfits(sciname[i])

      exposure[i] = sxpar(hdr, 'EXPOSURE')
      flavor[i] = strtrim(sxpar(hdr, 'FLAVOR'),2)
      cameras[i] = strtrim(sxpar(hdr, 'CAMERAS'),2)
      camcolor[i] = strmid(cameras[i],0,1)
      camspecid[i] = strmid(cameras[i],1,1)
      exptime[i] = sxpar(hdr, 'EXPTIME')

      if (i EQ 0) then begin
         platestr = string(sxpar(hdr, 'PLATEID'), format='(i4.4)')
         mjdstr = string(sxpar(hdr, 'MJD'), format='(i5.5)')
      endif
   endfor

   ;---------------------------------------------------------------------------
   ; LOOP THROUGH EACH SPECTROGRAPH (1 AND 2)
   ;---------------------------------------------------------------------------

   for ispec=1, 2 do begin

      ;----------
      ; Rank each exposure number by how good it is for spectro-photometry,
      ; and select the best one.

      allexpnum = exposure[ uniq(exposure, sort(exposure)) ]
      nall = n_elements(allexpnum)
      score = fltarr(nall)
      for iexp=0, nall-1 do begin
         iblue = (where(camcolor EQ 'b' AND long(camspecid) EQ ispec $
          AND exposure EQ allexpnum[iexp]))[0]
         ired = (where(camcolor EQ 'r' AND long(camspecid) EQ ispec $
          AND exposure EQ allexpnum[iexp]))[0]
         if (iblue NE -1 AND ired NE -1) then begin
            score[iexp] = 100 * (flavor[iblue] EQ 'smear') * exptime[iblue] $
                     + 1 * (flavor[iblue] EQ 'science') * exptime[iblue]
         endif
      endfor
      junk = max(score, ibest)
      if (score[ibest] EQ 0) then $
       message, 'No valid blue+red exposure for spectro-photometry'
      bestexpnum = allexpnum[ibest]
      splog, 'Select smear image as exposure #', bestexpnum, $
       ' for spectro-', ispec

      ibcalib = (where(camcolor EQ 'b' AND long(camspecid) EQ ispec $
       AND exposure EQ bestexpnum))[0]
      ircalib = (where(camcolor EQ 'r' AND long(camspecid) EQ ispec $
       AND exposure EQ bestexpnum))[0]

      ;----------
      ; Compute the flux-calibration function from the spectro-photometric
      ; stars on this best exposure.  This is in 10^(-17) erg/ADU.

      if (ispec EQ 1) then $
       fcalibfiles = fcalibprefix + ['-b1.fits','-r1.fits'] $
      else $
       fcalibfiles = fcalibprefix + ['-b2.fits','-r2.fits']

      myfluxcalib, [sciname[ibcalib], sciname[ircalib]], $
       fcalibfiles, colors=['b','r'], adderr=adderr

      ;----------
      ; Loop through each exposure and construct the flux-correction function
      ; for each fiber to make it agree with the specro-photo (smear) image.
      ; Make a correction file for each spectrograph, but that includes both
      ; blue and red.

      for iexp=0, nall-1 do begin
         iblue = (where(camcolor EQ 'b' AND long(camspecid) EQ ispec $
          AND exposure EQ allexpnum[iexp]))[0]
         ired = (where(camcolor EQ 'r' AND long(camspecid) EQ ispec $
          AND exposure EQ allexpnum[iexp]))[0]
         if (iblue NE -1) then bluefile = sciname[iblue] $
          else bluefile = ''
         if (ired NE -1) then redfile = sciname[ired] $
          else redfile = ''

         expstr = string(allexpnum[iexp], format='(i8.8)')
         corrfile = 'spFluxcorr-' + expstr + '-' + camspecid[iblue] + '.fits'

         myfluxcorr, sciname[ibcalib], sciname[ircalib], $
          bluefile, redfile, corrfile, adderr=adderr
      endfor

   endfor

   return
end
;------------------------------------------------------------------------------
