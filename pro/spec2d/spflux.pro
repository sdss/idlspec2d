; INPUTS:
;   sciname        - Name of science and smear reduced images, for only
;                    one of the spectrographs.
;   fcalibprefix   - Prefix for flux-calibration files.
;   adderr         - Additional error to add to the formal errors, as a
;                    fraction of the flux.

;------------------------------------------------------------------------------
pro spflux, sciname, fcalibprefix, adderr=adderr

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
      if (size(hdr,/tname) NE 'STRING') then $
        message, 'Header is not valid for '+sciname[i]

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
      smearscore = fltarr(nall)
      calibscore = fltarr(nall)
      for iexp=0, nall-1 do begin
         iblue = (where(camcolor EQ 'b' AND long(camspecid) EQ ispec $
          AND exposure EQ allexpnum[iexp]))[0]
         ired = (where(camcolor EQ 'r' AND long(camspecid) EQ ispec $
          AND exposure EQ allexpnum[iexp]))[0]

; WE REALLY WANT THESE SCORES TO REFLECT THE BEST S/N FOR SCIENCE OR SMEAR ???
         if (iblue NE -1 AND ired NE -1) then begin
            smearscore[iexp] = $
             100 * (flavor[iblue] EQ 'smear') * exptime[iblue] $
             + 1 * (flavor[iblue] EQ 'science') * exptime[iblue]
            calibscore[iexp] = $
             1 * (flavor[iblue] EQ 'smear') * exptime[iblue] $
             + 100 * (flavor[iblue] EQ 'science') * exptime[iblue]
         endif
      endfor

      ;----------
      ; Select best smear...

      junk = max(smearscore, ismear)
      if (smearscore[ismear] EQ 0) then begin 
         splog, 'ABORT: No valid blue+red exposure for spectro-photo smear'
         return
      endif
      splog, 'Select exposure ', exposure[ismear], ' for spectro-photo smear'

      smearexpnum = allexpnum[ismear]
      splog, 'Select smear image as exposure #', smearexpnum, $
       ' for spectro-', ispec
      ibsmear = (where(camcolor EQ 'b' AND long(camspecid) EQ ispec $
       AND exposure EQ smearexpnum))[0]
      irsmear = (where(camcolor EQ 'r' AND long(camspecid) EQ ispec $
       AND exposure EQ smearexpnum))[0]

      ;----------
      ; Select best calib...

      junk = max(calibscore, icalib)
      if (calibscore[icalib] EQ 0) then begin 
         splog, 'ABORT: No valid blue+red exposure for spectro-photo calib'
         return
      endif

      calibexpnum = allexpnum[icalib]
      splog, 'Select calib image as exposure #', calibexpnum, $
       ' for spectro-', ispec
      ibcalib = (where(camcolor EQ 'b' AND long(camspecid) EQ ispec $
       AND exposure EQ calibexpnum))[0]
      ircalib = (where(camcolor EQ 'r' AND long(camspecid) EQ ispec $
       AND exposure EQ calibexpnum))[0]

      ;------------------------------------------------------------------------
      ; Loop through each exposure and construct the flux-correction function
      ; for each fiber to make it agree with the specro-photo (smear) image.
      ; Make a correction file for each spectrograph, but that includes both
      ; blue and red.
      ;------------------------------------------------------------------------

      bluelist = 0
      redlist = 0
      corrlist = 0

      for iexp=0, nall-1 do begin
         iblue = (where(camcolor EQ 'b' AND long(camspecid) EQ ispec $
          AND exposure EQ allexpnum[iexp]))[0]
         ired = (where(camcolor EQ 'r' AND long(camspecid) EQ ispec $
          AND exposure EQ allexpnum[iexp]))[0]
         if (iblue NE -1) then bluefile = sciname[iblue] $
          else bluefile = ''
         if (ired NE -1) then redfile = sciname[ired] $
          else redfile = ''

         if (iblue NE -1 AND ired NE -1) then begin
           expstr = string(allexpnum[iexp], format='(i8.8)')
           corrfile = 'spFluxcorr-' + expstr + '-' + camspecid[iblue] + '.fits'

           if keyword_set(bluelist) then bluelist = [bluelist, bluefile] $
           else bluelist = bluefile 

           if keyword_set(redlist) then redlist = [redlist, redfile] $
           else redlist = redfile 

           if keyword_set(corrlist) then corrlist = [corrlist, corrfile] $
           else corrlist = corrfile 
         endif

      endfor

;       myfluxcorr, sciname[ibsmear], sciname[irsmear], $
;        bluelist, redlist, corrlist, adderr=adderr ; ???

       fluxcorr_new, sciname[ibsmear], sciname[irsmear], $
        bluelist, redlist, corrlist

      ;----------
      ; Compute the flux-calibration function from the spectro-photometric
      ; stars on this best exposure.  This is in 10^(-17) erg/ADU.

      if (ispec EQ 1) then $
       fcalibfiles = fcalibprefix + ['-b1.fits','-r1.fits'] $
      else $
       fcalibfiles = fcalibprefix + ['-b2.fits','-r2.fits']

      myfluxcalib, [sciname[ibcalib], sciname[ircalib]], $
       fcalibfiles, colors=['b','r'], adderr=adderr

   endfor

   return
end
;------------------------------------------------------------------------------
