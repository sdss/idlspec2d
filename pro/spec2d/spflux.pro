; INPUTS:
;   sciname        - Name of science and smear reduced images, for only
;                    one of the spectrographs.
;   fcalibprefix   - Prefix for flux-calibration files.
;   outdir         - optional output directory, default is current directory
;   adderr         - Additional error to add to the formal errors, as a
;                    fraction of the flux.
;
; Might be more rational to pass this routine spec-1 then spec-2 files.
;------------------------------------------------------------------------------
pro spflux, sciname, fcalibprefix, outdir=outdir, adderr=adderr

; Do we want to set this more rationally, or different for science vs. smear???
   sn2cut = 0.1

   if NOT keyword_set(outdir) then outdir=''

   ;----------
   ; Select the exposure number to use as the fiducial spectro-photo exposure.
   ; If there is a FLAVOR='smear' image, use that.

   nscience = n_elements(sciname)
   exposure = lonarr(nscience)
   flavor = strarr(nscience)
   cameras = strarr(nscience)
   camcolor = strarr(nscience)
   camspecid = strarr(nscience)
   framesn2 = fltarr(nscience)

   for i=0, nscience-1 do begin
      hdr = headfits(sciname[i])
      if (size(hdr,/tname) NE 'STRING') then begin
         splog, 'Skipping invalid or missing file '+sciname[i]
      endif else begin
         exposure[i] = sxpar(hdr, 'EXPOSURE')
         flavor[i] = strtrim(sxpar(hdr, 'FLAVOR'),2)
         cameras[i] = strtrim(sxpar(hdr, 'CAMERAS'),2)
         camcolor[i] = strmid(cameras[i], 0, 1)
         camspecid[i] = strmid(cameras[i], 1, 1)
         framesn2[i] = sxpar(hdr, 'FRAMESN2')
         if (NOT keyword_set(platestr)) then begin
            platestr = string(sxpar(hdr, 'PLATEID'), format='(i4.4)')
            mjdstr = string(sxpar(hdr, 'MJD'), format='(i5.5)')
         endif
      endelse
   endfor

   ;---------------------------------------------------------------------------
   ; LOOP THROUGH EACH SPECTROGRAPH (1 AND 2)
   ;---------------------------------------------------------------------------

   for ispec=1, 2 do begin

      qabort = 0

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

         ; Require that we have a total of 2 good frames (1 blue + 1 red)
         ; that meet the minimum S/N requirement.
         if (iblue NE -1 AND ired NE -1) then $
          qgoodexp = total([ framesn2[iblue] GT sn2cut, $
           framesn2[ired] GT sn2cut ]) EQ 2 $
         else $
          qgoodexp = 0

         if (qgoodexp) then begin
            smearscore[iexp] = $
             100 * (flavor[iblue] EQ 'smear') * framesn2[iblue] $
             + 1 * (flavor[iblue] EQ 'science') * framesn2[iblue]
            calibscore[iexp] = $
             1 * (flavor[iblue] EQ 'smear') * framesn2[iblue] $
             + 100 * (flavor[iblue] EQ 'science') * framesn2[iblue]
         endif

         splog, 'Exp #', allexpnum[iexp], ' score(smear)=', smearscore[iexp], $
          ' score(calib)=', calibscore[iexp]
      endfor

      ;----------
      ; Select best smear...

      junk = max(smearscore, ismear)
      if (smearscore[ismear] EQ 0) then begin 
         splog, 'ABORT: No valid blue+red exposure for spectro-photo smear'
         qabort = 1
      endif
 
      smearexpnum = allexpnum[ismear]
      splog, 'Select exposure ', smearexpnum, ' for spectro-photo smear'

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
         qabort = 1
      endif

if (NOT keyword_set(qabort)) then begin ; Do we abort on this spectrograph ???

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

         camid = -1
         if (iblue NE -1) then begin
           bluefile = sciname[iblue] 
           camid = camspecid[iblue]
         endif else bluefile = ''

         if (ired NE -1) then begin
           redfile = sciname[ired] 
           camid = camspecid[ired]
         endif else redfile = ''

         if (camid NE -1) then begin
           expstr = string(allexpnum[iexp], format='(i8.8)')
           corrfile = djs_filepath('spFluxcorr-'+expstr+'-'+camid+'.fits', $
              root_dir=outdir)

           if size(bluelist,/tname) EQ 'STRING' then $
              bluelist = [bluelist, bluefile] $
           else bluelist = bluefile 

           if size(redlist,/tname) EQ 'STRING'  then $
              redlist = [redlist, redfile] $
           else redlist = redfile 

           if size(corrlist,/tname) EQ 'STRING' then $
              corrlist = [corrlist, corrfile] $
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
        djs_filepath(fcalibfiles, root_dir=outdir), $
        colors=['b','r'], adderr=adderr

endif ; End-if for whether we abort on this spectrograph ???

   endfor

   return
end
;------------------------------------------------------------------------------
