;+
; NAME:
;   spfluxcorr_v5
;
; PURPOSE:
;   Compute flux-correction vectors for each CCD+exposure
;
; CALLING SEQUENCE:
;   spfluxcorr_v5, objname, [ adderr=, combinedir=, ] bestexpnum=
;
; INPUTS:
;   objname    - File names (including path) for spFrame files, all from
;                either spectro-1 or spectro-2, but not both!
;   bestexpnum - Exposure number for best exposure, to which all other
;                exposures are tied.
;
; OPTIONAL INPUTS:
;   adderr     - Additional error to add to the formal errors, as a
;                fraction of the flux; default to 0.03 (3 per cent).
;   combinedir - Directory for output files
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;   djs_filepath()
;   djs_reject()
;   mrdfits()
;   mwrfits
;   solve_poly_ratio
;   splog
;   soplot
;   splot
;
; INTERNAL SUPPORT ROUTINES:
;   fcorr_goodvector()
;
; REVISION HISTORY:
;   05-Feb-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
; Return 1 if the flux-correction vector appears to be in bounds.
function fcorr_goodvector, ymult1

   ymin = min(ymult1, max=ymax)

   return, ymin GT 0.2 AND ymax LT 5.
end
;------------------------------------------------------------------------------
pro spfluxcorr_v5, objname, adderr=adderr, combinedir=combinedir, $
 bestexpnum=bestexpnum

   common com_fcorr_lam, minlog, maxlog

   if (n_elements(adderr) EQ 0) then adderr = 0.03

   ; The following parameters are used for the first pass, which is used only
   ; in rejecting points.
   maxiter1 = 5
   sigrej = 2.5

   ; The following parameters are used for the final fits, where there is
   ; no more rejection and the errors are re-scaled in each iteration.
   maxpoly = 3
   nback = 0

   t0 = systime(1)

   ;----------
   ; Get the list of spectrograph ID and camera names

   nfile = n_elements(objname)
   camname = strarr(nfile)
   camcolor = strarr(nfile)
   expnum = lonarr(nfile)
   spectroid = lonarr(nfile)
   for ifile=0, nfile-1 do begin
      spframe_read, objname[ifile], hdr=hdr
      camname[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
      camcolor[ifile] = strmid(camname[ifile],0,1)
      spectroid[ifile] = strmid(camname[ifile],1,1)
      expnum[ifile] = sxpar(hdr, 'EXPOSURE')
   endfor

   explist = expnum[uniq(expnum, sort(expnum))]
   nexp = n_elements(explist)

   ;----------
   ; Read the plug-map from the first file.
   ; We only need this to know which are SKY fibers.

   spframe_read, objname[0], plugmap=plugmap
   isky = where(strmatch(plugmap.objtype,'SKY*'), nsky)

   ;----------
   ; Get the fiducial wavelength mapping from the "best" exposure

   ibest_b = (where(camcolor EQ 'b' AND expnum EQ bestexpnum))[0]
   ibest_r = (where(camcolor EQ 'r' AND expnum EQ bestexpnum))[0]
   spframe_read, objname[ibest_b], loglam=loglam1
   dims = size(loglam1, /dimens)
   npix = dims[0]
   nobj = dims[1]
   loglam = fltarr(npix, nobj, 2)
   loglam[*,*,0] = loglam1
   spframe_read, objname[ibest_r], loglam=loglam1
   loglam[*,*,1] = loglam1
   loglam1 = 0 ; clear memory

   ;----------
   ; Read all the spectra + errors, and re-sample to the fiducial wavelengths

   allflux = fltarr(npix,nobj,nfile)
   allivar = fltarr(npix,nobj,nfile)

   for ifile=0L, nfile-1 do begin
      splog, 'Reading + rebinning raw spectra ', objname[ifile]
      ; Should we also read in the mask and reject around bright sky, etc ?
      icolor = (camcolor[ifile] EQ 'b') ? 0 : 1

      ; Read the raw spectra for this file
      spframe_read, objname[ifile], objflux=objflux1, objivar=objivar1, $
       wset=wset1, loglam=loglam1, adderr=adderr

      ; Re-normalize the flux to ADU/(dloglam)
      binsz = 1.0d-4
      correct_dlam, objflux1, objivar1, wset1, dlam=binsz

      ; Apply the spectro-calib vector for this file
      calibfile = djs_filepath(string(camname[ifile], expnum[ifile], $
       format='("spFluxcalib-", a2, "-", i8.8, ".fits")'), $
       root_dir=combinedir)
      calibfile = (findfile(calibfile+'*'))[0]
      calibfac = mrdfits(calibfile, 0, /silent)
      minval = 0.05 * mean(calibfac)
      divideflat, objflux1, invvar=objivar1, calibfac, minval=minval

      if (expnum[ifile] EQ bestexpnum) then begin
         allflux[*,*,ifile] = objflux1
         allivar[*,*,ifile] = objivar1
      endif else begin
         for iobj=0L, nobj-1 do begin
            ; We have to call COMBINE1FIBER with ascending wavelengths...
            if (loglam1[0,iobj] GT loglam1[1,iobj]) then begin
               combine1fiber, reverse(loglam1[*,iobj]), $
                reverse(objflux1[*,iobj]), reverse(objivar1[*,iobj]), $
                newloglam=reverse(loglam[*,iobj,icolor]), $
                newflux=newflux1, newivar=newivar1
               allflux[*,iobj,ifile] = reverse(newflux1)
               allivar[*,iobj,ifile] = reverse(newivar1)
            endif else begin
               combine1fiber, loglam1[*,iobj], $
                objflux1[*,iobj], objivar1[*,iobj], $
                newloglam=loglam[*,iobj,icolor], $
                newflux=newflux1, newivar=newivar1
               allflux[*,iobj,ifile] = newflux1
               allivar[*,iobj,ifile] = newivar1
            endelse
         endfor
      endelse
   endfor

   ;----------
   ; Loop over each object, solving for the correction vectors

   ; Rescale the wavelengths to be in the range [0,1]
   minlog = min(loglam, max=maxlog)
   xarray = (loglam - minlog) / (maxlog - minlog)

   ymult = fltarr(npix,nobj,nfile) + 1.
   yadd = fltarr(npix,nobj,nfile)

   i1 = [ibest_b,ibest_r]
   for iobj=0L, nobj-1 do begin
      outmask = 0

      ; This first iteration loop allows generous fitting parameters,
      ; and is primarily to reject outlier points.
      iiter = 0L
      while (iiter LT maxiter1) do begin
         ; Loop over exposures
         sigvec = 0 * allflux[*,iobj,i1]
         for iexp=0L, nexp-1 do begin
            if (explist[iexp] EQ bestexpnum) then begin
               ymult1 = 0
               yadd1 = 0
            endif else begin
               i_b = where(camcolor EQ 'b' AND expnum EQ explist[iexp], ct1)
               i_r = where(camcolor EQ 'r' AND expnum EQ explist[iexp], ct2)
               i2 = [i_b,i_r]
               if (keyword_set(outmask)) then qgood = outmask $
                else qgood = 1
               thisivar = 0 * allivar[*,iobj,i1]
               indx = where(allivar[*,iobj,i2] GT 0 $
                AND allivar[*,iobj,i1] GT 0, ct)
               if (ct GT 0) then $
                thisivar[indx] = 1. / (1./(allivar[*,iobj,i1])[indx] $
                 + 1./(allivar[*,iobj,i2])[indx])
               solve_poly_ratio, xarray[*,iobj,*], $
                allflux[*,iobj,i2], allflux[*,iobj,i1], thisivar*qgood, $
                npoly=3, nback=2, yfit=yfit1, ymult=ymult1, yadd=yadd1

               ; SIGVEC1 = the returned sigma per pixel, even for masked pixels
               ; which were excluded from the fit
               sigvec1 = (yfit1 - (allflux[*,iobj,i1])[*]) * sqrt(thisivar)
               sigvec = sigvec > abs(sigvec1)
            endelse

            if (keyword_set(ymult1)) then begin
               ymult[*,iobj,i2] = ymult1
               yadd[*,iobj,i2] = yadd1
            endif
         endfor ; End loop over exposures

         if (djs_reject(sigvec, 0*sigvec, invvar=0*sigvec+1, outmask=outmask, $
          lower=sigrej, upper=sigrej)) $
          then iiter = maxiter1 $
         else $
          iiter = iiter + 1
      endwhile

      ; This second iteration rescales the errors.
      ; No more pixels will be rejected in this loop.
      ; Loop over exposures
      for iexp=0L, nexp-1 do begin
         if (explist[iexp] NE bestexpnum) then begin
            i_b = where(camcolor EQ 'b' AND expnum EQ explist[iexp], ct1)
            i_r = where(camcolor EQ 'r' AND expnum EQ explist[iexp], ct2)
            i2 = [i_b,i_r]

            ; Set default values in caes all fits are bad
            ymult[*,iobj,i2] = 1
            yadd[*,iobj,i2] = 0

            qgood =allivar[*,iobj,i2] GT 0 $
             AND allivar[*,iobj,i1] GT 0 $
             AND outmask GT 0
            igood = where(qgood, ct)
            npoly1 = 0
            if (ct GT 0) then begin
               lastchi2 = 0
               lastnpoly = 0
               inparams = 0
               ; Loop through adding more polynomial terms up to MAXPOLY.
               ; Replace the fluxing vectors with the new ones as long
               ; as the new ones are still good (no crazy values), and
               ; the chi^2 is significantly improved (by at least 5).
               while (npoly1 LT maxpoly) do begin
                  npoly1 = npoly1 + 1

                  ; Do a 1st pass re-fitting without scaling the errors.
                  thisivar = 0 * allivar[*,iobj,i1]
                  indx = where(allivar[*,iobj,i2] GT 0 $
                   AND allivar[*,iobj,i1] GT 0, ct)
                  if (ct GT 0) then $
                   thisivar[indx] = 1. / (1./(allivar[*,iobj,i1])[indx] $
                    + 1./(allivar[*,iobj,i2])[indx])
                  solve_poly_ratio, xarray[*,iobj,*], $
                   allflux[*,iobj,i2], allflux[*,iobj,i1], $
                   thisivar*qgood, npoly=npoly1, nback=0, $
                   ymult=ymult1, yadd=yadd1, acoeff=inparams

                  ; Now fit using the non-linear code that rescales the errors.
                  solve_poly_ratio, xarray[*,iobj,*], $
                   allflux[*,iobj,i2], allflux[*,iobj,i1], $
                   allivar[*,iobj,i2], allivar[*,iobj,i1] * qgood, $
                   npoly=npoly1, nback=nback, inparams=inparams, $
                   yfit=yfit1, ymult=ymult1, yadd=yadd1, acoeff=thiscoeff, $
                   totchi2=thischi2, status=status, perror=perror

                  if (fcorr_goodvector(ymult1) $
                   AND (npoly1 EQ 1 OR lastchi2-thischi2 GT 5.)) then begin
                     ymult[*,iobj,i2] = ymult1
                     yadd[*,iobj,i2] = yadd1
                     lastchi2 = thischi2
                     lastnpoly = npoly1
                     inparams = thiscoeff
                  endif
               endwhile

               ; Optional debugging plots
;               if (keyword_set(debug)) then begin
;                  print, 'STATUS = ', status
;                  print, 'Best-fit coeffs = ', acoeff
;                  print, 'Errors = = ', perror
;
;                  set_plot,'x'
;                  splot, 10^loglam[0:2047], smooth(bflux[0:2047], 9), $
;                   xrange=[3800,9200]
;                  soplot, 10^loglam[2048:4095], smooth(bflux[2048:4095], 9)
;                  soplot, 10^loglam[0:2048], smooth(yfit1[0:2047], 9), $
;                   color='red'
;                  soplot, 10^loglam[2048:4095], smooth(yfit1[2048:4095], 9), $
;                   color='red'
;                  cc = strupcase(get_kbrd(1))
;               endif

            endif
            splog, 'Fiber #', 320*(spectroid[0]-1)+iobj+1, $
             ' exposure #', explist[iexp], ' npoly=', lastnpoly
         endif
      endfor

   endfor

   ;----------
   ; Force the sky fibers to have no correction (meaning values of unity)

   if (nsky GT 0) then begin
      ymult[*,isky,*] = 1
      yadd[*,isky,*] = 0
   endif

   ;----------
   ; Write the output files

   for ifile=0L, nfile-1 do begin
      corrfile = djs_filepath(string(camname[ifile], expnum[ifile], $
       format='("spFluxcorr-", a2, "-", i8.8, ".fits")'), $
       root_dir=combinedir)
      mwrfits, ymult[*,*,ifile], corrfile, /create
      mwrfits, yadd[*,*,ifile], corrfile
      spawn, ['gzip','-f',corrfile], /noshell
   endfor

   splog, 'Time to compute fluxcorr vectors = ', systime(1)-t0, ' sec'

   return
end
;------------------------------------------------------------------------------
