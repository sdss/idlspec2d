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
;   mpfit()
;   mrdfits()
;   mwrfits
;   splog
;   soplot
;   splot
;
; INTERNAL SUPPORT ROUTINES:
;   fcorr_goodvector()
;   spfluxcorr_fn()
;   fcorr_chi_fn()
;   spfluxcorr_vectors()
;   spfluxcorr_solve2()
;   spfluxcorr_solve()
;
; REVISION HISTORY:
;   05-Feb-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
forward_function mpfit, fcorr_chi_fn

;------------------------------------------------------------------------------
; Return 1 if the flux-correction vector appears to be in bounds.
function fcorr_goodvector, ymult1

   ymin = min(ymult1, max=ymax)

   return, ymin GT 0.2 AND ymax LT 5.
end
;------------------------------------------------------------------------------
function spfluxcorr_fn, acoeff, ymult=ymult, yadd=yadd

   common com_fcorr_chi, npoly, nback, loglam, aflux, bflux, $
    aivar, bivar, aarr, barr

   ymult = acoeff[0] * aarr[*,0]
   for i=1, npoly-1 do ymult = ymult + acoeff[i] * aarr[*,i]
   if (nback EQ 0) then yadd = 0 $
    else yadd = acoeff[npoly:npoly+nback-1] ## barr
   yfit = ymult * aflux + yadd

   return, yfit
end
;------------------------------------------------------------------------------
; Return a vector of chi values
function fcorr_chi_fn, acoeff

   common com_fcorr_chi, npoly, nback, loglam, aflux, bflux, $
    aivar, bivar, aarr, barr

   yfit = spfluxcorr_fn(acoeff, ymult=ymult)

   ; The errors are rescaled at every function evaluation, but we
   ; only allow the errors to get smaller by up to a factor of 1e4,
   ; and we only allow them to get larger slowly (as the square root).
   ;  This should very strongly constrain the flux-corrrection vectors
   ; from going too small (or negative), or too large.
   qgood = aivar GT 0 AND bivar GT 0
   vmult = (ymult > 1e-4) * (ymult LE 1) + sqrt(ymult) * (ymult GT 1)
   totivar = qgood / ( 1./(aivar + (1-qgood)) + vmult^2/(bivar + (1-qgood)) )

   chivec2 = (bflux - yfit)^2 * totivar

   return, sqrt(chivec2)
end

;------------------------------------------------------------------------------
; Construct the polynomial or additive vectors
function spfluxcorr_vectors, loglam, npoly

   common com_fcorr_lam, minlog, maxlog

   if (npoly EQ 0) then return, 0

   npix = n_elements(loglam)
   xvector = (loglam[*] - minlog) / (maxlog - minlog)
   aarr = dblarr(npix, npoly)
   for ipoly=0, npoly-1 do aarr[*,ipoly] = xvector^ipoly

   return, aarr
end
;------------------------------------------------------------------------------
function spfluxcorr_solve2, loglam1, allflux1, allflux2, allivar1, allivar2, $
 npoly=npoly1, nback=nback1, ymult=ymult, yadd=yadd, $
 totchi2=totchi2, debug=debug

   common com_fcorr_chi, npoly, nback, loglam, aflux, bflux, $
    aivar, bivar, aarr, barr

   ; Pass arrays in the common block
   npoly = npoly1
   nback = nback1
   loglam = loglam1
   aflux = allflux1
   bflux = allflux2
   aivar = allivar1
   bivar = allivar2

   aarr = spfluxcorr_vectors(loglam, npoly)
   barr = spfluxcorr_vectors(loglam, nback)

   ; Call MPFIT to iterate on the solution for the template
   parinfo1 = {value: 0.D, fixed: 0, limited: [0b,0b], limits: [0.d0,0.d0], $
    tied: ''}
   parinfo = replicate(parinfo1, npoly+nback)
   parinfo[0].value = 1.D
   if (npoly GT 1) then parinfo[1:npoly-1].value = 1.D-6
   if (nback GT 0) then parinfo[npoly:npoly+nback-1].value = 1.D-6
   ftol = 1d-20
   gtol = 1d-20
   xtol = 1d-20

;   ; Add constraints that the multiplicative term must be in the
;   ; bounds set by LIMITS, at least at the end points.
;   limits = [0.1, 10]
;   xmin = min(xvector)
;   xmax = max(xvector)
;   if (npoly EQ 1) then begin
;      parinfo[0].limits = limits
;   endif else if (npoly EQ 2) then begin
;      tiepar = replicate(parinfo1, 2)
;      tiepar[0].tied = 'P(0) + ' + string(xmin) + ' * P(1)'
;      tiepar[1].tied = 'P(0) + ' + string(xmax) + ' * P(1)'
;      tiepar[0].limits = limits
;      tiepar[1].limits = limits
;      parinfo = [parinfo, tiepar]
;   endif else if (npoly EQ 3) then begin
;      tiepar[0].tied = 'P(0) + ' + string(xmin) + ' * P(1) + ' $
;       + string(xmin^2) + ' * P(2) * P(2)'
;      tiepar[1].tied = 'P(0) + ' + string(xmax) + ' * P(1) + ' $
;       + string(xmax^2) + ' * P(2) * P(2)'
;      tiepar[0].limits = limits
;      tiepar[1].limits = limits
;      parinfo = [parinfo, tiepar]
;   endif

   acoeff = mpfit('fcorr_chi_fn', parinfo=parinfo, perror=perror, $
    maxiter=maxiter, ftol=ftol, gtol=gtol, xtol=xtol, $
    niter=niter, status=status, /quiet)

   totchi2 = total( (fcorr_chi_fn(acoeff))^2 )

   yfit = spfluxcorr_fn(acoeff, ymult=ymult, yadd=yadd)

   if (keyword_set(debug)) then begin
      print, 'STATUS = ', status
      print, 'Best-fit coeffs = ', acoeff
      print, 'Errors = = ', perror

      set_plot,'x'
      splot,10^loglam[0:2047],smooth(bflux[0:2047], 9), xrange=[3800,9200]
      soplot,10^loglam[2048:4095],smooth(bflux[2048:4095], 9)
      soplot,10^loglam[0:2048],smooth(yfit[0:2047], 9), color='red'
      soplot,10^loglam[2048:4095],smooth(yfit[2048:4095], 9), color='red'
      cc = strupcase(get_kbrd(1))
   endif

   return, acoeff
end

;------------------------------------------------------------------------------
; Fit for BFLUX = polynomial * AFLUX + background.
; The errors are the simple quadrature sum from AIVAR and BIAVAR.
; SIGVEC = the returned sigma per pixel, even for masked pixels
;          which were excluded from the fit
; The additive term should be **per camera** ?
function spfluxcorr_solve, loglam, aflux, bflux, sqivar, mask=mask1, $
 sigvec=sigvec, npoly=npoly, nback=nback, ymult=ymult, yadd=yadd

   if (NOT keyword_set(npoly)) then npoly = 3
   if (n_elements(nback) EQ 0) then nback = 1
   if (keyword_set(mask1)) then mask = mask1 $
    else mask = 1
   wconstrain = 10 ; The weight of the constraints

   ; Set default return values
   ymult = 0
   madd = 0
   sigvec = 0 * aflux

   ; Construct the polynomial and additive vectors
   aarr = spfluxcorr_vectors(loglam, npoly)
   barr = spfluxcorr_vectors(loglam, nback)

   npix = n_elements(loglam)
   nterm = npoly + nback

   ; Compute the errors to use as the quadrature sum
   if (total(sqivar*mask GT 0) LT nterm) then return, 0

   ; Construct the empty matrices
   mmatrix = dblarr(npix+nterm, nterm)
   bvec = dblarr(npix+nterm)

   ; Put the measured values into the matrix
   for i=0, npoly-1 do $
    mmatrix[0:npix-1,i] = aflux[*] * sqivar * mask * aarr[*,i]
   for j=0, nback-1 do $
    mmatrix[0:npix-1,npoly+j] = barr[*,j]
   bvec[0:npix-1] = bflux[*] * sqivar * mask

   ; Put the constraints into the matrix: a[0] = 1
   mmatrix[npix,0] = 1 * wconstrain
   bvec[npix] = 1 * wconstrain

   ; Put the constraints into the matrix: a[1...] = 0
   if (npoly GT 1) then begin
      for i=1, npoly-1 do mmatrix[npix+i,i] = 1 * wconstrain
      bvec[npix+1:npix+npoly-1] = 0 * wconstrain
   endif

   ; Put the constraints into the matrix: b[0...] = 0
   if (nback GT 0) then begin
      for j=0, nback-1 do mmatrix[npix+npoly+j,npoly+j] = 1 * wconstrain
   endif

   ; Now invert the matrix
   mmatrixt = transpose(mmatrix)
   mm = mmatrixt # mmatrix
   if (nterm EQ 1) then begin
      mmi = 1.0 / (mm + (mm EQ 0))
   endif else begin
      svdc, mm, ww, uu, vv, /double
      mmi = 0 * vv
      ; The last term below is to protect against divide-by-zero
      ; in the degenerate case.
      for i=0L, nterm-1 do mmi[i,*] = vv[i,*] / (ww[i] + (ww[i] EQ 0))
      mmi = mmi ## transpose(uu)
   endelse

   acoeff = mmi # (mmatrixt # bvec)
   chi2 = total( (mmatrix # acoeff - bvec)^2, /double )

   ymult = acoeff[0] * aarr[*,0]
   for i=1, npoly-1 do ymult = ymult + acoeff[i] * aarr[*,i]

   yadd = 0 * aflux[*]
   if (nback GT 0) then yadd = yadd + acoeff[npoly:npoly+nback-1] ## barr

   yfit = ymult * aflux[*] + yadd
   sigvec = (yfit - bflux[*]) * sqivar

   return, yfit
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

   ; Set variables in common block
   minlog = min(loglam, max=maxlog)

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
               invsig = 0 * allivar[*,iobj,i1]
               indx = where(allivar[*,iobj,i2] GT 0 $
                AND allivar[*,iobj,i1] GT 0, ct)
               if (ct GT 0) then $
                invsig[indx] = 1. / sqrt(1./(allivar[*,iobj,i1])[indx] $
                 + 1./(allivar[*,iobj,i2])[indx])
               yfit1 = spfluxcorr_solve( loglam[*,iobj,*], $
                allflux[*,iobj,i2], allflux[*,iobj,i1], invsig, $
                mask=qgood, npoly=3, nback=2, $
                sigvec=sigvec1, ymult=ymult1, yadd=yadd1 )
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
               qcont = 1B
               lastchi2 = 0
               ; Add more polynomial terms to the fit as long as
               ; the vectors are still good (no crazy values), and
               ; the chi^2 is significantly improved (by at least 5).
               while (qcont AND npoly1 LT maxpoly) do begin
                  npoly1 = npoly1 + 1
                  thiscoeff = spfluxcorr_solve2(loglam[*,iobj,*], $
                   allflux[*,iobj,i2], allflux[*,iobj,i1], $
                   allivar[*,iobj,i2], allivar[*,iobj,i1] * qgood, $
                   npoly=npoly1, nback=nback, ymult=ymult1, yadd=yadd1, $
                   totchi2=thischi2)
                  if (fcorr_goodvector(ymult1) $
                   AND (npoly1 EQ 1 OR lastchi2-thischi2 GT 5.)) then begin
                     ymult[*,iobj,i2] = ymult1
                     yadd[*,iobj,i2] = yadd1
                     lastchi2 = thischi2
                  endif else begin
                     qcont = 0B
                  endelse
               endwhile
            endif
            splog, 'Fiber #', 320*(spectroid[0]-1)+iobj+1, $
             ' exposure #', explist[iexp], ' npoly=', npoly1
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

   return
end
;------------------------------------------------------------------------------
