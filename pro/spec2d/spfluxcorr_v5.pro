
;   objname    - These must be spFrame file names all from either spectro-1
;                or spectro-2, but not both!
;   adderr     - Additional error to add to the formal errors, as a
;                fraction of the flux; default to 0.03 (3 per cent).
;------------------------------------------------------------------------------
forward_function mpfit, fcorr_chi_fn

;------------------------------------------------------------------------------
; Return a single chi value
function fcorr_chi_fn, acoeff

   common com_fcorr_chi, npoly, nback, loglam, aflux, bflux, $
    aivar, bivar, aarr, barr

   ymult = acoeff[0] * aarr[*,0]
   for i=1, npoly-1 do ymult = ymult + acoeff[i] * aarr[*,i]

   ; The errors are rescaled at every function evaluation.
   qgood = aivar GT 0 AND bivar GT 0 AND ymult NE 0
   totivar = qgood / ( 1./(aivar + (1-qgood)) + ymult^2/(bivar + (1-qgood)) )

   chivec2 = (bflux - ymult * aflux)^2 * totivar

   return, sqrt(chivec2)
end

;------------------------------------------------------------------------------
function spfluxcorr_solve2, allloglam, allflux1, allflux2, allivar1, allivar2, $
 npoly=npoly1, nback=nback1, ymult=ymult, yadd=yadd, $
 debug=debug

   common com_fcorr_chi, npoly, nback, loglam, aflux, bflux, $
    aivar, bivar, aarr, barr

   ; Pass arrays in the common block
   npoly = npoly1
   nback = nback1
   loglam = allloglam
   aflux = allflux1
   bflux = allflux2
   aivar = allivar1
   bivar = allivar2

   ; Construct the polynomial and additive vectors
   npix = n_elements(loglam)
   nterm = npoly + nback
   minlog = min(loglam, max=maxlog)
   xvector = (loglam[*] - minlog) / (maxlog - minlog)
   aarr = dblarr(npix, npoly)
   for ipoly=0, npoly-1 do aarr[*,ipoly] = xvector^ipoly
   if (nback GT 0) then begin
      barr = dblarr(npix, nback)
      for iback=0, nback-1 do barr[*,iback] = xvector^ipoly
   endif

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

limits = [0.1, 10] ; ???
   ; Add constraints that the multiplicative term must be in the
   ; bounds set by LIMITS, at least at the end points.
   xmin = min(xvector)
   xmax = max(xvector)
   if (npoly EQ 1) then begin
      parinfo[0].limits = limits
   endif else if (npoly EQ 2) then begin
      tiepar = replicate(parinfo1, 2)
      tiepar[0].tied = 'P(0) + ' + string(xmin) + ' * P(1)'
      tiepar[1].tied = 'P(0) + ' + string(xmax) + ' * P(1)'
      tiepar[0].limits = limits
      tiepar[1].limits = limits
      parinfo = [parinfo, tiepar]
   endif else if (npoly EQ 3) then begin
      tiepar[0].tied = 'P(0) + ' + string(xmin) + ' * P(1) + ' $
       + string(xmin^2) + ' * P(2) * P(2)'
      tiepar[1].tied = 'P(0) + ' + string(xmax) + ' * P(1) + ' $
       + string(xmax^2) + ' * P(2) * P(2)'
      tiepar[0].limits = limits
      tiepar[1].limits = limits
      parinfo = [parinfo, tiepar]
   endif

   acoeff = mpfit('fcorr_chi_fn', parinfo=parinfo, perror=perror, $
    maxiter=maxiter, ftol=ftol, gtol=gtol, xtol=xtol, $
    niter=niter, status=status, /quiet)

;   print, 'STATUS = ', status
   print, 'Best-fit coeffs = ', acoeff
;   print, 'Errors = = ', perror

   ymult = acoeff[0] * aarr[*,0]
   for i=1, npoly-1 do ymult = ymult + acoeff[i] * aarr[*,i]
yadd = 0 ; ???

   yfit = ymult * aflux[*] + yadd
if (keyword_set(debug)) then begin ; ???
set_plot,'x'
splot,10^loglam[0:2047],smooth(bflux[0:2047], 9), xr=[3800,9200]
soplot,10^loglam[2048:4095],smooth(bflux[2048:4095], 9)
soplot,10^loglam[0:2048],smooth(yfit[0:2047], 9), color='red'
soplot,10^loglam[2048:4095],smooth(yfit[2048:4095], 9), color='red'
cc = strupcase(get_kbrd(1))
endif

   return, yfit
end

;------------------------------------------------------------------------------
; Fit for BFLUX = polynomial * AFLUX + background.
; The errors are the simple quadrature sum from AIVAR and BIAVAR.
; SIGVEC = the returned sigma per pixel, even for masked pixels
;          which were excluded from the fit
; The additive term should be **per camera** !!!???
function spfluxcorr_solve, loglam, aflux, bflux, sqivar, mask=mask1, $
 sigvec=sigvec, npoly=npoly, nback=nback, ymult=ymult, yadd=yadd, $
 debug=debug

   if (NOT keyword_set(npoly)) then npoly = 3
   if (n_elements(nback) EQ 0) then nback = 1
   if (keyword_set(mask1)) then mask = mask1 $
    else mask = 1
   constraint = 10 ; The weight of the constraints ???

   ; Set default return values
   ymult = 0
   madd = 0

   ; Construct the polynomial and additive vectors
   npix = n_elements(loglam)
   nterm = npoly + nback
   minlog = min(loglam, max=maxlog)
   xvector = (loglam[*] - minlog) / (maxlog - minlog)
   aarr = dblarr(npix, npoly)
   for ipoly=0, npoly-1 do aarr[*,ipoly] = xvector^ipoly
   if (nback GT 0) then begin
      barr = dblarr(npix, nback)
      for iback=0, nback-1 do barr[*,iback] = xvector^ipoly
   endif

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
   mmatrix[npix,0] = 1 * constraint
   bvec[npix] = 1 * constraint

   ; Put the constraints into the matrix: a[1...] = 0
   if (npoly GT 1) then begin
      for i=1, npoly-1 do mmatrix[npix+i,i] = 1 * constraint
      bvec[npix+1:npix+npoly-1] = 0 * constraint
   endif

   ; Put the constraints into the matrix: b[0...] = 0
   if (nback GT 0) then begin
      for j=0, nback-1 do mmatrix[npix+npoly+j,npoly+j] = 1 * constraint
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
print,acoeff

   ymult = acoeff[0] * aarr[*,0]
   for i=1, npoly-1 do ymult = ymult + acoeff[i] * aarr[*,i]

   yadd = 0 * aflux[*]
   if (nback GT 0) then yadd = yadd + acoeff[npoly:npoly+nback-1] ## barr

   yfit = ymult * aflux[*] + yadd
   sigvec = (yfit - bflux[*]) * sqivar

if (keyword_set(debug)) then begin ; ???
set_plot,'x'
splot,10^loglam[0:2047],djs_median(bflux[0:2047], $
 width=7,boundary='reflect'),xr=[3800,9200]
soplot,10^loglam[2048:4095],djs_median(bflux[2048:4095], $
 width=7,boundary='reflect')
soplot,10^loglam[0:2048],djs_median(yfit[0:2047], $
 width=7,boundary='reflect'),color='red'
soplot,10^loglam[2048:4095],djs_median(yfit[2048:4095], $
 width=7,boundary='reflect'),color='red'
cc = strupcase(get_kbrd(1))
endif

   return, yfit
end
;------------------------------------------------------------------------------
pro spfluxcorr_v5, objname, adderr=adderr, combinedir=combinedir, $
 bestexpnum=bestexpnum

   if (n_elements(adderr) EQ 0) then adderr = 0.03
   nfile = n_elements(objname)

   ;----------
   ; Get the list of spectrograph ID and camera names

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
; Should we also read in the mask and reject around bright sky, etc ???
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
print,'zzz',iobj
            ; We have to call COMBINE1FIBER with ascending wavelengths...
            if (loglam1[0,iobj] GT loglam[1,iobj]) then begin
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

   ymult = fltarr(npix,nobj,nfile) + 1.
   yadd = fltarr(npix,nobj,nfile)

maxiter1 = 5 ; ???
sigrej = 2.5 ; ???
   i1 = [ibest_b,ibest_r]
   for iobj=0L, nobj-1 do begin
      outmask = 0

      ; This first iteration loop allows generous fitting parameters,
      ; and is primarily to reject outlier points.
      iiter = 0L
      while (iiter LT maxiter1) do begin
splog,'Object', iobj, ' iter ', iiter
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
splog,'Object', iobj
      ; Loop over exposures
      for iexp=0L, nexp-1 do begin
         if (explist[iexp] EQ bestexpnum) then begin
            ymult1 = 0
            yadd1 = 0
         endif else begin
; We want to add more parameters as long as chi^2 is significantly improved???
            i_b = where(camcolor EQ 'b' AND expnum EQ explist[iexp], ct1)
            i_r = where(camcolor EQ 'r' AND expnum EQ explist[iexp], ct2)
            i2 = [i_b,i_r]
            qgood =allivar[*,iobj,i2] GT 0 $
             AND allivar[*,iobj,i1] GT 0 $
             AND outmask GT 0
            igood = where(qgood, ct)
            if (ct GT 0) then $
             yfit1 = spfluxcorr_solve2(loglam[*,iobj,*], $
              allflux[*,iobj,i2], allflux[*,iobj,i1], $
              allivar[*,iobj,i2], allivar[*,iobj,i1] * qgood, $
              npoly=2, nback=0, ymult=ymult1, yadd=yadd1)
;              npoly=2, nback=0, ymult=ymult1, yadd=yadd1, /debug) ; ???
               ; The errors are bounded not to get any larger, or
               ; smaller by more than a factor of 20 ???
;            mfac = ((ymult[*,iobj,i2])[igood] > 0.05) < 1
;            invsig = 0 * allivar[*,iobj,i1]
;            if (ct GT 0) then $
;             invsig[igood] = 1. / sqrt( 1./(allivar[*,iobj,i1])[igood] $
;              + mfac^2 / ((allivar[*,iobj,i2])[igood] ))
;            yfit1 = spfluxcorr_solve( loglam[*,iobj,*], $
;             allflux[*,iobj,i2], allflux[*,iobj,i1], invsig, $
;             mask=qgood, npoly=3, nback=1, $
;             sigvec=sigvec1, ymult=ymult1, yadd=yadd1, $
;             debug=(iiter EQ maxiter2) )
;            sigvec = sigvec > sigvec1
         endelse

         if (keyword_set(ymult1)) then begin
            ymult[*,iobj,i2] = ymult1
            yadd[*,iobj,i2] = yadd1
         endif
      endfor

   endfor

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
