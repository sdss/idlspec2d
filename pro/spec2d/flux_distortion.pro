;+
; NAME:
;   flux_distortion
;
; PURPOSE:
;   Compute the flux-distortion image for an entire plate.
;
; CALLING SEQUENCE:
;   corrimg = flux_distortion(objflux, objivar, andmask, ormask, plugmap=, $
;    loglam=, [ minflux=, minobj=, platefile=, plotfile=, hdr=, coeff= ] )
;
; INPUTS:
;   objflux    - Fluxes [NPIX,NFIBER]
;   objivar    - Inverse variances [NPIX,NFIBER]
;   andmask    - AND pixel mask [NPIX,NFIBER]
;   ormask     - OR pixel mask [NPIX,NFIBER]
;   plugmap    - Plug-map structure [NFIBER], where CALIBFLUX,CALIBFLUX_IVAR
;                will be used in preference over MAG for the photometric fluxes
;   loglam     - Wavelength vector in log10(Ang), which must be the same
;                for all objects [NPIX]
;
; OPTIONAL INPUTS:
;   minflux    - Minimum flux levels for objects to be used in the fit;
;                default to 5 nMgy in all three (gri) bands, corresponding
;                to 20.7-th mag.
;   minobj     - Minimum number of objects that have good fluxes in all
;                three gri-bands for computing the corrections; default to 50;
;                if fewer than this many, then CORRIMG is returned with all 1's
;   platefile  - If set, then read OBJFLUX and all the other inputs from
;                this spPlate file instead of using those inputs;
;                also, generate PostScript plots using the PLATESN procedure.
;   plotfile   - If set, then make a contour plot of the distortion
;                corrections to this PostScript file
;   hdr        - If set, then get the PLATE and MJD from this FITS header
;
; OUTPUTS:
;   corrimg    - Flux-distortion image by which OBJFLUX should be multiplied
;                [NPIX,NFIBER]
;
; OPTIONAL OUTPUTS:
;   coeff      - Best-fit coefficients for the distortion terms
;
; COMMENTS:
;   The the correction vectors are parameterized in terms of magnitude
;   (i.e. log-flux) that are achromatic with x, y, x^2, y^2, x*y,
;   where those are linear coordinates XFOCAL,YFOCAL from the plug-map.
;   There are also chromatic terms that scale as 1-(5070/wavelength)^2,
;   since that function gives an equal effect between 3900 and 5070 Ang
;   as between 5070 ang 9000 Ang.
;   There are also magnitude offsets as a function of spectrograph ID,
;   and a chromatic offset as a function of spectrograph ID.
;
;   In detail, the formulae is as follows (with 14 terms):
;     NEWFLUX = FLUX * [1 + a0*(SPECID EQ 1) + a1*(SPECID EQ 2)]
;               * exp{ a2*x + a3*y + a4*x*y + a5*x^2 + a6*y^2
;                + a7*x*LL + a8*y*LL 
;                + a9*LL*(SPECID EQ 1) + a10*LL*(SPECID EQ 2) }
;                + a11*x*y*LL + a12*x^2*LL + a13*y^2*LL
;   where x=XFOCAL, y=YFOCAL, LL = 1 - (5100 Ang/wavelength)^2
;
; EXAMPLES:
;
; BUGS:
;   At the moment, I've turned off the quadratic terms in the exponent
;   that also have LL, such as exp(x^2*LL).  That leaves 11 terms being used.
;
; PROCEDURES CALLED:
;   airtovac
;   djs_maskinterp()
;   djs_reject()
;   flux_distort_corrvec()
;   flux_distort_fn()
;   linterp
;   mpfit
;   platesn
;   readcol
;   splog
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   17-Feb-2004  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
forward_function mpfit, flux_distort_fn

;------------------------------------------------------------------------------
function flux_distort_corrvec, coeff, wavevec, thisplug

   lam0 = 5070.d0
   lwave2 = 1. - (lam0/wavevec)^2 ; Basically normalized to [0.5,2.0]
   xx = thisplug.xfocal / 320.d0
   yy = thisplug.yfocal / 320.d0
   specid = thisplug.spectrographid
   cvec = (1. + coeff[0] * (specid EQ 1) + coeff[1] * (specid EQ 2)) $
    * exp(coeff[2] * xx $
        + coeff[3] * yy $
        + coeff[4] * xx * yy $
        + coeff[5] * xx^2 $
        + coeff[6] * yy^2 $
        + coeff[7] * xx * lwave2 $
        + coeff[8] * yy * lwave2 $
        + coeff[9] * lwave2 * (specid EQ 1) $ 
        + coeff[10] * lwave2 * (specid EQ 2) $
        )
;        + coeff[11] * xx*yy * lwave2 $
;        + coeff[12] * xx^2 * lwave2 $
;        + coeff[13] * yy^2 * lwave2)

   return, cvec
end
;------------------------------------------------------------------------------
; Return a vector of all chi values.
function flux_distort_fn, coeff

   common com_flux_distort, trimflux, wavevec, fmask, calibflux, calibisig, $
    trimplug, outmask

   nobj = n_elements(trimplug)
   newflux = fltarr(nobj,3)
   for i=0L, nobj-1 do $
    for j=0, 2 do $
     newflux[i,j] = total(trimflux[*,i] * fmask[*,j] $
      * flux_distort_corrvec(coeff, wavevec, trimplug[i]), /double)

;   retval = (newflux / calibflux) - 1.
   retval = (newflux - calibflux) * calibisig
   if (keyword_set(outmask)) then for j=0, 2 do $
    retval[*,j] = retval[*,j] * outmask

   ; Re-caste from an array to a vector, since MPFIT needs a vector.
   return, retval[*]
end
;------------------------------------------------------------------------------
function flux_distortion, objflux, objivar, andmask, ormask, plugmap=plugmap, $
 loglam=loglam, minflux=minflux, minobj=minobj, platefile=platefile, $
 plotfile=plotfile, hdr=hdr, coeff=coeff

   common com_flux_distort, trimflux, wavevec, fmask, calibflux, calibisig, $
    trimplug, outmask

   if (NOT keyword_set(minobj)) then minobj = 50
   if (NOT keyword_set(minflux)) then minflux = 5.
   t0 = systime(1)

   if (keyword_set(platefile)) then begin
      objflux = mrdfits(platefile,0,hdr)
      objivar = mrdfits(platefile,1)
      andmask = mrdfits(platefile,2)
      ormask = mrdfits(platefile,3)
      plugmap = mrdfits(platefile,5)
      loglam = sxpar(hdr, 'COEFF0') $
       + dindgen(sxpar(hdr, 'NAXIS1')) * sxpar(hdr, 'COEFF1')
   endif

   dims = size(objflux, /dimens)
   npixobj = dims[0]
   nobj = dims[1]
   thisivar = skymask(objivar, andmask, ormask)

   wavevec = 10d0^loglam

   ;----------
   ; Read the three filter curves of interest

   flambda2fnu = wavevec^2 / 2.99792e18
   fmask = fltarr(npixobj, 3)
   ffiles = 'sdss_jun2001_' + ['g','r','i'] + '_atm.dat'
   for ifile=0, 2 do begin
      thisfile = filepath(ffiles[ifile], $
       root_dir=getenv('IDLUTILS_DIR'), subdirectory=['data','filters'])
      readcol, thisfile, fwave, fthru, /silent
      airtovac, fwave
      linterp, fwave, fthru, wavevec, fmask1
      fmask1 = fmask1 / total(fmask1)
      fmask[*,ifile] = fmask1 * flambda2fnu * 10^((22.5 + 48.6 - 2.5*17.)/2.5)
   endfor

   ;----------
   ; Deal with any old PLUGMAP structures without the linear flux units.

   if (tag_exist(plugmap,'CALIBFLUX') EQ 0) then begin
      plugmap = struct_addtags(plugmap, $
       replicate({calibflux:fltarr(5)},nobj))
      plugmap.calibflux = 10.^((22.5 - plugmap.mag) / 2.5)
   endif
   if (tag_exist(plugmap,'CALIBFLUX_IVAR') EQ 0) then begin
      plugmap = struct_addtags(plugmap, $
       replicate({calibflux_ivar:fltarr(5)},nobj))
   endif

   ;----------
   ; Select only objects that are not SKY or QSO, and
   ; with a positive flux in gri-bands,
   ; and at least 90% of pixels are good within the wavelengths of interest.
   ; Also, trim to only objects with known errors (according to CALIBFLUX_IVAR)
   ; if that's at least 80% of the otherwise-good objects.

   indx = where(wavevec GT 4000. AND wavevec LE 8300., nthis)
   fracgood = total(thisivar[indx,*] GT 0, 1) / nthis
   qivar = plugmap.calibflux_ivar[2] GT 0
   qtrim = strmatch(plugmap.objtype,'SKY*') EQ 0 $
    AND strmatch(plugmap.objtype,'QSO*') EQ 0 $
    AND plugmap.calibflux[1] GT minflux $
    AND plugmap.calibflux[2] GT minflux $
    AND plugmap.calibflux[3] GT minflux $
    AND fracgood GT 0.90
   if (total(qtrim AND qivar) GT 0.8*total(qtrim)) then begin
      splog, 'Trimming to ', 100*total(qtrim AND qivar)/total(qtrim), $
       '% of objects w/known photom errors from calibObj'
      qtrim = qtrim * qivar
   endif else begin
      splog,' No objects w/known photom errors from calibObj'
   endelse
   itrim = where(qtrim, ntrim)
   splog, 'Number of objects for fitting distortions = ', ntrim
   if (ntrim LT minobj) then begin
      splog, 'WARNING: Too few objects for fitting flux distortions ', ntrim
      return, 0 * objflux + 1
   endif
   calibflux = transpose(plugmap[itrim].calibflux[1:3])
   calibisig = sqrt(transpose(plugmap[itrim].calibflux_ivar[1:3]))

   ;----------
   ; Assign appropriate errors to all these points: either add 3% errors to
   ; objects that already have errors, or assign an error of 5%.

   qerr = calibisig GT 0
   i = where(qerr, ct)
   if (ct GT 0) then $
    calibisig[i] = 1. / sqrt(1./calibisig[i]^2 + (0.03 * calibflux[i])^2)
   i = where(qerr EQ 0, ct)
   if (ct GT 0) then $
    calibisig[i] = 1. / (0.05 * calibflux[i])

   trimflux = djs_maskinterp(objflux[*,itrim], thisivar[*,itrim] EQ 0, $
    iaxis=0, /const)
   trimplug = plugmap[itrim]

   ;----------
   ; Iterate the fit, rejecting outlier points.

   maxiter1 = 5
   maxiter2 = 50
   sigrej = 5.
   maxrej = ceil(0.05 * ntrim) ; Do not reject more than 5% of remaining objects
   npar = 11

   parinfo = {value: 0.d0, fixed: 0, limited: [0b,0b], limits: [0.d0,0.d0]}
   parinfo = replicate(parinfo, npar)
   parinfo.value = 1.d-3
   ftol = 1d-20
   gtol = 1d-20
   xtol = 1d-20
   iiter = 0L
   outmask = 0
   while (iiter LT maxiter1) do begin
      splog, 'Flux distortion fit iteration #', iiter
      splog, 'Initial chi^2=', total((flux_distort_fn(parinfo.value*0))^2)
      coeff = mpfit('flux_distort_fn', parinfo=parinfo, $
       maxiter=maxiter2, ftol=ftol, gtol=gtol, xtol=xtol, $
       niter=niter, status=status, /quiet)
      splog, 'MPFIT niter=', niter, ' status=', status

      ; For each object, set CHIVEC equal to the worst value in the 3 filters
      chiarr = abs(reform(flux_distort_fn(coeff), ntrim, 3))
      splog, 'Iteration #', iiter, ' chi^2=', total(chiarr^2)
      chivec = (chiarr[*,0] > chiarr[*,1]) > chiarr[*,2]

      ; Reject points w/out using the errors, but rather by simply compute
      ; the standard deviation of the results; this is the default behaviour
      ; for DJS_REJECT().
      if (djs_reject(chivec, 0*chivec, invvar=0*chivec+1, outmask=outmask, $
;      if (djs_reject(chivec, 0*chivec, outmask=outmask, $
       lower=sigrej, upper=sigrej, maxrej=maxrej)) $
       then iiter = maxiter1 $
      else $
       iiter = iiter + 1
      splog, 'Number of rejected objects = ', long(total(1-outmask))

      ; For the next iteration, start with the last best fit.
      parinfo.value = coeff
   endwhile

   for i=0, npar-1 do $
    splog, 'Parameter #', i, ' = ', coeff[i]

   corrimg = fltarr(npixobj, nobj)
   for i=0L, nobj-1 do $
    corrimg[*,i] = flux_distort_corrvec(coeff, wavevec, plugmap[i])

   minval = min(corrimg, max=maxval)
   sigval = stdev(corrimg)
   splog, 'Flux distortion min/max/sig = ', minval, maxval, sigval

   if (keyword_set(plotfile)) then begin
      xx = djs_laxisgen([641,641],iaxis=0) - 320
      yy = djs_laxisgen([641,641],iaxis=1) - 320
      thisplug = replicate(create_struct('XFOCAL', 0., 'YFOCAL', 0., $
       'SPECTROGRAPHID', 0), 641, 641)
      thisplug.xfocal = xx
      thisplug.yfocal = yy
      thisplug.spectrographid = 1 + (yy GT 0)

      dfpsplot, plotfile, /color

      pmulti = !p.multi
      ymargin = !y.margin
      yomargin = !y.omargin

      !p.multi = [0,2,3]
      !y.margin = [1,0]
      !y.omargin = [5,3]

      wcen = [4000., 5100., 9000.]

      for iplot=0, 2 do begin
         ; First skip the left-side plot...
         plot, [0,1], [0,1], /nodata, xstyle=5, ystyle=5
         if (keyword_set(hdr) AND iplot EQ 0) then begin
            platestr = string(sxpar(hdr,'PLATEID'), format='(i4)')
            mjdstr = string(sxpar(hdr,'MJD'), format='(i5)')
            xyouts, 1.0, 1.0, align=0.5, charsize=2, $
             'Flux Distortions PLATE=' + platestr + ' MJD=' + mjdstr
         endif

         if (iplot EQ 2) then begin
            xtitle = 'X [mm]'
            xtickname = ''
         endif else begin
            xtitle =''
            xtickname = strarr(20)+' '
         endelse

         cvec = flux_distort_corrvec(coeff, wcen[iplot], thisplug)
         mvec = -2.5*alog10(cvec) * (sqrt(xx^2 + yy^2) LT 320.)

         maxdiff = max(abs(mvec))
         if (maxdiff LT 0.10) then levels = 0.01*findgen(21) - 0.10 $
          else if (maxdiff LT 0.25) then levels = 0.025*findgen(21) - 0.25 $
          else if (maxdiff LT 0.50) then levels = 0.05*findgen(21) - 0.50 $
          else levels = 0.10*findgen(21) - 1.0
         c_colors = (levels GE 0) * djs_icolor('blue') $
          + (levels LT 0) * djs_icolor('red')
         contour, mvec, xx[*,0], transpose(yy[0,*]), /follow, $
          levels=levels, xrange=[-320,320], yrange=[-320,320], $
          /xstyle, /ystyle, c_colors=c_colors, c_charsize=1, $
          xtitle=xtitle, ytitle='Y [mm]', charsize=1, $
          title=title
         xyouts, -300, 280, string(wcen[iplot],format='(i4)') + ' Ang'
      endfor

      !p.multi = pmulti
      !y.margin = ymargin
      !y.omargin = yomargin
   endif

   if (keyword_set(platefile)) then begin
      platesn, objflux, objivar, $
       andmask, plugmap, loglam, hdr=hdr, plotfile='test1.ps'
      platesn, objflux*corrimg, objivar/corrimg^2, $
       andmask, plugmap, loglam, hdr=hdr, plotfile='test2.ps'
   endif

   splog, 'Time to compute distortions = ', systime(1)-t0, ' sec'

   return, corrimg
end
;------------------------------------------------------------------------------
