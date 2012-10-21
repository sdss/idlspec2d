;+
; NAME:
;   flux_distortion
;
; PURPOSE:
;   Compute the flux-distortion image for an entire plate.
;
; CALLING SEQUENCE:
;   corrimg = flux_distortion(objflux, objivar, andmask, ormask, plugmap=, $
;    loglam=, [ minflux=, minobj=, maxdelta=, platefile=, plotfile=, hdr=, $
;    coeff= ] )
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
;   minflux    - Minimum flux levels in ugriz-bands for objects to be used
;                in the fit; set to zero to ignore a band; default to
;                [0,0,2.5,0,0] to limit to r<21.5 objects
;                default to 5 nMgy in all three (gri) bands, corresponding
;                to 20.7-th mag.
;   minobj     - Minimum number of objects that have good fluxes in all
;                three gri-bands for computing the corrections; default to 50;
;                if fewer than this many, then CORRIMG is returned with all 1's
;   maxdelta   - Maximum peak deviation in flux distortion image allowed
;                per iteration from a change in any one parameter; default
;                to 0.03
;   platefile  - If set, then read OBJFLUX and all the other inputs from
;                this spPlate file instead of using those inputs;
;                also, generate PostScript plots using the PLATESN procedure.
;   plotfile   - If set, then make a contour plot of the distortion
;                corrections to this PostScript file
;   hdr        - If set, then get the PLATE and MJD from this FITS
;                header
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
;                + a11*LL^2*(SPECID EQ 1) + a12*LL^2*(SPECID EQ 2) }
;   where x=XFOCAL, y=YFOCAL, LL = 1 - (5100 Ang/wavelength)^2
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   airtovac
;   dfpsclose
;   dfpsplot
;   djs_maskinterp()
;   djs_oplot
;   djs_plot
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
   xx = thisplug.xfocal / 325.d0 ; This 325 is mm on the plate
   yy = thisplug.yfocal / 325.d0 ; This 325 is mm on the plate
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
        + coeff[11] * lwave2^2 * (specid EQ 1) $ 
        + coeff[12] * lwave2^2 * (specid EQ 2) $
        )
;        + coeff[11] * xx*yy * lwave2 $
;        + coeff[12] * xx^2 * lwave2 $
;        + coeff[13] * yy^2 * lwave2)

   return, cvec
end
;------------------------------------------------------------------------------
; Return a vector of all chi values.
function flux_distort_fn, coeff, nomask=nomask

   common com_flux_distort, trimflux, wavevec, fmask, calibflux, calibisig, $
    trimplug, outmask

   nobj = n_elements(trimplug)
   newflux = fltarr(nobj,3)
   for i=0L, nobj-1 do $
    for j=0, 2 do $
     newflux[i,j] = total(trimflux[*,i] * fmask[*,j] $
      * flux_distort_corrvec(coeff, wavevec, trimplug[i]))

;   retval = (newflux / calibflux) - 1.
   retval = (newflux - calibflux) * calibisig
   if (keyword_set(outmask) AND (keyword_set(nomask) EQ 0)) then $
    for j=0, 2 do $
     retval[*,j] = retval[*,j] * outmask

   ; Re-caste from an array to a vector, since MPFIT needs a vector.
   return, retval[*]
end
;------------------------------------------------------------------------------
function flux_distortion, objflux, objivar, andmask, ormask, plugmap=plugmap, $
 loglam=loglam, minflux=minflux1, minobj=minobj, maxdelta=maxdelta, $
 platefile=platefile, plotfile=plotfile, hdr=hdr, coeff=coeff

   common com_flux_distort, trimflux, wavevec, fmask, calibflux, calibisig, $
    trimplug, outmask

   if (NOT keyword_set(minobj)) then minobj = 10
   if (keyword_set(minflux1)) then minflux = minflux1 $
    else minflux = [0,0,2.5,0,0] ; Limit to r<21.5 only
   if (n_elements(maxdelta) EQ 0) then maxdelta = 0.03
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

   if (keyword_set(hdr)) then begin
      platestr = string(sxpar(hdr,'PLATEID'), format='(i4)')
      mjdstr = string(sxpar(hdr,'MJD'), format='(i5)')
      plottitle = 'PLATE=' + platestr + ' MJD=' + mjdstr
   endif

   dims = size(objflux, /dimens)
   npixobj = dims[0]
   nobj = dims[1]
   thisivar = skymask(objivar, andmask, ormask)

   wavevec = 10d0^loglam

   ;----------
   ; Set up for plots

   if (keyword_set(plotfile)) then begin
      dfpsplot, plotfile, /color
      pmulti = !p.multi
      ymargin = !y.margin
      yomargin = !y.omargin
      if (NOT keyword_set(plottitle)) then plottitle = ''
   endif

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
   ; Mask any flux values associated with catastrophically bad astrometry
   ; (ASB, 2012sep, 2012oct):
;   if tag_exist(plugmap, 'CALIB_STATUS') then plugmap.calibflux_ivar *= $
;      ((plugmap.calib_status AND sdss_flagval('CALIB_STATUS', 'ASTROMBAD')) eq 0)
   badflag = sdss_astrombad(plugmap.run, plugmap.camcol, plugmap.field)
   for iobj = 0L, nobj-1 do plugmap[iobj].calibflux_ivar[*] *= (badflag[iobj] eq 0)

   ;----------
   ; Explicitly select only SPECTROPHOTO_STD and REDDEN_STD stars
   ; on the first pointing offset (if there are multiple pointings),
   ; and with a positive flux in gri-bands,
   ; and at least 90% of pixels are good within the wavelengths of interest.
   ; Also, trim to only objects with known errors (according to CALIBFLUX_IVAR)
   ; if that's at least 80% of the otherwise-good objects.

   indx = where(wavevec GT 4000. AND wavevec LE 8300., nthis)
   fracgood = total(thisivar[indx,*] GT 0, 1) / nthis
   qivar = plugmap.calibflux_ivar[2] GT 0
   qtrim = (strmatch(plugmap.objtype,'SPECTROPHOTO_STD*') NE 0 $
    OR strmatch(plugmap.objtype,'REDDEN_STD*') NE 0) $
    AND plugmap.offsetid EQ 1 $
    AND fracgood GT 0.90
   for i=0, 4 do if (minflux[i] NE 0) then $
    qtrim *= (plugmap.calibflux[i] GT minflux[i])
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
   maxiter2 = 200
   sigrej = 5.
   maxrej = ceil(0.05 * ntrim) ; Do not reject more than 5% of remaining objects
   npar = 13

   parinfo = {value: 0.d0, fixed: 0, limited: [0b,0b], limits: [0.d0,0.d0], $
    mpmaxstep: 0.d0}
   parinfo = replicate(parinfo, npar)

   ; Compute the maximum step size per iteration for each parameter
   ; based upon MAXDPERSTEP.
   coeff_eps = 0.01
   corrimg = fltarr(npixobj, nobj)
   for ipar=0, npar-1 do begin
      coeff = dblarr(npar)
      coeff[ipar] = coeff_eps
      for i=0L, nobj-1 do $
       corrimg[*,i] = flux_distort_corrvec(coeff, wavevec, plugmap[i])
      thismin = min(corrimg, max=thismax)
      thisdiff = thismax - thismin
      parinfo[ipar].mpmaxstep = coeff_eps * maxdelta / thisdiff
   endfor

   ftol = 1d-20
   gtol = 1d-20
   xtol = 1d-20
   iiter = 0L
   outmask = 0

   while (iiter LT maxiter1) do begin
      splog, 'Flux distortion fit iteration #', iiter
      splog, 'Initial chi^2=', total((flux_distort_fn(parinfo.value*0))^2), $
       ' dof=', long(total(fmask NE 0) - n_elements(parinfo))

      ; Always start with tiny values, rather than values from the
      ; previous fit.  This is to prevent us from "walking away" from
      ; the correction solution in the case where some initially bad points
      ; effect the fit.
      parinfo.value = 0
      coeff = mpfit('flux_distort_fn', parinfo=parinfo, $
       maxiter=maxiter2, ftol=ftol, gtol=gtol, xtol=xtol, $
       niter=niter, status=status, /quiet)
      splog, 'MPFIT niter=', niter, ' status=', status

      ; Print the chi^2 after masking rejected points
      thischi2 = total((flux_distort_fn(coeff))^2)
      splog, 'Iteration #', iiter, ' chi^2=', thischi2

      ; For each object, set CHIVEC equal to the worst value in the 3 filters
      chiarr = abs(reform(flux_distort_fn(coeff, /nomask), ntrim, 3))
      chivec = (chiarr[*,0] > chiarr[*,1]) > chiarr[*,2]

      ; There was a bug in v5_0_0 of this code, before the /nomask option
      ; was implemented.  Rejected objects had chi values of zero, which
      ; meant that they would be put back in the fits.

      ; Reject points w/out using the errors, but rather by simply computing
      ; the standard deviation of the results; this is the default behaviour
      ; for DJS_REJECT().
;      if (djs_reject(chivec, 0*chivec, invvar=0*chivec+1, outmask=outmask, $
      ; Reject points using the errors.
      if (djs_reject(chivec, 0*chivec, outmask=outmask, $
       lower=sigrej, upper=sigrej, maxrej=maxrej)) $
       then iiter = maxiter1 $
      else $
       iiter = iiter + 1
      splog, 'Number of rejected objects = ', long(total(1-outmask))

      ; For the next iteration, start with the last best fit.
;      parinfo.value = coeff
   endwhile

   for i=0, npar-1 do $
    splog, 'Parameter #', i, ' = ', coeff[i]

   corrimg = fltarr(npixobj, nobj)
   for i=0L, nobj-1 do $
    corrimg[*,i] = flux_distort_corrvec(coeff, wavevec, plugmap[i])

   minval = min(corrimg, max=maxval)
   sigval = stddev(corrimg,/double)
   splog, 'Flux distortion min/max/sig = ', minval, maxval, sigval

   ;----------
   ; If multiple pointings on this plate, then deal with that here by
   ; rescaling the fluxes of each fiber according to its exposure time.

   igood = where(plugmap.fiberid GT 0)
   offsetlist = plugmap[igood].offsetid
   offsetlist = offsetlist[uniq(offsetlist,sort(offsetlist))]
   nlist = n_elements(offsetlist)
   if (nlist GT 1) then begin
      for ilist=1L, nlist-1L do begin
         indx1 = where(plugmap.offsetid EQ offsetlist[0] $
          AND strmatch(plugmap.objtype,'SKY*') EQ 0)
         indx2 = where(plugmap.offsetid EQ offsetlist[ilist] $
          AND strmatch(plugmap.objtype,'SKY*') EQ 0)
         spherematch, plugmap[indx1].ra, plugmap[indx1].dec, $
          plugmap[indx2].ra, plugmap[indx2].dec, 0.1/3600, i1, i2, d12
         nmatch = n_elements(i1) * (i1[0] NE -1)
         ratioarr = 0
         splog, nmatch, ' common objects between offset #', offsetlist[0], $
          ' and offset #', offsetlist[ilist]
         if (nmatch GT 0) then begin
            for imatch=0L, nmatch-1 do begin
               ; Solve for the ratio vector (YMULT) which is what we
               ; multiply the 2nd spectrum by to get this first.
               ; First, correct by the current flux-distortion image,
               ; which hardly matters if the spatial offsets are small.
               ; This is a quadratic (3 terms).
               minlog = min(loglam, max=maxlog)
               xvector = (loglam - minlog) / (maxlog - minlog)
               solve_poly_ratio, xvector, $
                objflux[*,indx2[i2[imatch]]] * corrimg[*,indx2[i2[imatch]]], $
                objflux[*,indx1[i1[imatch]]] * corrimg[*,indx1[i1[imatch]]], $
                thisivar[*,indx2[i2[imatch]]] / corrimg[*,indx2[i2[imatch]]]^2, $
                thisivar[*,indx1[i1[imatch]]] / corrimg[*,indx1[i1[imatch]]]^2, $
                npoly=3, nback=0, ymult=ymult
               ratioarr = keyword_set(ratioarr) ? [[ratioarr],[ymult]] $
                : ymult
            endfor
            if (keyword_set(ratioarr)) then begin
               if (size(ratioarr,/n_dimen) EQ 1) then nvec = 1 $
                else nvec = (size(ratioarr,/dimens))[1]
               if (nvec EQ 1) then ratiovec = ratioarr $
                else ratiovec = total(ratioarr,2) / nvec
;                else ratiovec = djs_median(ratioarr,2)
               splog, 'Scaling offset #', offsetlist[ilist], $
                ' by flux ratio vector with mean = ', $
                mean(ratiovec)
               splog, 'Scaling offset #', offsetlist[ilist], $
                ' RMS = ', 100.*stddev(ratioarr) / mean(ratiovec), ' %'
               for j=0L, n_elements(indx2)-1L do $
                corrimg[*,indx2[j]] = corrimg[*,indx2[j]] * ratiovec

               ; Make a plot
               if (keyword_set(plotfile)) then begin
                  !p.multi = [0,1,2]
                  xrange = minmax(10.^loglam)
                  djs_plot, xrange, minmax(ratioarr), $
                   /xstyle, /ynozero, /nodata, $
                   xtitle='Wavelength [Ang]', ytitle='Flux ratios', $
                   title = plottitle + ' Flux ratios for offset #' $
                    + strtrim(string(ilist),2)
                  for i=0L, nvec-1L do $
                   djs_oplot, 10.^loglam, ratioarr[*,i]
                  djs_oplot, 10.^loglam, ratiovec, color='red'
                  xyouts, total(!x.crange*[0.95,0.05]), $
                   total(!y.crange*[0.10,0.90]), $
                   strtrim(string(nvec),2) + ' non-SKY repeat objects'
                  djs_xyouts, total(!x.crange*[0.95,0.05]), $
                   total(!y.crange*[0.15,0.85]), 'Mean ratio vector', $
                   color='red'

                  plotratio = ratioarr $
                   / rebin(ratiovec,n_elements(ratiovec),nvec)
                  djs_plot, xrange, minmax(plotratio), $
                   /xstyle, /ynozero, /nodata, $
                   xtitle='Wavelength [Ang]', $
                   ytitle='Flux ratios / Mean vector'
                  for i=0L, nvec-1L do $
                   djs_oplot, 10.^loglam, plotratio[*,i]
                  djs_oplot, xrange, [1,1], color='red'
               endif
            endif
         endif
         if (nmatch EQ 0 OR keyword_set(ratioarr) EQ 0) then begin
            thisratio = plugmap[indx1[i1[0]]].sci_exptime $
             / plugmap[indx2[i2[0]]].sci_exptime
            splog, 'Scaling fluxes by exposure time ratio = ', thisratio
            corrimg[*,indx2] = corrimg[*,indx2] * thisratio
         endif
      endfor
   endif

   if (keyword_set(plotfile)) then begin
      ; Construct a fake plugmap across the image, with a sampling of
      ; every millimeter in x and y to evaluate these fits.
      xx = djs_laxisgen([651,651],iaxis=0) - 325
      yy = djs_laxisgen([651,651],iaxis=1) - 325
      thisplug = replicate(create_struct('XFOCAL', 0., 'YFOCAL', 0., $
       'SPECTROGRAPHID', 0), 651, 651)
      thisplug.xfocal = xx
      thisplug.yfocal = yy
      thisplug.spectrographid = 1 + (yy GT 0)

      !p.multi = [0,2,3]
      !y.margin = [1,0]
      !y.omargin = [5,3]

      wcen = [4000., 5100., 9000.]

      for iplot=0, 2 do begin
         ; First skip the left-side plot...
         plot, [0,1], [0,1], /nodata, xstyle=5, ystyle=5
         if (iplot EQ 0) then $
          xyouts, 1.0, 1.0, align=0.5, charsize=2, $
           plottitle + ' Flux Distortions'

         if (iplot EQ 2) then begin
            xtitle = 'X [mm]'
            xtickname = ''
         endif else begin
            xtitle =''
            xtickname = strarr(20)+' '
         endelse

         cvec = flux_distort_corrvec(coeff, wcen[iplot], thisplug)
         mvec = -2.5*alog10(cvec) * (sqrt(xx^2 + yy^2) LT 325.)

         maxdiff = max(abs(mvec))
         if (maxdiff LT 0.10) then levels = 0.01*findgen(21) - 0.10 $
          else if (maxdiff LT 0.25) then levels = 0.025*findgen(21) - 0.25 $
          else if (maxdiff LT 0.50) then levels = 0.05*findgen(21) - 0.50 $
          else levels = 0.10*findgen(21) - 1.0
         c_colors = (levels GE 0) * djs_icolor('blue') $
          + (levels LT 0) * djs_icolor('red')
         contour, mvec, xx[*,0], transpose(yy[0,*]), /follow, $
          levels=levels, xrange=[-325,325], yrange=[-325,325], $
          /xstyle, /ystyle, c_colors=c_colors, c_charsize=1, $
          xtitle=xtitle, ytitle='Y [mm]', charsize=1, $
          title=title
         xyouts, -300, 280, string(wcen[iplot],format='(i4)') + ' Ang'
      endfor

      !p.multi = pmulti
      !y.margin = ymargin
      !y.omargin = yomargin

      dfpsclose
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
