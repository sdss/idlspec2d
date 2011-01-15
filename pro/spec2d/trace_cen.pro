;+
; NAME:
;   trace_cen
;
; PURPOSE:
;   Find the fiber positions for the central row of an image.
;
; CALLING SEQUENCE:
;   xfiber = trace_cen( fimage, [ xstart=, ystart=, nmed=, nfiber=, nbundle=, $
;    fiberspace=, bundlespace=, xgood=, flux=, plottitle= ] )
;
; INPUTS:
;   fimage     - Image
;
; OPTIONAL INPUTS:
;   xstart     - Initial guess for X position of first fiber; default to 0
;   ystart     - Y position in image to search for initial X centers; default
;                to the central row
;   nmed       - Number of rows to median filter around YSTART;
;                default to 21
;   nfiber     - Number of fibers; default to 320
;   nbundle    - Number of bundles; default to 16
;   fiberspace - Fiber-to-fiber spacing in pix [NBUNDLE]
;   bundlespace- Extra spacing between bundles in pix [NBUNDLE-1]
;   plottitle  - If set, then create plot with this title
;
; OUTPUTS:
;   xfiber     - Vector of 320 X centers
;
; OPTIONAL OUTPUTS:
;   xgood      - Set to 1 for fibers that were actually found, 0 otherwise
;   flux       - Flux of each fiber
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_iterstat
;   mpfit
;
; REVISION HISTORY:
;   29-Aug-2009  Written by David Schlegel, LBL
;-
;------------------------------------------------------------------------------
; Cross-correlation routine with fractional shifts.
; Identical to C_CORRELATE if LAGS are all integer-valued.
function trace_correlate, yvec, fluxvec, lags
   nlag = n_elements(lags)
   npix = n_elements(yvec)
   cc = fltarr(nlag)
   ftmp = fluxvec - mean(fluxvec)
   ytmp = yvec - mean(yvec)
   for i=0L, nlag-1L do begin
      ytmp_shift = sshift(ytmp,lags[i])
      cc[i] = total(ytmp_shift * ftmp)
   endfor
   return, cc / (stddev(ytmp) * stddev(ftmp)) / npix
end
;------------------------------------------------------------------------------
function trace_param_to_xcen, param, name=name, $
 nfiber=nfiber, nbundle=nbundle

   i = where(name EQ 'xstart')
   xstart = param[i[0]]
   i = where(name EQ 'fiberspace', nbundle)
   fiberspace = param[i]
   if (nbundle GT 1) then begin
      i = where(name EQ 'bundlespace')
      bundlespace = param[i]
   endif else bundlespace = 0

   nperbundle = nfiber / nbundle
   bundlenum = lindgen(nfiber) / nperbundle ; bundle # for each fiber

   xfiber = fltarr(nfiber) + xstart
   for ibundle=0, nbundle-1 do begin
      if (ibundle EQ 0) then x0 = xstart $
       else x0 = xfiber[ibundle*nperbundle-1] + bundlespace[ibundle-1] $
        + fiberspace[ibundle-1]
      xfiber[ibundle*nperbundle:(ibundle+1)*nperbundle-1] = x0 $
       + findgen(nperbundle) * fiberspace[ibundle]
   endfor

   return, xfiber
end
;------------------------------------------------------------------------------
function trace_param_to_vec, param, name=name, $
 nfiber=nfiber, nbundle=nbundle, npix=npix

   i = where(name EQ 'psfsigma')
   psfsigma = param[i]
   i = where(name EQ 'flux')
   flux = param[i]

   nperbundle = nfiber / nbundle
   bundlenum = lindgen(nfiber) / nperbundle ; bundle # for each fiber

   xfiber = trace_param_to_xcen(param, name=name, $
    nfiber=nfiber, nbundle=nbundle)

   sigfiber = psfsigma[bundlenum]

   xvec = findgen(npix)
   fmodel = fltarr(npix)
   for ifiber=0, nfiber-1 do $
    if (flux[ifiber] NE 0) then $
     fmodel += flux[ifiber]/(sigfiber[ifiber]*sqrt(2.*!pi)) $
      * exp(-0.5  * ((xvec - xfiber[ifiber])/sigfiber[ifiber])^2)

   return, fmodel
end
;------------------------------------------------------------------------------
function fn_trace_model, params, fluxvec=fluxvec, _EXTRA=KeywordsForParamtovec

   fmodel = trace_param_to_vec(params, _EXTRA=KeywordsForParamtovec)

   chivec = fluxvec - fmodel ; ???

   return, chivec
end
;------------------------------------------------------------------------------
function trace_cen, fimage, xstart=xstart, ystart=ystart, nmed=nmed, $
 nfiber=nfiber, nbundle=nbundle, fiberspace=fiberspace, bundlespace=bundlespace, $
 xgood=xgood, flux=flux, plottitle=plottitle, fluxvec=fluxvec, fmodel=fmodel

   ; Need 1 parameter
   if (n_params() LT 1) then begin
      return, -1
   endif

   dims = size(fimage, /dimens)
   nx = dims[0]
   ny = dims[1]
   if (NOT keyword_set(ystart)) then ystart = ny/2
   if (NOT keyword_set(nmed)) then nmed = 11

   xstart_tol = 50.
   xstart_step = 0.25
   bundlespace_tol = 0.25 * fiberspace[0]
   bundlespace_step = 0.25
   scale_range = [0.9900, 1.0100, 0.0005] ; range of scalings +/- 1.0%
   psfsigma = 1.0 ; pix

   ; Make a sub-image copy of the image and error map
   ylo = ystart - (nmed-1)/2 > 0
   yhi = ystart + (nmed-1)/2 < ny-1
   subimg = fimage[*,ylo:yhi]

   ; Median filter along each column of the subimage
   fluxvec = djs_median(subimg, 2)

   ;----------
   ; Parameters as follows:
   ;   xstart
   ;   bundlespace[nbundle-1]
   ;   fiberspace[nbundle]
   ;   psfsigma[nbundle]
   ;   flux[nfiber]

   nparam = 3 * nbundle + nfiber
   parinfo = replicate({name: '', value: 0.D, fixed: 0B, $
    limited:[0,0], limits:[0.D,0]}, nparam)
   parinfo[0].name = 'xstart'
   parinfo[0].value = xstart
   parinfo[1:nbundle-1].name = 'bundlespace'
   parinfo[1:nbundle-1].value = bundlespace
   parinfo[1:nbundle-1].limits[0] = bundlespace - bundlespace_tol
   parinfo[1:nbundle-1].limits[1] = bundlespace + bundlespace_tol
   parinfo[1:nbundle-1].limited = [1,1]
   parinfo[nbundle:2*nbundle-1].name = 'fiberspace'
   parinfo[nbundle:2*nbundle-1].value = fiberspace
   parinfo[2*nbundle:3*nbundle-1].name = 'psfsigma'
   parinfo[2*nbundle:3*nbundle-1].value = psfsigma
   parinfo[3*nbundle:3*nbundle+nfiber-1].name = 'flux'
   parinfo[3*nbundle:3*nbundle+nfiber-1].value = total(fluxvec)/nfiber
   parinfo[3*nbundle:3*nbundle+nfiber-1].limits = [0,1e6] ; positive fluxes
   parinfo[3*nbundle:3*nbundle+nfiber-1].limited = [1,0]

   functargs = {name: parinfo.name, nfiber: nfiber, nbundle: nbundle, npix: nx}

   ;----------
   ; Find the best xstart using a simple cross-correlation
   ; Explore a small range of possible scalings of the fiber spacings

   nlag = 2*xstart_tol / xstart_step
   lags = - xstart_tol + xstart_step * findgen(nlag)
   cc_best = -1
   scale_best = 1
   lag_best = 0
   for scale=scale_range[0], scale_range[1], scale_range[2] do begin
      parinfo[nbundle:2*nbundle-1].value = scale * fiberspace
      fmodel = trace_param_to_vec(parinfo.value, _EXTRA=functargs)
      cc = trace_correlate(fmodel, fluxvec, lags)
      thiscc = max(cc, imax)
      if (thiscc GE cc_best) then begin
         cc_best = thiscc
         scale_best = scale
         lag_best = lags[imax]
      endif
   endfor
   parinfo[nbundle:2*nbundle-1].value = scale_best * fiberspace
   parinfo[0].value += lag_best
   splog, 'scale = ', scale_best
   splog, 'xstart = ', parinfo[0].value
;fmodel = trace_param_to_vec(parinfo.value, _EXTRA=functargs)
;splot,fluxvec & soplot,fmodel,color='red'

   ;----------
   ; For each bundle, find the best bundlespace
   ; Cross-correlate with only a model of that one bundle

   nperbundle = nfiber / nbundle
   bundlenum = lindgen(nfiber) / nperbundle ; bundle # for each fiber
   for ibundle=0, nbundle-1 do begin
      nlag = 2*bundlespace_tol / bundlespace_step
      lags = -bundlespace_tol + bundlespace_step * findgen(nlag)
      parinfo1 = parinfo
      izero = where(bundlenum NE ibundle)
      parinfo1[3*nbundle+izero].value = 0 ; zero out fluxes of other fibers
      fmodel1 = trace_param_to_vec(parinfo1.value, _EXTRA=functargs)
      cc = trace_correlate(fmodel1, fluxvec, lags)
      junk = max(cc, imax)
      splog, 'Bundle ', ibundle, ' shifted ', lags[imax]
      parinfo[ibundle].value += lags[imax]

      ; Make sure these are still in-bounds (round-off can send them out)
; The below not working if values are outside limits...???
      if (parinfo[ibundle].limited[0] EQ 1) then $
       parinfo[ibundle].value = parinfo[ibundle].value $
        > parinfo[ibundle].limits[0]
      if (parinfo[ibundle].limited[1] EQ 1) then $
       parinfo[ibundle].value = parinfo[ibundle].value $
        < parinfo[ibundle].limits[1]
;foo = trace_param_to_vec(parinfo1.value, _EXTRA=functargs)
;splot,fluxvec,xr=[240,400] & soplot,foo,color='red'

; Below can take parameters out-of-bounds from LIMITS !!!???
;      ; Following line is to prevent shifting start location of other bundles
;      if (ibundle LT nbundle-1) then $
;       parinfo[ibundle+1:nbundle-1].value -= lags[imax]
   endfor

   ;----------
   ; Loop through each bundle, solving for all the parameters
   ; of that bundle simultaneously using MPFIT.

   functargs_all = create_struct(functargs, {fluxvec: fluxvec})
   for ibundle=0, nbundle-1 do begin
      ; Create a fluxvec that subtracts the other bundle fits
      parinfo_sub = parinfo
      parinfo_sub[3*nbundle+where(bundlenum EQ ibundle)].value = 0
      fmodel_sub = trace_param_to_vec(parinfo_sub.value, _EXTRA=functargs)
      functargs_all = create_struct(functargs, {fluxvec: fluxvec-fmodel_sub})

      ; Set up to fix all parameters except for those for this bundle
      parinfo_new = parinfo
      parinfo_new[3*nbundle+where(bundlenum NE ibundle)].value = 0
      parinfo_new.fixed = 1B
      parinfo_new[ibundle].fixed = 0B ; Fit for the offset of the whole bundle
      parinfo_new[nbundle+ibundle].fixed = 0B ; Fit for the fiber spacing
      parinfo_new[2*nbundle+ibundle].fixed = 0B ; Fit for the PSF width
      ; Below fits for flux of these fibers
      parinfo_new[3*nbundle+ibundle*nperbundle+lindgen(nperbundle)].fixed = 0B
      niter = 0
      nfev = 0
      acoeff = mpfit('fn_trace_model', parinfo=parinfo_new, $
       functargs=functargs_all, $
       maxiter=100, niter=niter, nfev=nfev, status=status, /quiet, $
       errmsg=errmsg)
      if (status EQ 0) then $
       message, 'Invalid arguments to MPFIT: '+errmsg

;foo = trace_param_to_vec(parinfo_new.value, _EXTRA=functargs_all)
;foo2 = trace_param_to_vec(acoeff, _EXTRA=functargs_all)
;splot,fluxvec-fmodel_sub & soplot,foo2,color='red'

      ; Following line is to prevent shifting start location of other bundles
;      if (ibundle LT nbundle-1) then $
;       parinfo[ibundle+1:nbundle-1].value += $
;        (parinfo[ibundle].value - acoeff[ibundle])

      j = where(parinfo_new.fixed EQ 0)
      parinfo[j].value = acoeff[j] ; Copy over all the re-fit parameters
      splog, 'Bundle ', ibundle, ' bundlegap=', parinfo[ibundle].value, $
       ' niter=', niter
   endfor

   flux = parinfo[3*nbundle:3*nbundle+nfiber-1].value
   djs_iterstat, flux, sigrej=3.0, mean=mn
   xgood = flux GE 0.5*mn

   xfiber = trace_param_to_xcen(parinfo.value, _EXTRA=functargs)

   if (arg_present(fmodel) OR keyword_set(plottitle)) then $
    fmodel = trace_param_to_vec(parinfo.value, _EXTRA=functargs)

   if (keyword_set(plottitle)) then begin
      nplotrow = 5
      !p.multi = [0,1,nplotrow]
      xvec = findgen(nx)
      fmax = max(fmodel)
      yrange = [0,1.1*fmax]
      for j=0, nplotrow-1 do begin
         xrange=nx*[j,j+1]/float(nplotrow)
         if (j EQ 0) then title=plottitle $
          else title=''
         indx = where(xvec GE xrange[0] AND xvec LE xrange[1])
         djs_plot, xvec[indx], fluxvec[indx], $
          xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
          ymargin=[0.5,0.5]+[-0.4, 0.4]*(nplotrow-1-j), xmargin=[1, 1], $
          ticklen=0, title=title, charsize=2.5, $
          xtickname=replicate(' ',10), ytickname=replicate(' ',10)
         djs_oplot, xvec[indx], fmodel[indx], color='green'
         ibad = where(xfiber GE xrange[0] AND xfiber LE xrange[1] $
          AND xgood EQ 0, nbad)
         if (nbad GT 0) then begin
            djs_oplot, xfiber[ibad], 0*ibad+0.05*fmax, psym=2, color='red'
            for k=0, nbad-1 do $
             djs_xyouts, xfiber[ibad[k]], fmax, strtrim(ibad[k]+1,2) $
              + (k MOD 2 ? '' : '    '), $
              align=1.0, color='red', orient=90
         endif
      endfor
      !p.multi = [0,1,1]
   endif

   return, xfiber
end
;------------------------------------------------------------------------------
