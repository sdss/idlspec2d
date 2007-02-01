;+
; NAME:
;   spflux_fittemplates
;
; PURPOSE:
;   Fit atmospheric feature templates to a single exposure.
;
; CALLING SEQUENCE:
;   spflux_fittemplates, camera, loglam, flux, ivar, spline, helio_rv, $
;                        coeffs=coeffs, yfit=yfit, bkpts=bkpts,
;                        airmass=airmass, $
;                        doplot=doplot, xrange=xrange, yrange=yrange, title=title, $
;                        giffile=giffile
;
; INPUTS:
;   camera     - Name of the camera this exposure comes from ('r1',
;                etc.)
;   loglam     - wavelengths of standard star spectra [NWAVE]
;   flux       - normalized flux [NWAVE]
;   ivar       - normalized ivars [NWAVE]
;   spline     - bspline structure of flux calibration vector with
;                telluric features masked out 
;   helio_rv   - heliocentric-geocentric RV offset for the exposure.
;
; OPTIONAL INPUTS:
;   airmass         - optional airmass term, per pixel.
;   doplot          - Generate diagnostic plots.
;   xrange, yrange  - limit diagnostic plots to these ranges.
;   title           - title for diagnostic plot.
;   giffile         - name of multiframe GIF file, which has already
;                     been created.
;
; OPTIONAL OUTPUTS:
;   yfit       - The spline and all the fit templates, evaluated at
;                loglam.
;   coeffs     - The fit coefficients.
;   bkpts      - The base spline breakpoints.
;
; DATA FILES:
;   examples/spFluxtemplates-r1.fits, etc., which contain the templates.
;
; INTERNAL PROCEDURES:
;   fluxtempl_fit - an mpfitfun compatible fitting function.
;
; BUGS:
;   Needs diagnostics something awful.
;   All water templates should be fit together with one scale.
;   Give up on the optical depth fit; switch to PCAs.
;
function fluxtempl_fit, x, p, templates=templates, types=types

    nt = n_elements(templates[0, *])
    scales = p
                                ; Sum up the scaled templates
    sum = x * 0 + 1
    for i=0,nt-1 do begin 
       t = templates[*, i] - 1
       if types[i] EQ 2 then begin 
          tmin = min(t)         ; avoid the embarrassing (-tiny)^y
          t -= tmin    
          tsum = t^scales[i]
          tsum += tmin
       end else begin
          tsum = scales[i] * t
       end

       sum += tsum
    end

    return, sum
 end

pro spflux_fittemplates, camera, loglam, flux, ivar, spline, helio_rv, airmass=airmass, $
                         coeffs=coeffs, yfit=yfit, bkpts=bkpts, $
                         doplot=doplot, xrange=xrange, yrange=yrange, title=title, $
                         giffile=giffile

    config_dir = filepath('', root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
    templfiles = unixfind(config_dir, string('spFluxTemplates-', camera, '.fits'))

    templfile = templfiles[0]
    if (templfile eq '') then message, "could not find template file for camera ", camera

                                ; Evaluate the templates at the flux
                                ; wavelengths. Start with the base
                                ; spline.
    wsort = sort(loglam)
    sloglam = loglam[wsort]
    sflux = flux[wsort]
    sivar = ivar[wsort]
    sairmass = keyword_set(airmass) ? airmass[wsort] : 0
    basefit = bspline_valu(sloglam, spline, x2=sairmass)
    fullfit = basefit
    coeffs = [[1d, 0d]]

    fullspan = fltarr(n_elements(wsort)) + 1

                                ; Read all the templates, construct
                                ; full-wavelength vectors.
    xxx = mrdfits(templfile, 0, hdr)
    ntemplates = sxpar(hdr, "NTEMPLAT")
    keepspans = 0
    templivar = sivar * 0
    for i=0,ntemplates-1 do begin
       tset = mrdfits(templfile, i*2 + 1) 
       if not keyword_set(tset) then continue ; ??? - CPL

                                ; Shift the templates to barycentric
                                ; wavelengths. This mucks with the
                                ; tset internals, which is pretty
                                ; rude. But awful easy.
       tlam = 10^tset.fullbkpt
       tlam = tlam / (1 + helio_rv/299792.458)
       tset.fullbkpt[*] = alog10(tlam)
       splog, "HELIO_RV correction of ", helio_rv, " applied to templates."

       tspans = mrdfits(templfile, i*2 + 2, thdr)
       type = 1+sxpar(thdr, 'TEMPTYPE')
       nspans = (size(tspans,/dim))[0]

                                ; Put each span in as its own
                                ; template. We might be able to do
                                ; better by asserting that, say, all
                                ; the water features scale together.
       
       for t=0L, nspans-1 do begin
          wlimits = [tspans[t,0], tspans[t,1]]
          wmin = min(wlimits)
          wmax = max(wlimits)
          tpoints = where(sloglam GE wmin and sloglam LE wmax, npts)
          if npts eq 0 then $
             continue 
                                ; Note that while the splines are
                                ; defined outside the span, the values
                                ; would be completely nuts.
          templ = bspline_valu(sloglam[tpoints], tset)
          templivar[tpoints] = sivar[tpoints]
          fulltempl = fullspan
          fulltempl[tpoints] = templ

          alltemplates = keyword_set(alltemplates) ? [[alltemplates], [fulltempl]] : fulltempl
          alltypes = keyword_set(alltypes) ? [alltypes, type] : type
          keepspans += 1
       end

                                ; The redmost atmospheric region
                                ; extends right to the end of the
                                ; spectrum, and the last spline points
                                ; are suspect. Allow the last template
                                ; to tilt to correct this a little bit.
       if wmax GT alog10(9000) then begin
          templ = interpol([1.0,1.1],[wmin,wmax],sloglam(tpoints))

          fulltempl = fullspan
          fulltempl[tpoints] = templ

          alltemplates = [[alltemplates], [fulltempl]]
          alltypes = [alltypes, 1]
          keepspans += 1
       end
    end

                                ; We currently allow an exponential
                                ; optical depth type. If we got rid of
                                ; that we'd likely be better off
                                ; with PCA.
    fa = {templates:alltemplates, types:alltypes}
    p = replicate({value:0.0, fixed:0, limited:[0,0], limits:[0.0,0.0]}, keepspans)
    p.value = 1.0            ; template scale
    p.limited[0] = 1
                                ; iteratively fit the templates to the
                                ; whitened spectra, and the spline to
                                ; the templated spectra.
    maxiter = 50
    fittol = 1e-5
    coeffs = p.value * 0 - 999
    sfit = bspline_valu(sloglam, spline, x2=sairmass)
    for iter_i=0, maxiter-1 do begin

       wspec = sflux / sfit
       res = mpfitfun('fluxtempl_fit', sloglam, wspec, 1./sqrt(templivar), $
                      functargs=fa, yfit=tfit, parinfo=p, $
                      bestnorm=bestnorm, dof=dof, nfree=nfree, niter=niter, $
                      status=status, errmsg=errmsg, /quiet)
       p.value = res
       if status LE 0 then message, errmsg
       splog, 'iter: ', iter_i, ' fit: ', bestnorm/dof
       splog, 'scales: ', res
    
       outmask1 = 0
       x2 = 0

       tspec = sflux / tfit
       tspline = bspline_iterfit(sloglam, tspec, $
                                 invvar=sivar, lower=3, upper=3, nord=3, fullbkpt=bkpts, $
                                 maxrej=ceil(0.05*n_elements(XXXXXXX)), outmask=outmask1, $
                                 x2=x2, npoly=2*keyword_set(x2),requiren=1)
       if (max(tspline.coeff) EQ 0) then $
          message, 'B-spline fit failed!!'
       sfit = bspline_valu(sloglam, tspline)
       
       if max(abs(coeffs - res)) LT fittol then break 
       coeffs = res
    end

    if iter_i GE maxiter then $
       splog, "template fit did not converge! coeffs=", coeffs $
    else $
       splog, "template niter, coeffs:", iter_i, coeffs 
   

    fullfit = sfit * tfit
    yfit = fullfit
    resid = (sflux-fullfit) / fullfit

    if keyword_set(doplot) then begin
       wmin = min(loglam)
       wmax = max(loglam)

       stellarmask = spflux_masklines(loglam,/stellar)

       if not keyword_set(xrange) then xrange=10^[min(loglam),max(loglam)]
       if not keyword_set(yrange) then yrange=[0,max(fullfit)*1.25]

       resi = wmin + findgen(1000) / 999 * (wmax-wmin)
       resset = bspline_iterfit(loglam[wsort], resid, bkpt=resi)
       resspl = bspline_valu(resi, resset)

       if not keyword_set(title) then title=''
       djs_plot,10^loglam[wsort],flux[wsort],psym=3,yrange=yrange,xrange=xrange,title=title,xstyle=1
       djs_oplot,10^loglam[wsort],basefit,psym=3,color='yellow'
       djs_oplot,10^loglam[wsort],fullfit,psym=3,color='red'
       if keyword_set(bkpts) then $
          djs_oplot,10^bkpts,bspline_valu(bkpts, spline, x2=sairmass),psym=1,color='cyan'
       
                                ; Rail any offscale points so we can
                                ; see them.
       resid = resid < 0.195
       resid = resid > (-0.195)
       djs_plot,10^loglam[wsort],resid,psym=3,$
                yrange=[-0.2,0.2],xrange=xrange,xstyle=1,ystyle=1
       djs_oplot,10^loglam[wsort],0*loglam,psym=3,color='red'
       djs_oplot,10^resi,resspl, color='gray'
       if keyword_set(bkpts) then $
          djs_oplot,10^bkpts,bkpts*0,psym=1,color='cyan'

       if keyword_set(giffile) then begin
          tvlct,r,g,b,/get
          write_gif, giffile, tvrd(), r, g, b, /multiple
          print, "put GIF frame...."
       end
    end
 end
