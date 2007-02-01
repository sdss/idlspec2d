;+
; NAME:
;   spflux_maketemplates
;
; PURPOSE:
;   Generate atmospheric band templates
;
; CALLING SEQUENCE:
;   spflux_maketemplates, camera, loglam, flux, ivar, spline, helio_rv, $
;                        coeffs=coeffs, yfit=yfit, bkpts=bkpts, $
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
;   examples/r1Fluxtemplates.fits, etc., which contain the templates.
;
; INTERNAL PROCEDURES:
;   mkteltemp_getstds - gathers all the pieces which we need.
;   fluxtempl_fit - an mpfitfun compatible fitting function.
;
; Create templates for telluric features 
;  . Start with the largest available data set; 
;  . Fit a low-order spline through the non-telluric regions, 
;  . Divide that spline out,
;  . Select the one-or-more telluric regions as templates.
;
; We are only fitting to a few types of features:
;  1. a handful of water bands. These appear to scale nicely together
;  and linearly: this should be changed to a PCA.
;  2. the two oxygen bands. 
;
; OUTPUT:
;  One FITS file, with one HDU for each type of template. That HDU
;  has:
;    . a TMPLTYPE card
;    . the bspline structure for _all_ the segments ("spans") in the
;    template.
;    . the wavelength bounds for each span in the template. The spline
;    is only valid within those ranges.
;
pro mkteltemp_getstds, cam, plates, loglam, flux, ivar, files
    rootdir = getenv('SPECTRO_DATA')
    if rootdir eq '' then message, 'SPECTRO_DATA not defined!'

    for p=0,n_elements(plates)-1 do begin
       fcfiles = unixfind(string(rootdir, plates[p], format='(%"%s/%04d")'), $
                          string('spFluxcalib-',cam,'-*'), /debug)

       for fc=0,n_elements(fcfiles)-1 do begin
          stdparts,fcfiles[fc], fibers,cals,rawspecs,rawivars,flats,flativars,superflats,$
                   models,ll,framefile,helio_rv
          if helio_rv eq -1000 then continue

          mr = rawspecs / models
          mrivar = rawivars * models^2
          fa = spflux_mratio_flatten(ll, mr, mrivar)
          netspec = mr/fa

                                ; Shift the spectra _back_ to
                                ; geocentric wavelengths:
          splog, "undoing heliocentric correction of ", helio_rv, " to ", framefile
          helio_rv = -helio_rv
          lambda = 10^ll
          lambda = lambda / (1 + helio_rv/299792.458) ; directly from fitvacset.
          ll = alog10(lambda)
                                ; 
          for f=0,n_elements(fibers)-1 do begin
             spec1 = netspec[*,f]
             ivar1 = mrivar[*,f]

                                ; Need to normalize somewhere, since
                                ; we intend to combine standards from
                                ; different exposures.
             djs_iterstat, spec1, invvar=ivar1, mean=div1
             spec1 /= div1
             ivar1 /= sqrt(div1)

             if (fc+f eq 0) then begin
                loglam = ll[*,f]
                flux = spec1
                ivar = ivar1
                files = [fcfiles[fc]]
             end else begin
                loglam = [[loglam], [ll[*,f]]]
                flux = [[flux], [spec1]]
                ivar = [[ivar], [ivar1]]
                files = [files, fcfiles[fc]]
             end
          end
       end
    end
 end


pro spflux_maketellurictemplates, cam, plates, filename, doplot=doplot, giffile=giffile, debug=debug

    color = strmid(cam, 0, 1)
                                ; n.b. this converts loglam back to
                                ; geocentric.
    mkteltemp_getstds, cam, plates, geoll, flux, ivar, files
    isort = sort(geoll)
    nspec = (size(flux,/dim))[1]

                                ; Flatten out the aggregate flux
                                ; vector, built without fitting to the
                                ; telluric regions.
    basespline = spflux_basespline(color, geoll, flux, ivar, $
                                   outmask=outmask, bkpts=bkpts)
    bseval = bspline_valu(geoll[isort], basespline)
    flatflux = flux[isort] / bseval
    
                                ; Now fit to the telluric features.
                                ; Mask _in_ the telluric regions, and
                                ; do not widen the masked regions.

                                ; The oxygen bands, then the rest.
    nord=3
    o2mask = spflux_masklines(geoll[isort], hwidth=0., /o2, spans=o2spans)
    ii = where(ivar[isort] GT 0 AND o2mask EQ 1, o2sel)
    if (o2sel GT 0) then begin
       o2bkpts = bspline_bkpts(geoll[isort[ii]], everyn=nspec, nord=nord)
       o2set = bspline_iterfit(geoll[isort], flatflux * o2mask + (1-o2mask), $
                               invvar=ivar[isort] * o2mask, lower=3, upper=3, fullbkpt=o2bkpts, $
                               outmask=outmaskt, nord=nord, requiren=1)
    end else begin
       o2bkpts = 0
       o2set = 0
    end
                                ; Avoid the bitter end of the spectrum
    wmax = alog10(10^(max(geoll))-20)
    h2omask = spflux_masklines(geoll[isort], hwidth=0., /h2o, extrabands=[[wmax],[max(geoll)]], spans=h2ospans)
    ii = where(ivar[isort] GT 0 AND h2omask EQ 1)
    
                                ; Force the templates to zero slope at 1?
    h2obkpts = bspline_bkpts(geoll[isort[ii]], everyn=nspec, nord=nord)
    h2oset = bspline_iterfit(geoll[isort], flatflux * h2omask + (1-h2omask), $
                            invvar=ivar[isort] * h2omask, lower=3, upper=3, fullbkpt=h2obkpts, $
                            outmask=outmaskt, nord=nord, requiren=1)

    if o2sel GT 0 then begin
       sxaddpar, nhdr0, 'NTEMPLAT', 2
       mwrfits, 0, filename, nhdr0, /create
       sxaddpar, nhdr1, 'TMPLTYPE', 2
       mwrfits, o2set, filename, nhdr1
       mwrfits, o2spans, filename
    end else begin
       sxaddpar, nhdr0, 'NTEMPLAT', 1
       mwrfits, 0, filename, nhdr0, /create
    end

    
    for i=0,n_elements(fcfiles) do begin
       keyname = string('FILE',i,format='(a4,i03)')
       sxaddpar, nhdr0, keyname, files[i], 'input file for templates'
    end

    sxaddpar, nhdr2, 'TMPLTYPE', 1
    mwrfits, h2oset, filename, nhdr2
    mwrfits, h2ospans, filename

    if keyword_set(doplot) then begin
       !P.MULTI = [0,1,2]
       if color eq 'r' then begin
          for i=0,1 do begin
             wmin = o2spans[i,0]
             wmax = o2spans[i,1]

             sval = bspline_valu(geoll[isort], o2set)
             ssel = geoll[isort] GE wmin and geoll[isort] LE wmax
             resid = flatflux - (sval*ssel + 1-ssel)

             djs_plot,10^geoll[isort],flatflux,psym=3,yrange=[0,1.5],$
                      xrange=[10^wmin-10, 10^wmax+30]
             djs_oplot,10^geoll[isort],geoll*0+1,color='yellow'
             djs_oplot,10^geoll[isort],sval*o2mask,psym=3,color='red'

             djs_plot,10^geoll[isort],resid,psym=3,yrange=[-0.2,0.2],$
                      xrange=[10^wmin-10, 10^wmax+30]
             djs_oplot,10^geoll[isort],geoll*0,color='yellow'
             
             if !d.name eq 'X' and keyword_set(debug) then begin
                k = get_kbrd()
             end
             if keyword_set(giffile) then begin
                tvlct,r,g,b,/get ; Only works with 8-bit displays. I can't handle this issue right now. 
                write_gif, giffile, tvrd(), r, g, b, /multiple
                print, "put MPEG frame...."
             end
          end
       end
          
       for i=0,(size(h2ospans,/dim))[0]-1 do begin
          wmin = h2ospans[i,0]
          wmax = h2ospans[i,1]
          
          if color eq 'b' then begin
             if 10^wmin GT 6000 then break
          end else begin
             if 10^wmax LT 5800 then continue
          end
          
          sval = bspline_valu(geoll[isort], h2oset)
          ssel = geoll[isort] GE wmin and geoll[isort] LE wmax
          resid = flatflux - (sval*ssel + 1-ssel)
          
          djs_plot,10^geoll[isort],flatflux,psym=3,yrange=[0,1.5],$
                   xrange=[10^wmin-30, 10^wmax+30]
          djs_oplot,10^geoll[isort],geoll*0+1,color='yellow'
          djs_oplot,10^geoll[isort],sval*h2omask,psym=3,color='red'
          
          djs_plot,10^geoll[isort],resid,psym=3,yrange=[-0.2,0.2],$
                   xrange=[10^wmin-10, 10^wmax+30]
          djs_oplot,10^geoll[isort],geoll*0,color='yellow'
          
          if !d.name eq 'X' and keyword_set(debug) then begin
             k = get_kbrd()
          end
          if keyword_set(giffile) then begin
             tvlct,r,g,b,/get
             write_gif, giffile, tvrd(), r, g, b, /multiple
             print, "put MPEG frame...."
          end
       end
       !P.MULTI = [0,1,1]

    end

    if keyword_set(giffile) then begin
       write_gif,giffile,/close
    end

 end

