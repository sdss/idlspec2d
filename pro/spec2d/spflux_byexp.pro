function spflux_byexp,camera,loglam,flux,ivar, $
                      helio_rv, airmass=airmass, $
                      coeffs=coeffs,resids=resids, bkpts=bkpts, $
                      xrange=xrange,yrange=yrange, yfit=yfit, $
                      keywait=keywait,giffile=giffile

    !P.MULTI = [0,1,2]

    color = strmid(camera,0,1)

                                ; 1. Fit a provisional base spline
                                ; 2. Fit template splines using that spline.
                                ; 3. Pin the endpoints of the templates
                                ; 4. Refit a final base spline using
                                ;    the pinned endpoints.

    expsset = spflux_basespline(color, loglam, flux, ivar, bkpts=bkpts, airmass=airmass)
    spflux_fittemplates, camera, loglam, flux, ivar, expsset, helio_rv, $
                         coeffs=coeffs, yfit=yfit, bkpts=bkpts, airmass=airmass, $
                         /doplot, xrange=xrange,yrange=yrange, giffile=giffile

                                ; We need to return a bspline which
                                ; covers the full wavelength range. So
                                ; fit a new one. 
    wsort = sort(loglam)
    x2 = keyword_set(airmass) ? airmass[wsort] : 0
    e2 = bspline_iterfit(loglam[wsort], yfit, requiren=2, nbkpts=n_elements(yfit), $
                         x2=x2, npoly=2*keyword_set(airmass))
    v2 = bspline_valu(loglam[wsort], e2, x2=x2)

    if arg_present(resids) then begin
       rlli = sort(loglam)
       rll = loglam[rlli]
       rflux = flux[rlli]
       resids = (rflux-yfit) / yfit
    end

    if keyword_set(keywait) then begin
       key=get_kbrd()
       if key eq 'q' then return,0
    end

    !P.MULTI = [0,1,1]

    if keyword_set(giffile) then begin
       write_gif,giffile,/close
    end

    return, e2
 end

