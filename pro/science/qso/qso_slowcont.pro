
pro qso_slowcont, loglam, flux, invvar, range=range, model=model, mask=mask

    if (size(flux))[0] NE 1 then begin
         print, 'just 1d arrays for now'
         return
    endif

    good = where(invvar GT 0)
    if good[0] EQ -1 then return   
 
    spec_range=[min(loglam[good]),max(loglam[good])]
    if keyword_set(range) then begin
       range[0]= max([spec_range[0], range[0]])
       range[1]= min([spec_range[1], range[1]])
    endif else range=spec_range

    model = loglam * 0.0
    mask = long(model)

    inrange = where(loglam GE range[0] AND loglam LE range[1], nin)

    if inrange[0] EQ -1 then return

    maxiter = nin/5 > 10

    modelset = bspline_iterfit(loglam[inrange], flux[inrange], $
                  invvar=invvar[inrange], nord=4, yfit=yfit, $
                  maxiter=maxiter, upper=10.0, lower=3.0, /silent, $
                  outmask=outmask, everyn=20, maxrej=1, grow=5, /groupbadpix)

    diff = (flux[inrange] - yfit) * sqrt(invvar[inrange])
    close = where(diff LT 10.0 AND diff GT -2.0)
    if close[0] NE -1 then outmask[close] = 1

    modelset = bspline_iterfit(loglam[inrange], flux[inrange], $
                  invvar=invvar[inrange]*outmask, nord=4, yfit=yfit, $
                  maxiter=maxiter, upper=5.0, lower=2.0, /silent, $
                  outmask=outmask, everyn=15, maxrej=1, /groupbadpix, grow=3)

    diff = (flux[inrange] - yfit) * sqrt(invvar[inrange])
    close = where(diff LT 10.0 AND diff GT -2.0)
    if close[0] NE -1 then outmask[close] = 1

    modelset = bspline_iterfit(loglam[inrange], flux[inrange], $
                  invvar=invvar[inrange]*outmask, nord=4, yfit=yfit, $
                  maxiter=maxiter, upper=5.0, lower=2.0, /silent, $
                  outmask=outmask, everyn=15, maxrej=1, /groupbadpix, grow=1)

    if total(outmask EQ 0) LT maxiter then begin
      model[inrange] = yfit
      mask[inrange] = outmask NE 0 
    endif

    return
end

