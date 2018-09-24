function skysubtract_iter, objflux, objivar, plugsort, wset, objsub, objsubivar, $
 iskies=iskies, fibermask=fibermask, nord=nord, upper=upper, $
  lower=lower, maxiter=maxiter, pixelmask=pixelmask, thresh=thresh, $
   npoly=npoly, relchi2set=relchi2set, $
    novariance=novariance, tai=tai, nbkpt=nbkpt, newmask=newmask, sset=sset, niter=niter

    if (not keyword_set(niter)) then niter = 1

    objsub0 = objflux*1
    objsubivar0 = objivar*1

    for iter=0, niter-1 do begin
        skystruct = skysubtract(objsub0, objsubivar0, plugsort, wset, $
            objsub, objsubivar, iskies=iskies, nord=nord, $
            fibermask=fibermask, upper=upper, lower=lower, maxiter=maxiter, $
            pixelmask=pixelmask, thres=thres, npoly=npoly, relchi2set=relchi2set, $
            tai=tai, nbkpt=nbkpt,$
            newmask=newmask, sset=sset)
        if (NOT keyword_set(skystruct)) then return, 0
         
        objsub0=objsub
        objsubivar0=objsubivar
    endfor 

    return, skystruct
end
