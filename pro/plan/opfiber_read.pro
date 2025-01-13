function opfiber_read, opFiber_file, cartid, camname, mjd

    yanny_read,opFiber_file, opfiber_pars, hdr=hdr, /anonymous

    fiberparam=*opfiber_pars

    match=where((strtrim(cartid,2) eq strtrim(fiberparam.cartid,2)) and $
                (strtrim(CAMNAME,2) eq strtrim(fiberparam.CAMNAME,2) and $
                (long(mjd) ge long(fiberparam.mjd))),imatch)
    fiberparam = fiberparam[match]

    if (imatch gt 1) then begin
        i = where(fiberparam.mjd eq max(fiberparam.mjd),ni)
        fiberparam =fiberparam[i]
    endif
    if imatch eq 0 then begin
        splog, 'ERROR: No Matching opFiber FIBERPARAM in ', opFiber_file
        splog, /close
        exit
    endif
    

    
    if n_elements(fiberparam) gt 1 then $
        fiberparam = fiberparam[where(long(fiberparam.mjd) eq min(long(fiberparam.mjd)))]
    
    nfibers = fiberparam[0].nfiber
    nbundles = fiberparam[0].nbundles
    deadfibermask = intarr(nfibers)
    positions = fiberparam.deadfibermask
    valid = where(positions ne -1, n_valid)
    if n_valid gt 0 then begin
        positions = positions[valid]
        deadfibermask[positions - 1] = 1
    endif
    
    opfiber_par =  create_struct('cartid', fiberparam.cartid, $
                                 'camname', fiberparam.camname, $
                                 'mjd', fiberparam.mjd, $
                                 'xstart_tol', fiberparam.xstart_tol,$
                                 'nfiber', fiberparam.nfiber,$
                                 'nbundles',fiberparam.nbundles,$
                                 'fiberspace', fiberparam.fiberspace[0:nbundles-1],$
                                 'bundlegap', fiberparam.bundlegap[0:nbundles-1],$
                                 'bundlefibers', fiberparam.bundlefibers[0:nbundles-1],$
                                 'deadfibermask', deadfibermask )
;print_struct, opfiber_par
    splog, 'opfiber_par from:', strtrim(opfiber_par.mjd,2), $
          ' cartid:',strtrim(opfiber_par.cartid,2), $
          ' camname:', strtrim(opfiber_par.camname,2)
    return,opfiber_par
end
