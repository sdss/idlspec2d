pro rline_findew, spec, sigma=sigma, ew, ewinv, fopt, finv

    if NOT keyword_set(sigma) then sigma=1.0

    boxarea = 5*sigma 
    hpix = fix(boxarea-1)/2
    intarea = 2*hpix+1
    boxcar = replicate(1,intarea)
    if boxarea-intarea GT 0 then $
      boxcar = [0.5*(boxarea-intarea),boxcar, 0.5*(boxarea-intarea)]
    
 ;----------------------------------------------------------------------------
 ;
 ;  Do boxcar convolution first
 ;
 ;----------------------------------------------------------------------------

    mask = float(spec.finv LE 0)
    nspec = n_elements(spec)

    smask = convol(mask, boxcar, /edge_truncate)

    var = 1.0/(spec.finv + mask)
  
    sub = spec.flux - spec.model 
    zero = where(mask)
    if zero[0] NE -1 then sub[zero] = 0

    ew = convol(sub, boxcar, /edge_truncate)
    ewvar = convol(var, boxcar, /edge_truncate) 

    ewinv = ewvar * 0.0
    good = where(mask EQ 0 AND ewvar GT 0)
    if good[0] NE -1 then ewinv[good] = 1.0/ewvar[good]

 ;----------------------------------------------------------------------------
 ;
 ;  Now do optimized gaussian detection
 ;
 ;----------------------------------------------------------------------------

    gkernel = gauss_kernel(sigma) 

    finv = convol(spec.finv, gkernel^2, /edge_truncate)
    fnumerator = convol(sub*spec.finv, gkernel, /edge_truncate)

    smask = convol(mask, gkernel, /edge_truncate) 
    good = where(smask LT 0.4 AND finv GT 0)
    fopt = 0.0*finv

    if good[0] NE -1 then fopt[good] = fnumerator[good] / finv[good]
 
    return 
end 
