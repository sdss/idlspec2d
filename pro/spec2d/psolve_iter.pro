; BUGS:
;  -- Differentiate between FLUX of a star, and the multiplier for the PSF.
;  -- What to do about stars where max shift exceeded?  Toss them from fits???
;  -- The fit PSF can have its centroid drift away from 0,0.
;     We must somehow enforce that it have a flux-weighted center at 0,0,
;     preferably for the PSF evaluated anywhere on the image - is that possible?
;  -- Initial iteration should use only the most compact objects (small Neff?)
;  -- Should add+delete objects after each iteration.
;  -- Call to FIND should do a better job setting HMIN.
;  -- Distribute the stars intelligently across the image.
;  -- Dynamically choose the number of polynomial terms to use.
;  -- Be sure to reject stars with bad chi^2 (at least after first iter)
;  -- Reject stars that overlap non-stellar objects
;  -- Limit the total number of objects to prevent enormous CPU times
;------------------------------------------------------------------------------
; Shift the object centers such that if we rederive the PSF eigen-templates,
; they should all have flux-weighted centroids at 0,0.
pro psolve_shift_centers, psfimg, objs

   kdims = size(psfimg, /dimens)
   nhalfx = (kdims[0]-1)/2
   nhalfy = (kdims[1]-1)/2
   if (size(psfimg,/n_dimen) EQ 2) then nmap = 1 $
    else nmap = kdims[2]

   ; Compute the flux-weighted center of each PSF eigen-template
   xnumer = fltarr(nmap)
   ynumer = fltarr(nmap)
   xdenom = fltarr(nmap)
   ydenom = fltarr(nmap)
   xaxis = djs_laxisgen(kdims[0:1], iaxis=0) - nhalfx
   yaxis = djs_laxisgen(kdims[0:1], iaxis=1) - nhalfy
   for imap=0, nmap-1 do begin
      xnumer[imap] = total(xaxis * psfimg[*,*,imap])
      ynumer[imap] = total(yaxis * psfimg[*,*,imap])
      xdenom[imap] = total(psfimg[*,*,imap])
      ydenom[imap] = total(psfimg[*,*,imap])
   endfor

   nobj = n_elements(objs)
   for iobj=0L, nobj-1L do begin
      objs[iobj].xcen = objs[iobj].xcen $
       + total(xnumer * objs[iobj].pvalue[0:nmap-1]) $
       / total(xdenom * objs[iobj].pvalue[0:nmap-1])
      objs[iobj].ycen = objs[iobj].ycen $
       + total(ynumer * objs[iobj].pvalue[0:nmap-1]) $
       / total(ydenom * objs[iobj].pvalue[0:nmap-1])
   endfor

   return
end
;------------------------------------------------------------------------------
function psolve_solve_matrix, mmatrix, bvec, chi2=chi2

   mmatrixt = transpose( mmatrix )
   atcinva = mmatrixt # mmatrix
   atcinvb = mmatrixt # bvec

;   if (n_elements(atcinva) EQ 1) then begin
;      mm = mmatrixt # mmatrix
;      mmi = 1.0 / (mm + (mm EQ 0))
;      acoeff = mmi # (mmatrixt # bvec)
;   endif else begin
;      svdc, atcinva, ww, uu, vv, /double
;      acoeff = svsol(uu, ww, vv, atcinvb, /double)
;   endelse

; Should use SVD instead of this...???
   mmi = invert(atcinva)
   acoeff = mmi # atcinvb

   chi2 = total( (mmatrix # acoeff - bvec)^2, /double )

   return, acoeff
end
;------------------------------------------------------------------------------
function psolve_kern, nkx, nky, dx, dy
   eps     = 1.0e-5 ; Smallest fractional shift allowed.
   dampfac = 3.25   ; Damping factor for gaussian.

   xmid = (nkx - 1) / 2
   ymid = (nky - 1) / 2
   
   ; Do not sinc shift if very nearly a full-pixel shift
   if (abs(dx MOD 1) GT eps) then begin
      xvec = findgen(nkx) - xmid + dx
      xkern = exp(-(xvec/dampfac)^2) * sin(!pi*xvec) / (!pi*xvec)
   endif else begin
      xkern = fltarr(nkx)
      xkern[xmid] = 1.
   endelse
   if (abs(dy MOD 1) GT eps) then begin
      yvec = findgen(nky) - ymid + dy
      ykern = exp(-(yvec/dampfac)^2) * sin(!pi*yvec) / (!pi*yvec)
   endif else begin
      ykern = fltarr(nky)
      ykern[ymid] = 1.
   endelse

   kern = xkern # ykern

   return, reform(kern, nkx*nky)
end
;------------------------------------------------------------------------------
function psolve_kern_dxy, nkx, nky, dx, dy, psfpix, psfimg

   if (size(psfimg,/n_dimen) NE 2) then $
    message, 'PSFIMG must be 2-dimensional in this function!'

   eps     = 1.0e-5 ; Smallest fractional shift allowed.
   dampfac = 3.25   ; Damping factor for gaussian.

   iy = psfpix / nkx
   ix = psfpix - iy * nkx

   ; Do not sinc shift if very nearly a full-pixel shift
   if (abs(dx MOD 1) GT eps) then begin
      xvec = findgen(nkx) - ix + dx
      dfdx1 = exp(-(xvec/dampfac)^2) * (cos(!pi*xvec)/xvec $
       - (2.*xvec/dampfac^2 + 1./xvec)*sin(!pi*xvec)/(!pi*xvec))
   endif else begin
      dfdx1 = fltarr(nkx)
   endelse
   if (abs(dy MOD 1) GT eps) then begin
      yvec = findgen(nky) - iy + dy
      dfdy1 = exp(-(yvec/dampfac)^2) * (cos(!pi*yvec)/yvec $
       - (2.*yvec/dampfac^2 + 1./yvec)*sin(!pi*yvec)/(!pi*yvec))
   endif else begin
      dfdy1 = fltarr(nky)
   endelse

   return, [total(dfdx1 * psfimg[*,iy]), total(dfdy1 * psfimg[ix,*])]
end
;------------------------------------------------------------------------------
function fn_psolve_fxy, params, thisimg=thisimg, thisivar=thisivar, $
 thispsf=thispsf

   chivec = sqrt((thisimg - params[0] * sshift2d(thispsf,params[1:2]))^2 $
    * thisivar)

   return, reform(chivec,n_elements(chivec))
end
;------------------------------------------------------------------------------
; Resolve for the flux + position of one star.
function psolve_mpsolvefxy1, image, ivar, psfimg, obj1, maxshift=maxshift

   dims = size(psfimg, /dimens)
   nkx = dims[0]
   nky = dims[1]
   if (size(psfimg,/n_dimen) EQ 2) then nmap = 1 $
    else nmap = dims[2]

   nn = nkx * nky ; nn=number of pixels to use in the image
   if (n_elements(obj1) NE 1) then $
    message, 'OBJ1 must contain exactly 1 object'

   xcen0 = round(obj1.xcen)
   ycen0 = round(obj1.ycen)
   xcen1 = obj1.xcen - xcen0
   ycen1 = obj1.ycen - ycen0

   thisimg = image[xcen0-nkx/2:xcen0+nkx/2,ycen0-nky/2:ycen0+nky/2]
   thisivar = ivar[xcen0-nkx/2:xcen0+nkx/2,ycen0-nky/2:ycen0+nky/2]

   thispsf = 0
   for imap=0, nmap-1 do $
    thispsf = thispsf + obj1.pvalue[imap] * psfimg[*,*,imap]

   newobj = obj1

   if (maxshift GT 0) then begin
      ;----------
      ; Fix flux + centers (non-linear fit using MPFIT)

      parinfo = replicate({value:0.D, limited:[0,0], $
       limits:[0.D,0], step: 0.D, mpmaxstep: 0.D}, 3)
      parinfo.value = [obj1.flux, xcen1, ycen1]
      parinfo[1:2].limited = [1,1]
      ; Arbitrarily set the fitting range to [-1,1] from the fixed value
      ; if the center gets fixed
      parinfo[1].limits = xcen1 + [-1,1]*(maxshift + (maxshift EQ 0))
      parinfo[2].limits = ycen1 + [-1,1]*(maxshift + (maxshift EQ 0))
      parinfo.step = [0, 0.01, 0.01]
      parinfo.mpmaxstep = [0, 0.10, 0.10]

; What about masked pixels???  How many iterations???
      functargs = {thisimg: thisimg, thisivar: thisivar, thispsf: thispsf}
      acoeff = mpfit('fn_psolve_fxy', parinfo.value, $
       parinfo=parinfo, functargs=functargs, $
       maxiter=100, niter=niter, status=status, perror=perror, /quiet)
      if (status LT 0) then perror = [-1,-1,-1] ; Case where fit failed
      chi2 = total(fn_psolve_fxy(acoeff, thisimg=thisimg, thisivar=thisivar, $
       thispsf=thispsf)^2)

      newobj.fitstatus = status
      newobj.fititer = niter
      newobj.flux = acoeff[0]
      newobj.xcen = xcen0 + acoeff[1]
      newobj.ycen = ycen0 + acoeff[2]
      newobj.flux_err = perror[0]
      newobj.xcen_err = perror[1]
      newobj.ycen_err = perror[2]
      newobj.chi2 = chi2
      newobj.dof = nn - 3
      if (acoeff[0] LE 0) then newobj.fitstatus = 2
   endif else begin
      ;----------
      ; Fit flux only (linear fit)

      newobj.chi2 = computechi2(thisimg[*], sqrt(thisivar[*]), $
       (sshift2d(thispsf,[xcen1,ycen1]))[*], acoeff=acoeff, dof=dof, var=var)
      newobj.fitstatus = 1
      newobj.fititer = 0
      newobj.flux = acoeff
      newobj.flux_err = sqrt(var)
      newobj.xcen_err = 0
      newobj.ycen_err = 0
      newobj.dof = dof
   endelse

   indx = where(thisivar GT 0, ct)
   if (ct GT 0) then begin
      newobj.neff = djs_neff(thisimg[indx], 1./sqrt(thisivar[indx]), nerr=nerr)
      newobj.nerr = nerr
   endif else begin
      newobj.neff = 0
      newobj.nerr = 0
   endelse

   return, newobj
end
;------------------------------------------------------------------------------
; Resolve for the flux + position of one star.
function psolve_solvefxy1, image, ivar, psfimg, obj1, maxshift=maxshift

   dims = size(psfimg, /dimens)
   nkx = dims[0]
   nky = dims[1]
   if (size(psfimg,/n_dimen) EQ 2) then nmap = 1 $
    else nmap = dims[2]

   nn = nkx * nky ; nn=number of pixels to use in the image
   if (n_elements(obj1) NE 1) then $
    message, 'OBJ1 must contain exactly 1 object'

   xcen0 = round(obj1.xcen)
   ycen0 = round(obj1.ycen)
   xcen1 = obj1.xcen - xcen0
   ycen1 = obj1.ycen - ycen0

   mmatrix = dblarr(nn, 3)
   bvec = dblarr(nn)
   for imap=0L, nmap-1 do begin
      mmatrix[*,0] = mmatrix[*,0] + obj1.pvalue[imap] $
       * convol(psfimg[*,*,imap], $
       (reform(psolve_kern(nkx,nky,xcen1,ycen1),nkx,nky))[1:nkx-2,1:nky-2], $
       /center, /edge_truncate) ; This is equivalent to a sinc shift
      thispsf = obj1.pvalue[imap] * psfimg[*,*,imap]
      for i=0L, nn-1L do begin
         mmatrix[i,[1,2]] = mmatrix[i,[1,2]] + $
          psolve_kern_dxy(nkx, nky, xcen1, ycen1, i, thispsf)
      endfor
   endfor

   thisisig = reform( $
    sqrt( ivar[xcen0-nkx/2:xcen0+nkx/2,ycen0-nky/2:ycen0+nky/2] ), nkx*nky)
   thisimg = reform( image[xcen0-nkx/2:xcen0+nkx/2,ycen0-nky/2:ycen0+nky/2], $
    nkx*nky)
   for j=0, 2 do $
    mmatrix[*,j] = mmatrix[*,j] * thisisig
   bvec = thisimg * thisisig

   acoeff = psolve_solve_matrix(mmatrix, bvec, chi2=chi2)

   newobj = obj1
   newobj.fitstatus = 0

   ; Do not allow any shifts larger than MAXSHIFT
   xdiff = acoeff[1] - maxshift
   if (abs(xdiff) GT maxshift) then begin
      acoeff[1] = (acoeff[1] > (-maxshift)) < maxshift
      newobj.fitstatus = 1
   endif
   ydiff = acoeff[2] - maxshift
   if (abs(ydiff) GT maxshift) then begin
      acoeff[2] = (acoeff[2] > (-maxshift)) < maxshift
      newobj.fitstatus = 1
   endif

   newobj.flux = acoeff[0]
   newobj.xcen = newobj.xcen + acoeff[1]
   newobj.ycen = newobj.ycen + acoeff[2]
   newobj.chi2 = chi2
   newobj.dof = nn - 3
   if (acoeff[0] LT 0 OR newobj.fitstatus NE 0) then begin
      ; If MAXSHIFT is exceeded or flux is negative, then fit only for the flux
      acoeff = psolve_solve_matrix(mmatrix[*,0], bvec, chi2=chi2)
      newobj.flux = acoeff[0]
      newobj.chi2 = chi2
      newobj.dof = nn - 1
      newobj.fitstatus = 1
   endif
   if (acoeff[0] LE 0) then newobj.fitstatus = 2

   indx = where(thisisig GT 0)
   newobj.neff = djs_neff(thisimg[indx], 1./thisisig[indx], nerr=nerr)
   newobj.nerr = nerr

   return, newobj
end
;------------------------------------------------------------------------------
function psolve_solvepsf, image, ivar, pixmask, psfpix, objs

   splog, 'Number of PSF stars = ', n_elements(objs)
   nn = max(pixmask)+1 ; nn=number of pixels to use in the image
   nz = max(psfpix)+1 ; nz=number of PSF pixels to solve for
   pndim = size(psfpix, /n_dimen)
   pdim = size(psfpix, /dimens)
   nkx = pdim[0]
   nky = pdim[1]
   if (pndim EQ 2) then nmap = 1 $
    else nmap = pdim[2]
   dims = size(image, /dimens)

   ; Construct the M matrix

   mmatrix = dblarr(nn, nz)
   for iobj=0L, n_elements(objs)-1 do begin
      xint = round(objs[iobj].xcen)
      yint = round(objs[iobj].ycen)
      dx = objs[iobj].xcen - xint
      dy = objs[iobj].ycen - yint
      thiskern = psolve_kern(nkx, nky, -dx, -dy)
      for ix=0L,nkx-1 do begin
         for iy=0L,nky-1 do begin
            ; If the object is too close to the edge of the image,
            ; then we know there will be no overlapping object pixels.
;            thismask = pixmask[xint-nkx+ix+1:xint+ix,yint-nky+iy+1:yint+iy]
            x1 = xint - nkx + ix + 1
            x2 = xint + ix
            y1 = yint - nky + iy + 1
            y2 = yint + iy
            if (x1 GE 0 AND x2 LT dims[0] $
             AND y1 GE 0 AND y2 LT dims[1]) then $
             thismask = pixmask[x1:x2,y1:y2] $
            else $
             thismask = -1L
            for imap=0, nmap-1 do begin
               ipix = psfpix[ix,iy,imap]
               if (ipix NE -1) then begin
                  jj = where(thismask NE -1, ct)
                  if (ct GT 0) then $
                   mmatrix[thismask[jj],ipix] += thiskern[jj] $
                    * objs[iobj].pvalue[imap] * objs[iobj].flux
               endif
            endfor
         endfor
      endfor
   endfor

   ii = where(pixmask NE -1)
   isig = sqrt(ivar[ii])
   bvec = double( image[ii] * isig )
   for i=0L, nz-1 do mmatrix[*,i] *= isig

   psfimg = fltarr(pdim)
   res = psolve_solve_matrix(mmatrix, bvec, chi2=chi2)
   for i=0L, nz-1 do psfimg[where(psfpix EQ i)] = res[i]

   dof = nn - nz - 3*n_elements(objs)
   splog, 'Chi2/dof = ', chi2, '/', dof, ' = ', chi2/dof

   return, psfimg
end
;------------------------------------------------------------------------------
; Select the best (brightest) objects to use as PSF stars
pro psolve_selectbest, objs, nselect=nselect

   objs.bestmask = 0B

   igood = where(objs.goodmask EQ 1, ngood)
   if (ngood LE nselect) then begin
      objs[igood].bestmask = 1B
   endif else begin
      flux = objs[igood].flux
      fluxlim = flux[(reverse(sort(flux)))[nselect-1]]
      ibest = where(flux GE fluxlim)
      objs[igood[ibest]].bestmask = 1B
   endelse

   splog, 'Selected ', fix(total(objs.bestmask)), ' of ', $
    fix(total(objs.goodmask)), ' good stars'

   return
end
;------------------------------------------------------------------------------
; Return the footprint of the PSF, where =1 for pixels used, =0 otherwise

function psolve_footprint, psfpix
   if (size(psfpix,/n_dimen) EQ 2) then kern = psfpix GE 0 $
    else kern = total(psfpix GE 0, 3) GT 0
   return, kern
end
;------------------------------------------------------------------------------
; Measure the 1st moment of all objects
pro psolve_moments, image, ivar, objs, psfpix, skyimg

   pdim = size(psfpix, /dimens)
   nkx = pdim[0]
   nky = pdim[1]
   nhalfx = (pdim[0]-1)/2
   nhalfy = (pdim[1]-1)/2

   xcen0 = round(objs.xcen)
   ycen0 = round(objs.ycen)
   xcen1 = objs.xcen - xcen0
   ycen1 = objs.ycen - ycen0

   diffimg = image - skyimg

   kern = psolve_footprint(psfpix)
   rimg = shift(dist(nkx, nky), nhalfx, nhalfy) * kern

   for iobj=0L, n_elements(objs)-1 do begin
      xint = xcen0[iobj]
      yint = ycen0[iobj]
      thisimg = diffimg[xint-nhalfx:xint+nhalfx,yint-nhalfy:yint+nhalfy]
      thisimg = sshift2d(thisimg, [-xcen1[iobj], -ycen1[iobj]])
      thismask = ivar[xint-nhalfx:xint+nhalfx,yint-nhalfy:yint+nhalfy] NE 0
      objs[iobj].rmoment = total(rimg * thisimg * thismask) $
       / total(thisimg * thismask)
   endfor

   return
end
;------------------------------------------------------------------------------
; Mark any objects as saturated if they have a bad pixel within RADIUS pix
; of the center.
pro psolve_saturated, objs, ivar, radius=radius

   width = 2*fix(radius) + 1
   mask = convol((ivar EQ 0), intarr(width,width)+1, $
    /center, /edge_truncate)

   objs.qsatur = mask[objs.xcen,objs.ycen] NE 0
   return
end
;------------------------------------------------------------------------------
; Construct an image of which pixels to use in the image.
; Unused pixels are given a value of -1, and others are enumerated
; with each used pixel having a unique integer from [0,N-1].
function psolve_mask, objs, ivar, psfpix, itouch=itouch

   ;----------
   ; Construct the mask of which pixels to use in the PSF fit

   dims = size(ivar, /dimens)
   pixmask = fltarr(dims)

   kern = psolve_footprint(psfpix)

   ibest = where(objs.bestmask EQ 1)
   populate_image, pixmask, objs[ibest].xcen, objs[ibest].ycen
   pixmask = convol(pixmask, kern, /center, /edge_truncate)

   pixmask = pixmask AND (ivar GT 0)
   ii = where(pixmask, nn)
   retmask = lonarr(dims) - 1L
   retmask[ii] = lindgen(nn)

   ;----------
   ; Decide which objects touch the pixel mask (may not be all
   ; the so-called PSF stars, and may include some other objects)

   kern = psolve_footprint(psfpix)
   qtouch = bytarr(n_elements(objs))
   pdim = size(psfpix, /dimens)
   nhalfx = (pdim[0]-1)/2
   nhalfy = (pdim[1]-1)/2
   for i=0L, n_elements(objs)-1 do begin
      xint = round(objs[i].xcen)
      yint = round(objs[i].ycen)
      qtouch[i] = total(pixmask[xint-nhalfx:xint+nhalfx, $
       yint-nhalfx:yint+nhalfx] * kern) GT 0
   endfor
   itouch = where(qtouch, ntouch)

   return, retmask
end
;------------------------------------------------------------------------------
function psolve_specsky, image, ivar, xpad=xpad, ypad=ypad

   dims = size(image, /dimens)
   skyimg = fltarr(dims)

   for i=0, dims[0]-1 do $
    skyimg[i,*] = djs_median(reform(image[i,*],dims[1]), $
     width=99, boundary='reflect')

   return, skyimg
end
;------------------------------------------------------------------------------
function psolve_sky, image, ivar, xpad=xpad, ypad=ypad

   dims = size(image, /dimens)
   skyimg = fltarr(dims)

   diffimg = image
   for i=0, 2 do begin
      skyvecx = $
       djs_avsigclip(diffimg[xpad:dims[0]-xpad-1,ypad:dims[1]-ypad-1], $
       2, inmask=(ivar EQ 0))
      skyimg[xpad:dims[0]-xpad-1,ypad:dims[1]-ypad-1] += $
       skyvecx # (fltarr(1,dims[1]-2*ypad)+1)
      diffimg = image - skyimg

      skyvecy = $
       djs_avsigclip(diffimg[xpad:dims[0]-xpad-1,ypad:dims[1]-ypad-1], $
       1, inmask=(ivar EQ 0))
      skyimg[xpad:dims[0]-xpad-1,ypad:dims[1]-ypad-1] += $
       (fltarr(dims[0]-2*xpad)+1) # skyvecy
      diffimg = image - skyimg
   endfor

   return, skyimg
end
;------------------------------------------------------------------------------
pro psolve_find, image, ivar, objs, hmin=hmin

   ;----------
   ; Run DAO find

   if (NOT keyword_set(objs)) then begin
;      sharplim = [0.0, 1.0]
;      roundlim = [-1.0, 1.0]
      sharplim = [-1e10, 1e10]
      roundlim = [-1e10, 1e10]
      fwhm = 3. ; ???
      if (NOT keyword_set(hmin)) then return

      find, image, xcen, ycen, flux, sharp, round, hmin, fwhm, roundlim, $
       sharplim, /silent
      splog, 'Number of objects found = ', n_elements(xcen)
      if (NOT keyword_set(xcen)) then return

      objs = psolve_obj_struct(n_elements(xcen))
      objs.xcen = xcen
      objs.ycen = ycen
      objs.flux = flux
   endif

   psolve_saturated, objs, ivar, radius=3

   return
end
;------------------------------------------------------------------------------
; XPAD,YPAD are what the image is already padded with...???
; MAXSHIFT is a required keyword
pro psolve_iter, image, ivar, objs, psfpix, psfimg, skyimg, $
 xpad=xpad, ypad=ypad, $
 nselect=nselect, npoly=npoly, niter=niter, maxshift=maxshift, fixpsf=fixpsf, $
 filename=filename

   if (keyword_set(filename)) then restore, filename

   dims = size(image, /dimens)
   t0 = systime(1)

   ;----------
   ; Solve for the sky

   qfitsky = keyword_set(skyimg) EQ 0

   if (qfitsky) then $
    skyimg = psolve_specsky(image, ivar, xpad=xpad, ypad=ypad)

   ;----------
   ; Measure the 1st moment of all objects

   psolve_moments, image, ivar, objs, psfpix, skyimg

   ;---------------------------------------------------------------------------
   ; MAIN ITERATION LOOP
   ;---------------------------------------------------------------------------

   for iiter=0L, niter-1L do begin
      splog, 'ITERATION #', iiter+1, ' of ', niter

      ;----------
      ; Select the objects to use in the PSF fit

      if (keyword_set(nselect)) then objs.goodmask = 1B
      objs.goodmask = objs.goodmask AND objs.qsatur EQ 0 $
       AND objs.xcen GT xpad+maxshift AND objs.xcen LT dims[0]-xpad-1-maxshift $
       AND objs.ycen GT ypad+maxshift AND objs.ycen LT dims[1]-ypad-1-maxshift
      igood = where(objs.goodmask, ngood)
      if (ngood EQ 0) then $
       message, 'No good stars?!'

      if (total(objs.chi2) EQ 0) then begin
         ; If no chi^2 values computed yet for individual objects,
         ; then select only those objects whose radial moments are
         ; not more than 1.20 times larger than the 25-th percentile value.
         isort = igood[sort(objs[igood].rmoment)]
         rmom_cut = 1.20 * objs[isort[0.25*ngood]].rmoment
         splog, 'Selecting stars with RMOMENT < ', rmom_cut
         objs.goodmask = objs.goodmask $
          AND objs.rmoment LT rmom_cut
      endif else begin
         rchi2 = objs.chi2 / (objs.dof > 1)
         isort = igood[sort(rchi2[igood])]
         rchi2_cut = 1.20 * rchi2[isort[0.25*ngood]]
         splog, 'Selecting stars with RCHI2 < ', rchi2_cut
         objs.goodmask = objs.goodmask $
          AND rchi2 LT rchi2_cut
      endelse

      if (keyword_set(nselect)) then $
       psolve_selectbest, objs, nselect=nselect $
      else $
       splog, 'Keeping previously-selected best stars'

      ;----------
      ; Compute the polynomial terms for each object based upon its location

; ???
;      xnorm = (objs.xcen - xpad - 0.5 * dims[0]) / float(dims[0])
;      ynorm = (objs.ycen - ypad - 0.5 * dims[1]) / float(dims[1])
      rc_scale = 1.e-3 ; Hard-wired scale factor for ypos/xpos coefficients
      xnorm = (objs.xcen - xpad) * rc_scale
      ynorm = (objs.ycen - ypad) * rc_scale
      for ix=0, npoly[0]-1 do $
       for iy=0, npoly[1]-1 do $
        objs.pvalue[ix+iy*npoly[0]] = xnorm^ix * ynorm^iy

      ;----------
      ; Construct the mask of which pixels to use in the PSF fit, and solve.

      splog, 'Solving PSF image'
      t1 = systime(1)
      if (keyword_set(fixpsf) EQ 0 OR keyword_set(psfimg) EQ 0) then begin
         pixmask = psolve_mask(objs, ivar, psfpix, itouch=itouch)
         psfimg = psolve_solvepsf(image - skyimg, ivar, pixmask, psfpix, $
          objs[itouch])
; ???
;         splog, 'Shifting PSF eigen-image centers and re-solving'
;         psolve_shift_centers, psfimg, objs
;         pixmask = psolve_mask(objs, ivar, psfpix, itouch=itouch)
;         psfimg = psolve_solvepsf(image - skyimg, ivar, pixmask, psfpix, $
;          objs[itouch])
      endif

      ;----------
      ; Construct the fake image

      t2 = systime(1)
      fakeimg = skyimg
      psolve_addstars, fakeimg, psfimg, objs

      ;----------
      ; Loop through all objects, refitting their fluxes+centroids

      splog, 'Refitting stars with maxshift=',maxshift
      irefit = lindgen(n_elements(objs)) ; Refit all objects
;      irefit = itouch ; Refit only PSF stars

      t2 = systime(1)
      diffimg = image - fakeimg
      for i=0L, n_elements(irefit)-1L do begin
         print, i+1, n_elements(irefit), string(13b), $
          format='("Star ",i," of ",i,a1,$)'
         j = irefit[i]
         psolve_addstars, diffimg, psfimg, objs[j]

;         objs[j] = psolve_solvefxy1(diffimg, ivar, psfimg, objs[j], $
;          maxshift=maxshift)
         objs[j] = psolve_mpsolvefxy1(diffimg, ivar, psfimg, objs[j], $
          maxshift=maxshift)

         obj1 = objs[j]
         obj1.flux = -obj1.flux
         psolve_addstars, diffimg, psfimg, obj1
      endfor
      print

      ;----------
      ; Re-solve for the sky image

      t3 = systime(1)
      if (qfitsky) then $
       skyimg = psolve_specsky(diffimg+skyimg, ivar, xpad=xpad, ypad=ypad)
      t4 = systime(1)

      splog, 'Iteration #', iiter, ': Time to fit PSF = ', t2-t1, ' sec'
      splog, 'Iteration #', iiter, ': Time to refit stars = ', t3-t2, ' sec'
      splog, 'Iteration #', iiter, ': Time to fit sky = ', t4-t3, ' sec'
   endfor

   splog, 'Total time = ', t4-t0, ' sec'

   if (keyword_set(filename)) then $
    save, filename=filename, psfimg, skyimg, objs

   return
end
;------------------------------------------------------------------------------
