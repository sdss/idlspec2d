;------------------------------------------------------------------------------
function computez, objflux, objivar, starflux, starivar, zoffset=zoffset, $
 zmin=zmin, zmax=zmax, chi2=chi2best

; ??? STARIVAR is not used anywhere!!!

   ;---------------------------------------------------------------------------
   ; Check dimensions of input arrays

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npixobj = dims[0]
   if (ndim EQ 1) then nobj = 1 $
    else if (ndim EQ 2) then nobj = dims[1] $
    else message, 'OBJFLUX is neither 1-D or 2-D'

   ndim = size(starflux, /n_dimen)
   dims = size(starflux, /dimens)
   npixstar = dims[0]
   if (ndim EQ 1) then nstar = 1 $
    else if (ndim EQ 2) then nstar = dims[1] $
    else message, 'STARFLUX is neither 1-D or 2-D'

;   if total(abs(size(starflux, /dimens)-size(starivar, /dimens))) NE 0 $
;    OR size(starflux, /n_dimen) NE size(starivar, /n_dimen) THEN  $
;    message, 'Dimensions of STARFLUX and STARIVAR do not match'

   if total(abs(size(objflux, /dimens)-size(objivar, /dimens))) NE 0 $
    OR size(objflux, /n_dimen) NE size(objivar, /n_dimen) THEN  $
    message, 'Dimensions of OBJFLUX and OBJIVAR do not match'

   pixoffset = long(zoffset)
   lags = - lindgen(zmax-zmin+1) + pixoffset
   nlag = n_elements(lags)

   chi2arr = fltarr(nobj, nlag)

   ;---------------------------------------------------------------------------

   sqivar = sqrt(objivar)
   bvec = objflux * sqivar

   for ilag=0, nlag-1 do begin
junk = min(chi2arr[iobj,*], ilag)
      i1 = (lags[ilag] < npixstar-1)
      if (i1 LT 0) then i0 = -i1 $
       else i0 = 0
      i1 = i1 > 0
      i2 = npixstar-1 < (npixobj+i1-i0-1)
      ncomp = i2 - i1 + 1
      mmatrixt = starflux[i1:i2,*] 
      for i=0, nstar-1 do $
       mmatrixt[*,i] = mmatrixt[*,i] * sqivar[i0:i0+ncomp-1]
      mmatrix = transpose( mmatrixt )
      mm = mmatrix # mmatrixt
      mmi = invert(mm, /double)

      for iobj=0, nobj-1 do begin
         acoeff = mmi # (mmatrix # bvec[i0:i0+ncomp-1,iobj])
         chi2arr[iobj,ilag] = $
          total( (mmatrixt # acoeff - bvec[i0:i0+ncomp-1,iobj])^2 )
 yfit = transpose(acoeff # mmatrix) / sqivar[i0:i0+ncomp-1] ; ???
 plot, smooth(objflux[i0:i0+ncomp-1,iobj],7)
 djs_oplot, yfit+5, color='red'
      endfor
print, 'LAG # ', ilag
   endfor

; Should actually fit to this peak with gaussian or parabola???
   zbest = fltarr(nobj)
   chi2best = fltarr(nobj)
   for iobj=0, nobj-1 do begin
      chi2best[iobj] = min(chi2arr[iobj,*], imin)
      zbest[iobj] = pixoffset - lags[imin]
   endfor
stop

   return, zbest
end

