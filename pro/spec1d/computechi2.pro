;------------------------------------------------------------------------------
function computechi2, objflux, sqivar, starflux, $
 acoeff=acoeff, dof=dof, yfit=yfit

   ndim = size(starflux, /n_dimen)
   if (ndim EQ 1) then nstar = 1 $
    else nstar = (size(starflux, /dimens))[1]

   bvec = objflux * sqivar
   mmatrix = starflux
   for i=0L, nstar-1 do $
    mmatrix[*,i] = mmatrix[*,i] * sqivar
   mmatrixt = transpose( mmatrix )
   mm = mmatrixt # mmatrix

   ; Use SVD to invert the matrix
;   mmi = invert(mm, /double)
   if (nstar EQ 1) then begin
      mmi = 1.0 / mm
   endif else begin
      svdc, mm, ww, uu, vv, /double
      mmi = 0 * vv
      for i=0L, nstar-1 do mmi[i,*] = vv[i,*] / ww[i]
      mmi = mmi ## transpose(uu)
   endelse

   acoeff = mmi # (mmatrixt # bvec)
   chi2 = total( (mmatrix # acoeff - bvec)^2 )

   if (arg_present(yfit)) then $
    yfit = acoeff ## starflux
;    yfit = transpose(acoeff # mmatrixt) / sqivar
   if (arg_present(dof)) then $
    dof = total(sqivar NE 0) - nstar

   return, chi2
end
;------------------------------------------------------------------------------
