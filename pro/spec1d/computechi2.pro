;------------------------------------------------------------------------------
function computechi2, objflux, sqivar, starflux, $
 acoeff=acoeff, dof=dof, yfit=yfit

   ndim = size(starflux, /n_dimen)
   if (ndim EQ 1) then nstar = 1 $
    else nstar = (size(starflux, /dimens))[1]

   bvec = objflux * sqivar
   mmatrix = starflux * (sqivar # replicate(1,nstar))

;   ---------  the line above is about twice as fast --------------
;   for i=0L, nstar-1 do $
;     mmatrix[*,i] = mmatrix[*,i] * sqivar

   mmatrixt = transpose( mmatrix )
   mm = mmatrixt # mmatrix

   ; Use SVD to invert the matrix
;   mmi = invert(mm, /double)
   if (nstar EQ 1) then begin
      ; The last term below is to protect against divide-by-zero
      ; in the degenerate case.
      mmi = 1.0 / (mm + (mm EQ 0))
   endif else begin
      svdc, mm, ww, uu, vv, /double
      mmi = 0 * vv
      ; The last term below is to protect against divide-by-zero
      ; in the degenerate case.
      for i=0L, nstar-1 do mmi[i,*] = vv[i,*] / (ww[i] + (ww[i] EQ 0))
      mmi = mmi ## transpose(uu)
   endelse

   acoeff = mmi # (mmatrixt # bvec)
   chi2 = total( (mmatrix # acoeff - bvec)^2, /double )

   if (arg_present(yfit)) then $
    yfit = acoeff ## starflux
;    yfit = transpose(acoeff # mmatrixt) / sqivar
   if (arg_present(dof)) then $
    dof = total(sqivar NE 0) - nstar

   return, chi2
end
;------------------------------------------------------------------------------
