;+
; NAME:
;   computez
;
; PURPOSE:
;   Compute relative redshift of object(s) vs. eigen-templates.
;
; CALLING SEQUENCE:
;   zstruct = computez(objflux, objivar, starflux, starivar, [nfind=, $
;    zoffset=, zmin=, zmax= ]
;
; INPUTS:
;   objflux    - Object fluxes [NPIXOBJ,NOBJ]
;   objivar    - Object inverse variances [NPIXOBJ,NOBJ]
;   starflux   - Eigen-template fluxes [NPIXSTAR,NTEMPLATE]
;   starivar   - Eigen-template inverse variances [NPIXSTAR,NTEMPLATE]
;
; OPTIONAL INPUTS:
;   nfind      - Number of solutions to find per object; default to 1.
;   zoffset    - Offset between all objects and templates, in pixels.
;                A value of 10 indicates that STARFLUX begins ten pixels
;                after OBJFLUX, i.e. OBJFLUX[i+10] = STARFLUX[i] for the
;                case when the relative redshift should be zero.  If the
;                wavelength coverage of the templates is larger, then the
;                value of ZOFFSET will always be negative.
;   zmin       - The smallest redshift to consider.
;   zmax       - The largest redshift to consider.
;
; OUTPUTS:
;   zstruct    - Output structure [NOBJECT,NFIND] with the following elements:
;                z : The relative redshift.
;                zerr : Error in the redshift, based upon the local quadratic
;                       fit to the chi^2 minimum. 
;                chi2 : Fit value for the best (minimum) chi^2
;                dof : Number of degrees of freedom, equal to the number of
;                      pixels in common between the object and templates
;                      minus the number of templates.
;                theta : Mixing angles [NTEMPLATE].  These are computed at the
;                        nearest integral redshift, e.g. at ROUND(ZOFFSET).
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/TEMPLATEFILES
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Jul-2000  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function computez, objflux, objivar, starflux, starivar, nfind=nfind, $
 zoffset=zoffset, zmin=zmin, zmax=zmax

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

   if (n_elements(zoffset) GT 1) then $
    message, 'ZOFFSET must be a scalar'
   pixoffset = round(zoffset)

   if (NOT keyword_set(nfind)) then nfind = 1
   if (n_elements(zoffset) EQ 0) then zoffset = 0
   if (n_elements(zmin) EQ 0) then $
    zmin = -2 * ((npixobj < npixstar) + 1) + pixoffset
   if (n_elements(zmax) EQ 0) then $
    zmax = zmin + 2 * ((npixobj < npixstar) - 1)

   if (n_elements(zmin) GT 1) then $
    message, 'ZMIN must be a scalar'
   if (n_elements(zmax) GT 1) then $
    message, 'ZMAX must be a scalar'

   lags = - lindgen(zmax-zmin+1) + pixoffset - long(zmin) ; must be integers
   nlag = n_elements(lags)

   chi2arr = fltarr(nobj, nlag)
   dofarr = fltarr(nobj, nlag)

   ;---------------------------------------------------------------------------
   ; Create outuput structure

   zstruct1 = create_struct( $
    'z'    , 0.0, $
    'zerr' , 0.0, $
    'chi2' , 0.0, $
    'dof'  ,  0L, $
    'theta', fltarr(nstar) )
   zstruct = replicate(zstruct1, nobj, nfind)

   ;---------------------------------------------------------------------------

   sqivar = sqrt(objivar)
   bvec = objflux * sqivar
   objmask = objivar NE 0
   starmask = starivar NE 0

   for iobj=0, nobj-1 do begin
      for ilag=0, nlag-1 do begin
;junk = min(chi2arr[iobj,*], ilag)
         i1 = (lags[ilag] < npixstar-1)
         if (i1 LT 0) then i0 = -i1 $
          else i0 = 0
         i1 = i1 > 0
         i2 = npixstar-1 < (npixobj+i1-i0-1)
         ncomp = i2 - i1 + 1
         mmatrixt = starflux[i1:i2,*]
         for i=0, nstar-1 do $
          mmatrixt[*,i] = mmatrixt[*,i] $
           * sqivar[i0:i0+ncomp-1,iobj] * starmask[i1:i2]
         mmatrix = transpose( mmatrixt )
         mm = mmatrix # mmatrixt

         ; Use SVD to invert the matrix
;         mmi = invert(mm, /double)
         if (nstar EQ 1) then begin
            mmi = 1.0 / mm
         endif else begin
            svdc, double(mm), ww, uu, vv
            mmi = 0 * vv
            for i=0, nstar-1 do mmi[i,*] = vv[i,*] / ww[i]
            mmi = mmi ## transpose(uu)
         endelse

         acoeff = mmi # (mmatrix # (bvec[i0:i0+ncomp-1,iobj] * starmask[i1:i2]))
         chi2arr[iobj,ilag] = $
          total( (mmatrixt # acoeff $
           - bvec[i0:i0+ncomp-1,iobj] * starmask[i1:i2])^2 )
         dofarr[iobj,ilag] = $
          total(objmask[i0:i0+ncomp-1,iobj] * starmask[i1:i2]) - nstar
; yfit = transpose(acoeff # mmatrix) / sqivar[i0:i0+ncomp-1] ; ???
; yfit1 = transpose(acoeff[0] # mmatrix[0,*]) / sqivar[i0:i0+ncomp-1] ; ???
; splot, smooth(objflux[i0:i0+ncomp-1,iobj],7)
; soplot, yfit, color='red'
; soplot,yfit1,color='green'

         ; Burles counter of lag number...
         print, format='("Lag ",i5," of ",i5,a1,$)', $
          ilag, nlag, string(13b)
      endfor

if (nlag EQ 1) then $
 zstruct[iobj,0].theta = acoeff

   endfor

   ;-----
   ; Fit this chi2 minimum with a parabola

if (n_elements(lags) GT 1) then begin
   for iobj=0, nobj-1 do begin
      xpeak1 = find_npeaks(-chi2arr[iobj,*], lags, nfind=nfind, $
       minsep=3, width=3, ypeak=ypeak1, xerr=xerr1, npeak=npeak)
      zstruct[iobj,0:npeak-1].z = transpose([ zoffset - xpeak1 ])
      zstruct[iobj,0:npeak-1].zerr = transpose([ xerr1 ])
      zstruct[iobj,0:npeak-1].chi2 = transpose([ -ypeak1 ])
      for ipeak=0, npeak-1 do begin
         junk = min(abs(lags-xpeak1[ipeak]), ilag)
         zstruct[iobj,ipeak].dof = dofarr[iobj,ilag]
;         zstruct[iobj,ipeak].theta = ???
      endfor
   endfor
endif

   return, zstruct
end
;------------------------------------------------------------------------------
