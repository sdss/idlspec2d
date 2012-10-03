;+
; NAME:
;   sn_median
;
; PURPOSE:
;   Compute median spectroscopic S/N in the SDSS filters
;
; CALLING SEQUENCE:
;   snmed = sn_median(loglam, objflux, objivar [, sn_all=sn_all])
;
; INPUTS:
;   loglam     - Log10 wavelength [NPIX]
;   objflux    - Flux [NPIX,NOBJ]
;   objivar    - Inverse variance in same units as OBJFLUX [NPIX,NOBJ]
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   snmed      - Median S/N in each of the SDSS filters per pix [5,NOBJ]
;
; OPTIONAL OUTPUTS:
;   sn_all     - Median S/N across all bands [NOBJ]
;
; COMMENTS:
;   The median S/N is computed for only good data points that fall
;   within some hard-wired wavelength ranges defined by those
;   wavelengths at which the SDSS filters have half the transmission
;   of the mean transmission between those wavelengths.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   16-Feb-2011  Written by D. Schlegel, LBL
;   03-Oct-2012  Added sn_all keyword (A. Bolton, Utah)
;-
;------------------------------------------------------------------------------
function sn_median, loglam, objflux, objivar, sn_all=sn_all

   w1 = [3280, 4030, 5600, 6900, 8200]
   w2 = [3850, 5330, 6760, 8160, 9860]

   w1all = min(w1)
   w2all = max(w2)

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = dims[1]
   snmed = fltarr(5,nobj)

   if arg_present(sn_all) then sn_all = fltarr(nobj)

   for iobj=0L, nobj-1L do begin
      for j=0, 4 do begin
         indx = where(loglam GE alog10(w1[j]) AND loglam LT alog10(w2[j]) $
          AND objivar[*,iobj] GT 0, ct)
         if (ct GT 1) then $
          snmed[j,iobj] = median( sqrt(objivar[indx,iobj])*objflux[indx,iobj] )
      endfor
      if arg_present(sn_all) then begin
         indx = where(loglam GE alog10(w1all) AND loglam LT alog10(w2all) $
          AND objivar[*,iobj] GT 0, ct)
         if (ct GT 1) then $
          sn_all[iobj] = median( sqrt(objivar[indx,iobj])*objflux[indx,iobj] )
      endif
   endfor

   return, snmed
end
;------------------------------------------------------------------------------
