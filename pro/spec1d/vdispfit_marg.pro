;+
; NAME:
;   vdispfit_marg
;
; PURPOSE:
;   Compute velocity dispersions for objects using the marginalized likelihood functions
;
; CALLING SEQUENCE:
;   vdispfit_marg, vdisp_lnl, vdisp_marg=vdisp_marg, vdisp_err_marg=vdisp_err_marg
;
; INPUTS:
;   vdisp_lnl        - Natural logarithm of the likelihoods; array of [NSIG, NOBJECT]
;
; OUTPUTS:
;   vdisp_marg       - Velocity dispersions in km/s; array of [NOBJECT]
;   vdisp_err_marg   - Velocity dispersion errors in km/s; array of [NOBJECT]
;
; PROCEDURES CALLED:
;   find_nminima
;
; REVISION HISTORY:
;   Dec. 2012 Written by Yiping Shu, Utah
;------------------------------------------------------------------------------

pro vdispfit_marg, vdisp_lnl, vdisp_marg=vdisp_marg, vdisp_err_marg=vdisp_err_marg

   nobj=(size(vdisp_lnl))[2]
   vdisp_marg=fltarr(nobj)
   vdisp_err_marg=fltarr(nobj)

   ; Construct the velocity dispersion baseline, consistent with vdispfit.pro
   dsig=25.0
   bigsig=findgen(35)*dsig
   
   for iobj=0L, nobj-1 do begin
      chi2arr=-2.0*reform(vdisp_lnl[*, iobj]) ; Compute the chi square array from likelihood
      junk=min(chi2arr, imin)
      ; Use the 3 points nearest the minimum of the chi square curve to compute
      ; a more accurate velocity dispersion
      ; Following the same strategy as vdispfit.pro
      if (imin gt 0) then $
         sigma=find_nminima(chi2arr, bigsig, width=1.5*dsig, ypeak=minchi2, xerr=sig_err, $
         errcode=errcode) $
      else $
         sigma=find_nminima([chi2arr[1], chi2arr], [-bigsig[1], bigsig], width=1.5*dsig, $
         ypeak=minchi2, xerr=sig_err, errcode=errcode)
      if (sigma ge max(bigsig)) then sig_err=-3L
      vdisp_marg[iobj]=sigma > 0.0
      ; Set vdisp_err_marg to the error-code if it is non-zero
      vdisp_err_marg[iobj]=sig_err*(errcode EQ 0)+errcode
   endfor
end
