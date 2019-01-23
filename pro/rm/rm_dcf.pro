;+
; NAME:
;   rm_dcf
;
; PURPOSE:
;   Compute the discrete correlation function (Edelson & Krolik 1988)
;
; CALLING SEQUENCE:
;   rm_dcf,t1,LC1,t2,LC2
;
; INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; REVISION HISTORY:
;   28-Apr-2014  Written by Yue Shen, Carnegie
;-
;------------------------------------------------------------------------
pro rm_dcf, t1, LC1, t2, LC2, err1=err1,err2=err2,tau_grid=tau_grid, $
     rm_corr_err=rm_corr_err, DCF=DCF,eDCF=eDCF,UDCFij=UDCFij, $
     npair=npair,indp_npair=indp_npair

   if not keyword_set(tau_grid) then tau_grid=-50.+findgen(41)*5
   dtau=median(tau_grid-shift(tau_grid,1))
   ngrid=n_elements(tau_grid)   

   ; Default is to remove data observed on the same day (i.e., with correlated errors)
   if n_elements(rm_corr_err) eq 0 then rm_corr_err=1

   n1=n_elements(LC1) & n2=n_elements(LC2)

   if n_elements(err1) eq 0 then err1=0.
   if n_elements(err2) eq 0 then err2=0.
   mean1=mean(LC1) & sig1=stddev(LC1)
   mean2=mean(LC2) & sig2=stddev(LC2)

   ; Collect the pairs [eq.3 of EK 1988]
   UDCFij = ( (LC1 - mean1) # (LC2 - mean2) ) / sqrt((sig1^2-err1^2)*(sig2^2-err2^2))
   dtij = dblarr(n1,n2)
   for i=0L, n2-1 do dtij[*,i] = t2[i] - t1

   ; Bin the pairs
   DCF=dblarr(ngrid) & eDCF=dblarr(ngrid)-1.
   npair=lonarr(ngrid) & indp_npair=lonarr(ngrid)
   for i=0L, ngrid-1 do begin

      if keyword_set(rm_corr_err) then $
       ind = where(dtij ge tau_grid[i]-0.5*dtau and dtij lt tau_grid[i]+0.5*dtau $
         and abs(dtij) gt 1d-5, mm) else $
       ind = where(dtij ge tau_grid[i]-0.5*dtau and dtij lt tau_grid[i]+0.5*dtau, mm)
      npair[i]=mm

      if mm gt 0 then begin

         DCF[i] = total(UDCFij[ind])/double(mm)

         ; Determine the independent data points in LC1
         indd_col = ind mod n1
         ttt = uniq(indd_col, sort(indd_col))
         mm1 = n_elements(ttt)
         indp_npair[i]=mm1
        
         ; Estimate the error of the DCF in this bin
         if mm1 gt 1 then $
         eDCF[i] = sqrt(total( (UDCFij[ind] - DCF[i])^2 ) ) / sqrt( (mm - 1.)*(mm1 - 1.)  )

      endif 


   endfor

end
