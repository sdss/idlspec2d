; Procedure to fit 'template.fit' or simluations of velocity dispersions

pro fitsimul, filename, starflux, starsig, starwave,nn0,nn, czmin, czmax, returntest


for i=nn0,nn do begin 
  
  a = mrdfits(filename, i)
  asig = a*0.0
  meanstar = mean(a)

  sn = 0.0
  for n=0,8 do begin
    sn = sn + 10.0
    asig[*,*,n,*] = meanstar/sn
  endfor

  b = reform(a, 3918, 10*9*16)
  bsig = reform(asig, 3918, 10*9*16)  
  bwave= b*0
  nstar= (size(b))[2]
  for k=0,nstar-1 do begin 
  bwave[*,k] = starwave[*,i]
  endfor 
  	
  bredshifts=0*findgen((size(b))[2])	

  mveldisp, b, bsig, bwave, starflux[*,i], starsig[*,i], starwave[*,i], result, redshifts=bredshifts, czmin= czmin, czmax= czmax, klo_cut=0.016, khi_cut=0.23, maxsig=6,sigmastep=0.4


  teststruct = { ANSWER, templateno: 0L, sigma:0.0, sn:0.0, $
                 measured_cc:fltarr(16), measured_cc_err:fltarr(16), $
                 measured_quot_diff:fltarr(16), measured_quot_diff_err:fltarr(16), $
                 measured_quotient:fltarr(16), measured_quotient_err:fltarr(16),$
                 measured_diff:fltarr(16), measured_diff_err:fltarr(16)}

  fulltest = replicate(teststruct,90)

  fulltest.templateno = i
  count = 0
  sn = 0.0
  for n=0,8 do begin
    sn = sn+ 10.0
    res = 80.0
    for nn=0,9 do begin
      res = res + 20.0
      fulltest[count].sigma = res
      fulltest[count].sn = sn
      fulltest[count].measured_cc = result[count+lindgen(16)*90].sigma_cc
      fulltest[count].measured_cc_err = result[count+lindgen(16)*90].sigma_cc_err
      fulltest[count].measured_quot_diff = result[count+lindgen(16)*90].sigma_quot_diff
      fulltest[count].measured_quot_diff_err = result[count+lindgen(16)*90].sigma_quot_diff_err
      fulltest[count].measured_quotient = result[count+lindgen(16)*90].sigma_quotient
      fulltest[count].measured_quotient_err = result[count+lindgen(16)*90].sigma_quotient_err
      fulltest[count].measured_diff = result[count+lindgen(16)*90].sigma_diff
      fulltest[count].measured_diff_err = result[count+lindgen(16)*90].sigma_diff_err
      count = count + 1
    endfor
  endfor 

if (size(returntest,/tname) EQ 'STRUCT') then begin
returntest = [returntest, fulltest]
endif else begin 
returntest = fulltest 
endelse

endfor

nametest=string('veldisp_simul_aug00_'+string(nn0,format='(i2.2)')+$
                 '_'+string(nn,format='(i2.2)')+'.fit')
mwrfits,returntest, nametest


return
end
