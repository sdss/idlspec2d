
; Procedure to fit 'template.fit' or simluations of velocity dispersions

pro fitsimul, filename, star

for i=0,31 do begin 
 
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
  veldisp, b, bsig, star[*,0], star[*,1], result

  teststruct = { ANSWER, templateno: 0L, sigma:0.0, sn:0.0, measured:fltarr(16)}

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
      fulltest[count].measured = result[count+lindgen(16)*90].sigma_diff*70.0
      count = count + 1
    endfor
  endfor

  if (i EQ 0) then returntest = fulltest $
  else returntest = [returntest, fulltest]

endfor

  mwrfits, returntest, 'simulsample.fit'
  stop
return
end
