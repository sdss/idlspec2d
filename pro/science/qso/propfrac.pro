
;
;   Calculate the line of sight proper distance between z1 and z2
;   with LambdaCDM (0.3,0.7) and Standard CDM cosmology (1.0,0.0)
;
pro propfrac, z1, z2, ldist=ldist, sdist=sdist

;
;  These two are proper distances good to z=10
;
     lambdaCDM = [   0.0013107657d, 2.2809960d, -0.87851526d, -3.5277213d, $
                     6.1501481d, -4.2997680d, 1.3070295d, -0.10737073d]

     standardCDM = [   0.0010117448d, 2.2989689d, -3.9564659d, 4.4823482d, $
                     -3.6410577d, 2.0945059d, -0.76402620d, 0.13079970d]


     lambda_dist1 = poly(alog10(1+z1), lambdaCDM)
     standard_dist1 = poly(alog10(1+z1), standardCDM)
   
     lambda_dist2 = poly(alog10(1+z2), lambdaCDM)
     standard_dist2 = poly(alog10(1+z2), standardCDM)
  
     x = lambda_dist1/lambda_dist2
     ldist = x * (1-x)
   
     x = standard_dist1/standard_dist2
     sdist = x * (1-x)
  
     return
end 

