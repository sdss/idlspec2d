pro monster, loglam, skyimage

     nsky = (size(skyimage))[2]
     image = skyimage

     for j=0,nsky-1 do image[*,j] = median(image[*,j],21)
     
     filter = exp(-0.5*((loglam - 3.814)/0.004)^2)
       
     con = total(filter * image,1)/total(filter,1)
     con2 = total(shift(filter,300) * image,1)/total(filter,1)

     djs_iterstat, con, median=median, sigma=sigma
     djs_iterstat, con2, median=median, sigma=sigma2

;     printf, 1, strmid(allrframes[i],12), mjd, ' ', name, median, sigma, $
;         sigma2, format='(a,i6,a,a,3f8.3)'

;     print, strmid(allrframes[i],12), mjd, ' ', name, median, sigma, sigma2, $
;          format='(a,i6,a,a,3f8.3)'

      if sigma GT 3 * sigma2 then $
        splog, 'WARNING: Red Monster has been spotted' , sigma/sigma2

 return
end

