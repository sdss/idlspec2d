;------------------------------------------------------------------------------
pro psolve_addstars, image, psfimg, objs

   xcen0 = round(objs.xcen)
   ycen0 = round(objs.ycen)
   xcen1 = objs.xcen - xcen0
   ycen1 = objs.ycen - ycen0

   kdims = size(psfimg, /dimens)
   nhalfx = (kdims[0]-1)/2
   nhalfy = (kdims[1]-1)/2
   if (size(psfimg,/n_dimen) EQ 2) then nmap = 1 $
    else nmap = kdims[2]

   for istar=0L, n_elements(objs)-1 do begin
      xint = xcen0[istar]
      yint = ycen0[istar]
      for imap=0L, nmap-1 do begin
         thispsf = objs[istar].flux * objs[istar].pvalue[imap] $
          * sshift2d(psfimg[*,*,imap], [xcen1[istar], ycen1[istar]])
         image[xint-nhalfx:xint+nhalfx,yint-nhalfy:yint+nhalfy] = $
          image[xint-nhalfx:xint+nhalfx,yint-nhalfy:yint+nhalfy] + thispsf
      endfor
   endfor

   return
end
;------------------------------------------------------------------------------
