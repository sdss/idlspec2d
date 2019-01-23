; given # of panels, return the xypositions in the plot

pro plot_layout, npanel, xypos=xypos, omargin=omargin, pmargin=pmargin, $
     nrow=nrow

   if ~keyword_set(omargin) then omargin=[0.1, 0.1, 0.98, 0.98] ; overall page margin
   if ~keyword_set(pmargin) then pmargin=[0.02,0.04]

   if npanel eq 1 then xypos=omargin else begin

      pmargin_x=pmargin[0]  ; margin between panels
      pmargin_y=pmargin[1]
      xtot=omargin[2]-omargin[0] & ytot=omargin[3]-omargin[1]
      if ~keyword_set(nrow) then nrow=round(sqrt(npanel))
      ncol=ceil(float(npanel)/nrow)
      dx0=xtot/ncol & dy0=ytot/nrow
      dx=(xtot - (ncol-1)*pmargin_x)/ncol
      dy=(ytot - (nrow-1)*pmargin_y)/nrow
      xypos=dblarr(4,npanel)

      for j=0, npanel-1 do begin
         xypos[0,j]=omargin[0] + (j mod ncol)*dx0
         xypos[2,j]=xypos[0,j]+dx
         xypos[1,j]=omargin[3] - (j/ncol)*dy0 - dy
         xypos[3,j]=xypos[1,j]+dy
      endfor
   endelse
end
