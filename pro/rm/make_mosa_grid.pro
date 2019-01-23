;+
; NAME:
;
; PURPOSE:
;   Make dither grid for the Mayall MOSA observation of the RM field

pro make_mosa_grid, psplot=psplot, quad=quad, absoff=absoff

   ; This is the RM field center
   RA_cen=213.704D & DEC_cen=+53.08334D
   ; 14:14:49.0+53.05:00.0
   
   ; This is the pointing center of 4 quadrants
   RA_quad=dblarr(4) & DEC_quad=dblarr(4)
   RA_quad_off=[-0.75,-0.75,0.75,0.75]
   DEC_quad_off=[0.75,-0.75,-0.75,0.75]
   theta = DEC_cen + DEC_quad_off

   ; Compute the offset of each grid point (32 in total)
   ; from the central pointing (RM field center)
   ; Two dither positions (a/b) for each grid point
   RA_grid = dblarr(32, 2)
   DEC_grid = dblarr(32,2)
   RA_grid_abs = dblarr(32,2)
   DEC_grid_abs = dblarr(32,2)
   
   inext = 0
   d_off = 30.D ; arcsec, the chip gaps are 20.8"
   name = string(indgen(32)+1,format='(i2.2)')
   if not keyword_set(quad) then begin
      for irow=0L, 5L do begin

         for jcol=0L,5L do begin

            ra_off = (jcol-3L + 0.5)*0.5*3600.D
            dec_off = (-irow+3L - 0.5)*0.5*3600.D
         
            if irow eq 0 or irow eq 5 then begin ; top and bottom fields
               if jcol ne 0 and jcol ne 5 then begin
               
                  RA_grid[inext,0]=ra_off+d_off
                  RA_grid[inext,1]=ra_off-d_off
               
                  DEC_grid[inext,0]=dec_off+d_off
                  DEC_grid[inext,1]=dec_off-d_off
               
                  inext = inext+1L
               endif
            
            endif else begin
         
               RA_grid[inext,0]=ra_off+d_off
               RA_grid[inext,1]=ra_off-d_off
               
               DEC_grid[inext,0]=dec_off+d_off
               DEC_grid[inext,1]=dec_off-d_off   
         
               inext = inext+1L 
            endelse
            
         endfor
      endfor
   endif else begin ; divide the 3x3deg square into 4 quadrants
   
      for i=0L, 3L do begin
         ; First find the pointing center of this quadrant   
         RA_quad[i] = RA_quad_off[i]/cos(theta[i]/180.*!PI) + RA_cen
         DEC_quad[i] = (DEC_quad_off[i] + DEC_cen)
         
         print, i+1, dec2hms(ra_quad[i]/15.,/double), ' ', dec2hms(dec_quad[i],/double)
      
         if i eq 0 or i eq 1 then sign1 = 1 else sign1 = -1
         if i eq 1 or i eq 2 then sign2 = -1 else sign2 = 1
         
         case i of
               0: begin
                  step_ra=[0,1,1,1,0,-1,-1,0]
                  step_dec=[1,1,0,-1,-1,-1,0,0]
               end
               1: begin
                  step_ra=[0,1,1,1,0,-1,-1,0]
                  step_dec=[-1,-1,0,1,1,1,0,0]
               end
               2: begin
                  step_ra=[0,-1,-1,-1,0,1,1,0]
                  step_dec=[-1,-1,0,1,1,1,0,0]
               end
               3: begin
                  step_ra=[0,-1,-1,-1,0,1,1,0]
                  step_dec=[1,1,0,-1,-1,-1,0,0]
               end
         endcase
         ; now assign the 8 sub-points in each quadrant
         for j=0L,7L do begin
    
           RA_grid[inext,0]=step_ra[j]*0.5*3600.D + d_off
           RA_grid[inext,1]=step_ra[j]*0.5*3600.D - d_off
           DEC_grid[inext,0]=step_dec[j]*0.5*3600.D + d_off
           DEC_grid[inext,1]=step_dec[j]*0.5*3600.D - d_off
           ;print,inext+1,RA_grid[inext,0],DEC_grid[inext,0]
      
           ; assign absolute offset from the previous position
           if j eq 0 then begin
              ra_grid_abs[inext,0]=ra_grid[inext,0]
              ra_grid_abs[inext,1]= -2.*d_off
              dec_grid_abs[inext,0]=dec_grid[inext,0]
              dec_grid_abs[inext,1]= -2.*d_off
           endif else begin
              ra_grid_abs[inext,0]=ra_grid[inext,0]-ra_grid[inext-1,1]
              ra_grid_abs[inext,1]= -2.*d_off
              dec_grid_abs[inext,0]=dec_grid[inext,0]-dec_grid[inext-1,1]
              dec_grid_abs[inext,1]= -2*d_off
              
           endelse
      
           inext = inext+1L
         endfor
        
      endfor
   
   endelse
   
   ; Read in the targets   
   file=getenv('IDLRM_DIR')+'/etc/RM_targets.fits'
   target=mrdfits(file,1)
   ind=where(target.tiled eq 1)
   target=target[ind]
   obj_raoff = (target.ra - ra_cen)*cos(target.dec/180.*!PI)
   obj_decoff = (target.dec - dec_cen)
   
   if keyword_set(psplot) then begin
     figfile=getenv('IDLRM_DIR')+'/etc/mosa_offsets.eps'
     if keyword_set(quad) then figfile=getenv('IDLRM_DIR')+'/etc/mosa_offsets_quad.eps'
     begplot,name=figfile,/color,/encap,xsize=8,ysize=8,/cmyk
   endif
   
   plot, [0],[0],xrange=[1.6,-1.6],yrange=[-1.6,1.6],/xsty,/ysty $
      , xtitle='RA Offset [deg]', ytitle='DEC Offset [deg]'
   
   oplot, obj_raoff,obj_decoff,psym=1
   
   box=2130. ; arcsec; 8192x0.26"
   color1 = fsc_color('red')
   color2 = fsc_color('blue')
   textcolor=fsc_color(['black','blue','cyan','red'])
   
   for i=0L, 31L do begin
      
      if not keyword_set(quad) then begin
         add_off_ra = 0.D & add_off_dec = 0.D
      endif else begin
         add_off_ra = (RA_quad[i/8]-RA_cen)*cos(theta[i/8]/180.*!PI)
         add_off_dec = DEC_quad[i/8]-DEC_cen
         ;print, add_off_ra, add_off_dec
      endelse
      
      oplot, [ra_grid[i,0]-0.5*box,ra_grid[i,0]-0.5*box]/3600.D +add_off_ra, $
             [dec_grid[i,0]-0.5*box,dec_grid[i,0]+0.5*box]/3600.D +add_off_dec,color=color1
      oplot, [ra_grid[i,0]+0.5*box,ra_grid[i,0]+0.5*box]/3600.D +add_off_ra, $
             [dec_grid[i,0]-0.5*box,dec_grid[i,0]+0.5*box]/3600.D +add_off_dec,color=color1
      oplot, [ra_grid[i,0]-0.5*box,ra_grid[i,0]+0.5*box]/3600.D +add_off_ra, $
             [dec_grid[i,0]-0.5*box,dec_grid[i,0]-0.5*box]/3600.D +add_off_dec,color=color1
      oplot, [ra_grid[i,0]-0.5*box,ra_grid[i,0]+0.5*box]/3600.D +add_off_ra, $
             [dec_grid[i,0]+0.5*box,dec_grid[i,0]+0.5*box]/3600.D +add_off_dec,color=color1
   
      oplot, [ra_grid[i,1]-0.5*box,ra_grid[i,1]-0.5*box]/3600.D +add_off_ra, $
             [dec_grid[i,1]-0.5*box,dec_grid[i,1]+0.5*box]/3600.D +add_off_dec,color=color2
      oplot, [ra_grid[i,1]+0.5*box,ra_grid[i,1]+0.5*box]/3600.D +add_off_ra, $
             [dec_grid[i,1]-0.5*box,dec_grid[i,1]+0.5*box]/3600.D +add_off_dec,color=color2
      oplot, [ra_grid[i,1]-0.5*box,ra_grid[i,1]+0.5*box]/3600.D +add_off_ra, $
             [dec_grid[i,1]-0.5*box,dec_grid[i,1]-0.5*box]/3600.D +add_off_dec,color=color2
      oplot, [ra_grid[i,1]-0.5*box,ra_grid[i,1]+0.5*box]/3600.D +add_off_ra, $
             [dec_grid[i,1]+0.5*box,dec_grid[i,1]+0.5*box]/3600.D +add_off_dec,color=color2
   
      xyouts, 0.5*(ra_grid[i,0]+ra_grid[i,1])/3600.D + add_off_ra, $
         0.5*(dec_grid[i,0]+dec_grid[i,1])/3600.D + add_off_dec, name[i], $
         align=0.5, color=textcolor[i/8]
   endfor
   
   if keyword_set(psplot) then endplot
   
   ;NOW output an ascii file with one absolute offset (in arcsec) per line
   filename = getenv('IDLRM_DIR')+'/etc/mosa_offsets'
   if keyword_set(quad) then filename = getenv('IDLRM_DIR')+'/etc/mosa_offsets_quad'
   if keyword_set(absoff) then filename = filename + '_absoff'
   openw,lun,filename,/get_lun
   
   fmt='(A0, " ", i5, " ", i5, " ", f7.1)'
   fmt1='(A0, " ", i5, " ", i5, " ", f8.3, " ", f8.4)'
   for i=0L, 31L do begin
   
      ; If require absolute offset, then subtract the previous offset
      if keyword_set(absoff) then begin
         ra_grid = ra_grid_abs & dec_grid = dec_grid_abs
      endif
   
      if not keyword_set(quad) then begin
         printf, lun, format=fmt, 'RM'+string(i+1,format='(i2.2)')+'a', $
             ra_grid[i,0], dec_grid[i,0], 2000.
         printf, lun, format=fmt, 'RM'+string(i+1,format='(i2.2)')+'b', $
             ra_grid[i,1], dec_grid[i,1], 2000.  
       endif else begin       
         printf, lun, format=fmt1, 'RM'+string(i+1,format='(i2.2)')+'a', $
             ra_grid[i,0], dec_grid[i,0], ra_quad[i/8],dec_quad[i/8]
         printf, lun, format=fmt1, 'RM'+string(i+1,format='(i2.2)')+'b', $
             ra_grid[i,1], dec_grid[i,1], ra_quad[i/8],dec_quad[i/8]
       endelse
   endfor
   
   close,lun
   free_lun,lun
   
end