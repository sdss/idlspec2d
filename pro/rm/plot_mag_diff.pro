;+
; NAME:
;  plot_mag_diff
;
; PURPOSE:
;  Plot the difference between synthetic and calibflux magnitudes
;   for different RM epochs
;---------------------------------------------------------------------------

pro plot_mag_diff, shuffle=shuffle, plotdata=plotdata

   if n_elements(plotdata) eq 0 then plotdata=0L
   if n_elements(shuffle) eq 0 then shuffle = 1L

   platelist = [7338,7338,7338,7339,7339]
   mjdlist = [56660,56664,56669,56683,56686]
   seeing = [1.6,3.,3.,1.6,2.]
   sn2 = [17.9,10.9,3.3,22.4,21.9]
   ncoadd = [7,6,6,9,8]
   nepoch = n_elements(platelist)

   std_fiber = rm_findstd(platelist,mjdlist,synflux=synflux,calibflux=calibflux)

   ; get airmass info
   airmass = dblarr(nepoch)
   for i=0L, nepoch - 1L do begin
     rm_expinfo,platelist[i],mjdlist[i],expinfo=expinfo
     ; keep only coadded exposures
     ind = where(expinfo.coadded[0,*] eq 1)
     expinfo = expinfo[ind]
     airmass[i] = mean(expinfo.airmass)
   endfor


   plotfile = 'std_specphoto_dep.ps'
   set_plot, 'ps'
   plotfile = getenv('IDLRM_DIR')+'/misc/'+plotfile
   begplot,name=plotfile,/color,/landscape

   !p.multi = [0,3,3]
   ;multipanel,row=3,col=3,/first,position=p,margin=[0.05,0.05,0.05,0.05]

   colors = fsc_color(['opposite','green','red','cyan','blue','magenta'])
   csize=2.0 & symsize=0.1 & csize2 = 0.8

   ; first plot synmag-calibmag as function of SN^2
   magdiff = -alog10(synflux/calibflux)*2.5
   title = ['u band','g band','r band','i band','z band']
   yrange = [-0.15, 0.15]
   string2 = string(ncoadd,format='(i0)')
   string2[0]='exp:'+string2[0]
   for iband=1,3 do begin

     if iband eq 1 then noerase = 0 else noerase = 1

     plot,[0],[0],/nodata,xrange=[0,25],yrange=yrange,/xsty,/ysty, $
       xtitle = 'SN2', ytitle='SynMag - CalibMag',charsize=csize, $
       title = title[iband]

     ;multipanel, /advance, pos=p

     xpos = 1 & ypos = 0.12
     for iepoch=0,nepoch - 1L do begin

        ind = where(synflux[iepoch,*,iband] gt 0, ngood)
        magdiff1 = magdiff[iepoch,ind,iband]
        xarr = replicate(sn2[iepoch],ngood)
        if keyword_set(shuffle) then xarr = xarr + randomu(seed,ngood)*1 - 0.5
        if keyword_set(plotdata) then oplot, xarr, magdiff1,psym=2, $
          symsize=symsize, color=colors[iepoch]
        oploterror,[sn2[iepoch]],[mean(magdiff1)],[stddev(magdiff1)],$
          psym=symcat(6),color=colors[iepoch],errcolor=colors[iepoch]
        oploterror,[sn2[iepoch]],[mean(magdiff1)],[median(abs(magdiff1))],$
          psym=symcat(6),color=colors[iepoch],errcolor=colors[iepoch],errstyle=1

        if iband eq 1 then begin
           xyouts, xpos, ypos, string(mjdlist[iepoch],format='(i0)'),$
             color=colors[iepoch], charsize=csize2
           xyouts, xpos, -ypos, string2[iepoch],$
             color=colors[iepoch], charsize=csize2
           xpos = xpos + 4.5
        endif
     endfor     
   endfor

   ; then plot synmag - calibmag as function of seeing
   for iband=1,3 do begin

     plot,[0],[0],/nodata,xrange=[1,3.5],yrange=yrange,/xsty,/ysty, $
       xtitle = 'seeing', ytitle='SynMag - CalibMag',charsize=csize, $
       title = title[iband]  ;, pos=p, /noerase

     ;multipanel, /advance, pos=p

     for iepoch=0,nepoch - 1L do begin

        ind = where(synflux[iepoch,*,iband] gt 0, ngood)
        magdiff1 = magdiff[iepoch,ind,iband]
        xarr = replicate(seeing[iepoch],ngood) 
        if keyword_set(shuffle) then xarr = xarr + randomu(seed,ngood)*0.1 - 0.5*0.1
        if keyword_set(plotdata) then oplot, xarr, magdiff1,psym=2, $
          symsize=symsize, color=colors[iepoch]
        oploterror,[seeing[iepoch]],[mean(magdiff1)],[stddev(magdiff1)],$
          psym=symcat(6),color=colors[iepoch],errcolor=colors[iepoch]
        oploterror,[seeing[iepoch]],[mean(magdiff1)],[median(abs(magdiff1))],$
          psym=symcat(6),color=colors[iepoch],errcolor=colors[iepoch],errstyle=1

     endfor
   endfor

   ; then plot synmag - calibmag as function of airmass
    for iband=1,3 do begin

     plot,[0],[0],/nodata,xrange=[1.05,1.3],yrange=yrange,/xsty,/ysty, $
       xtitle = 'avg airmass', ytitle='SynMag - CalibMag',charsize=csize, $
       title = title[iband] ;, pos=p, /noerase

     ;multipanel, /advance, pos=p

     for iepoch=0,nepoch - 1L do begin

        ind = where(synflux[iepoch,*,iband] gt 0, ngood)
        magdiff1 = magdiff[iepoch,ind,iband]
        xarr=replicate(airmass[iepoch],ngood)
        if keyword_set(shuffle) then xarr = xarr + randomu(seed,ngood)*0.01 - 0.5*0.01
        if keyword_set(plotdata) then oplot, xarr, magdiff1,psym=2, $
          symsize=symsize, color=colors[iepoch]
        oploterror,[airmass[iepoch]],[mean(magdiff1)],[stddev(magdiff1)],$
          psym=symcat(6),color=colors[iepoch],errcolor=colors[iepoch]
        oploterror,[airmass[iepoch]],[mean(magdiff1)],[median(abs(magdiff1))],$
          psym=symcat(6),color=colors[iepoch],errcolor=colors[iepoch],errstyle=1

     endfor
   endfor


   ;multipanel, /off
   !p.multi = 0
   endplot

end
