;+
; NAME:
;  plot_exp_mratio
;  rm_checkmratio
;  rm_checkmratio_1exp
;  rm_diag_mratio_epoch
; PURPOSE:
;    Check the mratio distribution of good std for different Epochs and its dependence
;      on properties such as airmass, focal plane position, etc.
;
;
pro plot_exp_mratio, expname, calibdir=calibdir,overplot=overplot $
       , ratioscat=ratioscat, loglam=loglam, _extra=extra $
       , plotmratio=plotmratio, noplot=noplot, xrange = xrange $
       , channel=channel, perdiff=perdiff, nprox=nprox,selfexclud=selfexclud $
       , airmass=airmass, xyfit=xyfit, xfocal=xfocal, yfocal=yfocal
       
  ; plot the distribution of mratio in a given exposure
  ; only use good fluxing stars

   if not keyword_set(calibdir) then $
      calibdir='/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/recalib/test20/'

   filename=calibdir+expname
   kindx=mrdfits(filename,2,/silent)
    
   ; get mratio and flatarr for the current exposure
   struct_out=mrdfits(filename,4,/silent)
   
   ; is this a blue ccd?
   if strmatch(expname,'*-b*') then begin
      xrange=[3600,6300]
      yrange=[0,3]
      channel = 'blue'
   endif
   ; is this a red ccd?
   if strmatch(expname,'*-r*') then begin
      xrange=[6300,10000.]
      yrange=[0,3]
      channel = 'red'
   endif

   
   if not keyword_set(overplot) and not keyword_Set(noplot) then begin
      plot,[0],[0],/nodata,xrange=xrange,yrange=yrange,/xsty $
         ,title=expname    
   endif
   
   
   thisloglam=reform(struct_out.thisloglam)
   thismratio=reform(struct_out.thismratio)
   thisflatarr=reform(struct_out.thisflatarr)
   
   dims=size(thisloglam)
   npix=dims[1] & nobj=dims[2]
   
   
   for i=0L, nobj-1L do begin
   
      ind_norm=where(10.^thisloglam[*,i] gt 5000. and $
         10.^thisloglam[*,i] le 5100.)
      if not keyword_set(noplot) then begin
      
         if keyword_Set(plotmratio) then $   
          oplot, 10.^thisloglam[*,i],thismratio[*,i]/median(thismratio[ind_norm,i]) $
            , psym=3,_extra=extra
         
         ; oplot, 10.^thisloglam[*,i],thisflatarr[*,i]/median(thisflatarr[ind_norm,i]) $
         ;   , _extra=extra
         oplot, 10.^thisloglam[*,i],thisflatarr[*,i], _extra=extra
      endif 
   endfor
   
   ; compute the scatter in flatarr for all std stars in this exposure
   npix_new = floor((alog10(xrange[1])-alog10(xrange[0]))/1d-4)
   loglam = alog10(xrange[0]) + findgen(npix_new)*1d-4
   ratioscat = dblarr(npix_new)
   perdiff = dblarr(npix_new,nobj)

   flatarr = fltarr(npix_new,nobj)
   for i=0L,nobj-1L do flatarr[*,i]=interpol(thisflatarr[*,i], thisloglam[*,i], loglam)

   ; Using the nearest neighbors method to recompute flatarr
   if keyword_set(nprox) or keyword_set(xyfit) then begin

      flatarr_new = 0*flatarr

      if keyword_set(nprox) then begin
         nprox_ori = nprox
         nprox = nprox < nobj
   
         ; Here is the trick we do the localization for each std
         for i=0L,nobj-1L do begin
            ; Use the input ngood airmass, and find the nearest nprox std in airmass
            air_diff = abs(airmass - airmass[i])
            if (not keyword_set(selfexclud)) then $
              ind_prox = (sort(air_diff))[0:nprox-1] $
            else begin
              ind_prox = (sort(air_diff))[1:nprox]
            endelse

            if nprox ge 2 then flatarr_new[*,i] = mean(flatarr[*,ind_prox],dim=2) $
              else flatarr_new[*,i] = flatarr[*,ind_prox]

         endfor
         nprox = nprox_ori
      endif

      if keyword_set(xyfit) then begin
         ; at each wavelength pixel, fit a B-spline on xfocal, using yfocal as the 2nd dimension 
         everyn = nobj/3
         isort = sort(xfocal)
         nord = 3
         bkpt = 0
         fullbkpt = bspline_bkpts(xfocal[isort],everyn=everyn,bkpt=bkpt, nord=nord)
         for ipix=0L, npix_new - 1L do begin

            ; Fit a low-order B-spline to xfocal,using yfocal as a 2nd parameter
            ; This method has some issues with objects not sampled by the std locations
            ;sset = bspline_iterfit(xfocal[isort],reform(flatarr[ipix,isort]),lower=0, upper=0 $
            ;   , invvar=replicate(1,nobj) $
            ;   , fullbkpt=fullbkpt,outmask=outmask1,nord=nord,x2=yfocal[isort],npoly=2 $
            ;   , requiren=(everyn-1)>1)
            ;flatarr_new[ipix,isort] = bspline_valu(xfocal[isort], sset, x2=yfocal[isort])
 
            ; this is the actual xyfit method used
            pp_guess = [median(flatarr[ipix,*]), 0.,0.,0.,0. ]
            xarr = [[xfocal],[yfocal]]
            xyfit_para = mpfitfun('xy_polyfit',xarr, flatarr[ipix,*], $
                  replicate(1.,nobj), pp_guess, /quiet, status=status)
            ;splog, 'xy_polyfit Status=',status, ' PP=', xyfit_para
            flatarr_new[ipix,*] = xy_polyfit(xarr, xyfit_para)
           
            ;yrange1=[0.5,2.5] 
            ;plot, xfocal[isort], flatarr[ipix,isort], psym=5,yrange=yrange1
            
            ;oplot, xfocal[isort], bspline_valu(xfocal[isort], sset, x2=yfocal[isort]) $
            ;   ,color=fsc_color('red') 
            ;xyouts, -200, 1.8, 10.^loglam[ipix],charsize=1.5
            ;xarr = -300. + findgen(301)*2
            ;yarr = replicate(-200., 300)
            ;oplot,xarr,bspline_valu(xarr,sset,x2=yarr),color=fsc_color('cyan')
            ;message, 'stop'
            ;pause

         endfor

      endif

   endif

   ; normalize flatarr
   if (not keyword_set(nprox)) and (not keyword_set(xyfit)) then begin
      mean_flatarr = mean(flatarr,dim=2)
      flatarr_norm = 0*flatarr
      for i=0L,nobj-1L do $
       flatarr_norm[*,i] = mean_flatarr/flatarr[*,i]
   endif else begin
      flatarr_norm = 0*flatarr
      for i=0L,nobj-1L do $
       flatarr_norm[*,i] = flatarr_new[*,i]/flatarr[*,i]
   endelse

   perdiff = flatarr_norm - 1.

   for ipix=0L, npix_new - 1L do begin
   
      ratioscat[ipix] = stddev( flatarr_norm[ipix,*] )
      
   endfor
   
end

pro rm_checkmratio, calibname, calibdir=calibdir, topdir=topdir,title=title,nprox=nprox,$
     selfexclud=selfexclud,xyfit=xyfit

   if not keyword_set(calibdir) then $
      calibdir='/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/recalib/test5/'

   if not keyword_set(topdir) then topdir='/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/'

   nexp=n_elements(calibname)
   airmass=dblarr(nexp)
   sn=dblarr(nexp)

   for i=0L,nexp - 1L do begin
      spname=repstr(calibname[i],'Fluxcalib','Frame')
    
      print, topdir+spname, calibdir+calibname[i]
      kindx=mrdfits(calibdir+calibname[i],2,/silent)
      indx=where(kindx.qgood eq 1)

      if strmatch(calibname[i],'*1-*') then ccd=1 else ccd=2
      if ccd eq 1 then igood=kindx[indx].fiberid - 1L $
       else igood=kindx[indx].fiberid - 500L - 1L
      
      spframe_read,topdir+spname,igood,plugmap=plugmap,hdr=hdr1,objflux=objflux,objivar=objivar
      sn_all = objflux*sqrt(objivar)
      sn[i] = median(sn_all)
      get_tai, hdr1, tai_beg, tai, tai_end
      airmass_all = tai2airmass(plugmap.ra, plugmap.dec, tai=tai)
      xfocal = plugmap.xfocal & yfocal = plugmap.yfocal
      airmass[i] = mean(airmass_all) 
 
      expname=calibname[i]
      plot_exp_mratio,expname,calibdir=calibdir,ratioscat=ratioscat, loglam=loglam $
          , /noplot, xrange=xrange, channel=channel
   
      if n_elements(final_ratio) eq 0 then final_ratio=ratioscat $
        else final_ratio = [[final_ratio], [ratioscat]]

      if keyword_set(nprox) then begin ; the nprix nearest neighbor approach
         plot_exp_mratio,expname,calibdir=calibdir,ratioscat=ratioscat_nprox, loglam=loglam $
          , /noplot, xrange=xrange, channel=channel,airmass=airmass_all,nprox=nprox,selfexclud=selfexclud
         if n_elements(final_ratio_nprox) eq 0 then final_ratio_nprox=ratioscat_nprox $
           else final_ratio_nprox = [[final_ratio_nprox], [ratioscat_nprox]]
      endif
  
      if keyword_set(xyfit) then begin; fit positional-dependent polynomial to the stars
         plot_exp_mratio,expname,calibdir=calibdir,ratioscat=ratioscat_nprox, loglam=loglam $
          , /noplot, xrange=xrange, channel=channel,airmass=airmass_all,/xyfit $
          , xfocal=xfocal,yfocal=yfocal
         if n_elements(final_ratio_nprox) eq 0 then final_ratio_nprox=ratioscat_nprox $
           else final_ratio_nprox = [[final_ratio_nprox], [ratioscat_nprox]] 
      endif

   endfor

   min=min(airmass,max=max)
   
   colors=cgcolor(['opposite','red','green','cyan','blue','magenta','gold','dark green','Olive','dark gray'])
   yrange=[0,0.3]
   plot,[0],[0],xrange=xrange,yrange=yrange,/xsty,/ysty,/nodata, $
     xtitle='Wavelength',ytitle='normalized stddev in fluxing vector [100%]', $
     title=title
  
   for i=0L, nexp - 1L do begin 
      oplot, 10.^loglam,final_ratio[*,i],color=colors[i]
      if keyword_set(nprox) or keyword_set(xyfit) then $
          oplot,10.^loglam,final_ratio_nprox[*,i],color=colors[i],line=2
   endfor

   xyouts, xrange[0]+(xrange[1]-xrange[0])*0.1, 0.9*yrange[1],$
     channel+'cam '+string(ccd,format='(i0)'),charsize=1.0

   if keyword_Set(nprox) then xyouts, xrange[0]+(xrange[1]-xrange[0])*0.1, 0.08*yrange[1],$
     'Nprox='+string(nprox,format='(i0)'),charsize=1.0
   if keyword_set(xyfit) then begin
     ; xyouts, xrange[0]+(xrange[1]-xrange[0])*0.1, 0.08*yrange[1],$
     ;  'xyfit',charsize=1.0
     items=['BOSS pipeline', 'custom']
     legend, items, line=[0,2],box=0,pos=[xrange[0]+(xrange[1]-xrange[0])*0.1, 0.12*yrange[1]]
   endif

   xpos = xrange[0]+(xrange[1]-xrange[0])*0.7
   ypos = 0.9*yrange[1]
   xyouts, xpos, ypos, 'airmass',charsize=1.0
   isort = sort(airmass)
   items=string(airmass[isort],format='(f4.2)')
   legend, items,color=colors[isort],textcolor=colors[isort],pos=[xpos,ypos],box=0,charsize=1.0
   isort = sort(sn)
   items=string(sn[isort],format='(f4.1)') 
   xpos = xpos - (xrange[1]-xrange[0])*0.3
   xyouts, xpos, ypos, '    SN',charsize=1.0
   legend, items,color=colors[isort],textcolor=colors[isort],pos=[xpos,ypos],box=0,charsize=1.0  

   items=strmid(calibname,15,8)
   xpos = xpos - (xrange[1]-xrange[0])*0.35
   legend, items,color=colors,textcolor=colors,pos=[xpos,ypos],box=0,charsize=1.0

end

; this routine plot the actual mratio for each good std in a single exposure
pro rm_checkmratio_1exp, calibname, calibdir=calibdir, topdir=topdir,title=title, $
      noanno=noanno,nprox=nprox,selfexclud=selfexclud,xyfit=xyfit

   if not keyword_set(calibdir) then $
      calibdir='/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/recalib/test5/'

   if not keyword_set(topdir) then topdir='/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/'

   if n_elements(calibname) ne 1 then message, 'rm_checkmratio_1exp only works for a single exp'

   spname=repstr(calibname,'Fluxcalib','Frame')

   print, topdir+spname, calibdir+calibname
   kindx=mrdfits(calibdir+calibname,2,/silent)
   indx=where(kindx.qgood eq 1, ngood)

   if strmatch(calibname,'*1-*') then ccd=1 else ccd=2
   if ccd eq 1 then igood=kindx[indx].fiberid - 1L $
   else igood=kindx[indx].fiberid - 500L - 1L

   spframe_read,topdir+spname,igood,plugmap=plugmap,hdr=hdr1,objflux=objflux,objivar=objivar
   sn_all = objflux*sqrt(objivar)
   get_tai, hdr1, tai_beg, tai, tai_end
   ; airmass, xfocal and yfocal for ngood std
   airmass = tai2airmass(plugmap.ra, plugmap.dec, tai=tai)
   xfocal = plugmap.xfocal & yfocal = plugmap.yfocal

   expname=calibname
   plot_exp_mratio,expname,calibdir=calibdir,ratioscat=ratioscat, loglam=loglam $
     , /noplot, xrange=xrange, channel=channel, perdiff=perdiff

   if keyword_set(nprox) then $
     plot_exp_mratio,expname,calibdir=calibdir,ratioscat=ratioscat_nprox, loglam=loglam $
     , /noplot, xrange=xrange, channel=channel, perdiff=perdiff_nprox,airmass=airmass $
     , nprox=nprox,selfexclud=selfexclud

   if keyword_set(xyfit) then $
     plot_exp_mratio,expname,calibdir=calibdir,ratioscat=ratioscat_nprox, loglam=loglam $
     , /noplot, xrange=xrange, channel=channel, perdiff=perdiff_nprox,/xyfit $
     , xfocal=xfocal,yfocal=yfocal

   yrange=[-0.5,0.5]
   xmargin = !x.margin
   ; print, xmargin
   !x.margin = [7, 14]
   ang=string(197B)
   plot,[0],[0],xrange=xrange,yrange=yrange,/xsty,/ysty,/nodata, $
     xtitle='Wavelength ['+ang+']',ytitle='mratio deviation [100%]', $
     title=title

   cgloadct,25, ncolors=ngood
   colors = indgen(ngood)
   for i=0L, ngood - 1L do begin 
      oplot, 10.^loglam, perdiff[*,i],color=colors[i]
      if keyword_set(nprox) or keyword_set(xyfit) then $
         oplot, 10.^loglam, perdiff_nprox[*,i],color=colors[i],line=2
   endfor
   if keyword_Set(nprox) then xyouts, xrange[0]+(xrange[1]-xrange[0])*0.7, 0.8*yrange[0],$
     'Nprox='+string(nprox,format='(i0)'),charsize=1.0,color=cgcolor('opposite',256L)
   if keyword_set(xyfit) and strmatch(expname,'*-b*') then begin
      ; xyouts, xrange[0]+(xrange[1]-xrange[0])*0.7, 0.8*yrange[0],$
      ; 'xyfit',charsize=1.0,color=cgcolor('opposite',256L)
      items=['BOSS pipeline', 'custom']
      legend, items, line=[0,2],box=0,pos=[0.35,0.68],charsize=1.5, /norm, $
        color=cgcolor('opposite',256L),textcolor=cgcolor('opposite',256L),spacing=1,thick=5
   endif

   if not keyword_set(noanno) then begin
      xpos = 0.75 & ypos = 0.95
      xyouts, xpos, ypos, ' airmass',charsize=1.0,/norm, color=cgcolor('opposite',256L)
      isort = sort(airmass)
      items=string(airmass[isort],format='(f4.2)')
      legend, items,color=colors[isort],textcolor=colors[isort],pos=[xpos,ypos],box=0,charsize=1.0,/norm
      isort = sort(xfocal)
      items=string(xfocal[isort],format='(f0.1)')
      xpos = xpos + 0.075
      xyouts, xpos, ypos, '   xfocal',charsize=1.0,/norm, color=cgcolor('opposite',256L)
      legend, items,color=colors[isort],textcolor=colors[isort],pos=[xpos,ypos],box=0,charsize=1.0,/norm
      isort = sort(yfocal)
      items=string(yfocal[isort],format='(f0.1)')
      xpos = xpos + 0.075
      xyouts, xpos, ypos, '   yfocal',charsize=1.0,/norm, color=cgcolor('opposite',256L)
      legend, items,color=colors[isort],textcolor=colors[isort],pos=[xpos,ypos],box=0,charsize=1.0,/norm
   endif
   cgloadct,0

   !x.margin = xmargin
end

pro rm_diag_mratio_epoch, platelist, mjdlist, oneexp=oneexp,nprox=nprox,selfexclud=selfexclud $
       , xyfit=xyfit, nexp=nexp

; set /oneexp to plot the mratio for individual good std in each exposure

   if n_elements(platelist) eq 0 then begin
      ;platelist=[7338,7338,7338,7339,7339]
      ;mjdlist=[56660,56664,56669,56683,56686]
      target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
      fibermap = mrdfits(target_file,1,/silent)
      indd = where(fibermap[0].plate ne 0, nepoch)
      platelist = (fibermap[0].plate)[indd]
      mjdlist = (fibermap[0].mjd)[indd]
   endif

   ; default is to exclude the std itself in fluxing
   if n_elements(selfexclud) eq 0 then selfexclud = 1L

   nepoch = n_elements(platelist)
   calibdir=[replicate('/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/recalib/test20/',3), $
             replicate('/data3/quasar/yshen/spectro/bossredux/v5_6_0/7339/recalib/test20/',3)]
   topdir=[replicate('/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/',3), $
           replicate('/data3/quasar/yshen/spectro/bossredux/v5_6_0/7339/',3)]
   tag = string(platelist,format='(i4.4)')+'-'+string(mjdlist,format='(i5.5)')

   figfile = getenv('IDLRM_DIR')+'/misc/diag_mratio.ps'
   if keyword_set(oneexp) then figfile = getenv('IDLRM_DIR')+'/misc/diag_mratio_1exp.ps'
   if keyword_set(xyfit) then figfile = repstr(figfile,'.ps', '_xyfit.ps')
   begplot,name=figfile,/color
   !p.multi = [0,2,3]
   if keyword_set(oneexp) then !p.multi = [0,1,2]

   for iepoch = 0L, nepoch - 1L do begin

      rm_expinfo,platelist[iepoch],mjdlist[iepoch],expinfo=expinfo
      ind=where(expinfo.coadded[0,*] eq 1, ncoadd)
      if keyword_set(nexp) then ncoadd=nexp ; only plot the first nexp exposures
     
      calibdir1=calibdir[iepoch] & topdir1=topdir[iepoch]

      ; This plots the dispersion in mratio for all exps
      if not keyword_set(oneexp) then begin
         ; this is blue ccd1 [coadded exposures]
         spname=reform((expinfo[ind].name)[0,*])
         calibname=repstr(spname,'Frame','Fluxcalib') + '.gz'
         title = tag[iepoch] + ', Ncoadd='+string(ncoadd,format='(i0)')
         rm_checkmratio, calibname, calibdir=calibdir1, topdir=topdir1,title=title,$
           nprox=nprox,selfexclud=selfexclud,xyfit=xyfit

         ; this is red ccd1 [coadded exposures]
         spname=reform((expinfo[ind].name)[2,*])
         calibname=repstr(spname,'Frame','Fluxcalib') + '.gz'
         title = tag[iepoch] + ', Ncoadd='+string(ncoadd,format='(i0)')
         rm_checkmratio, calibname, calibdir=calibdir1, topdir=topdir1,title=title,$
           nprox=nprox,selfexclud=selfexclud,xyfit=xyfit

         ; this is blue ccd2 [coadded exposures]
         spname=reform((expinfo[ind].name)[1,*])
         calibname=repstr(spname,'Frame','Fluxcalib') + '.gz'
         title = tag[iepoch] + ', Ncoadd='+string(ncoadd,format='(i0)')
         rm_checkmratio, calibname, calibdir=calibdir1, topdir=topdir1,title=title,$
           nprox=nprox,selfexclud=selfexclud,xyfit=xyfit

         ; this is red ccd2 [coadded exposures]
         spname=reform((expinfo[ind].name)[3,*])
         calibname=repstr(spname,'Frame','Fluxcalib') + '.gz'
         title = tag[iepoch] + ', Ncoadd='+string(ncoadd,format='(i0)')
         rm_checkmratio, calibname, calibdir=calibdir1, topdir=topdir1,title=title,$
           nprox=nprox,selfexclud=selfexclud,xyfit=xyfit


      endif else begin
         
         ; this for each exposure
         for iexp=0L, ncoadd - 1L do begin            
            ; blue ccd1
            spname=reform((expinfo[ind].name)[0,iexp])
            calibname=repstr(spname,'Frame','Fluxcalib') + '.gz'
            title = tag[iepoch] + ', ' + strmid(spname,8,11)
            rm_checkmratio_1exp,calibname[0],calibdir=calibdir1, topdir=topdir1,$
              title=title,nprox=nprox,selfexclud=selfexclud,xyfit=xyfit
            ; red ccd1
            spname=reform((expinfo[ind].name)[2,iexp])
            calibname=repstr(spname,'Frame','Fluxcalib') + '.gz'
            title = tag[iepoch] + ', ' + strmid(spname,8,11)
            rm_checkmratio_1exp,calibname[0],calibdir=calibdir1, topdir=topdir1,$
              title=title, /noanno,nprox=nprox,selfexclud=selfexclud,xyfit=xyfit

            ; blue ccd2
            spname=reform((expinfo[ind].name)[1,iexp])
            calibname=repstr(spname,'Frame','Fluxcalib') + '.gz'
            title = tag[iepoch] + ', ' + strmid(spname,8,11)
            rm_checkmratio_1exp,calibname[0],calibdir=calibdir1, topdir=topdir1,$
              title=title,nprox=nprox,selfexclud=selfexclud,xyfit=xyfit
            ; red ccd2
            spname=reform((expinfo[ind].name)[3,iexp])
            calibname=repstr(spname,'Frame','Fluxcalib') + '.gz'
            title = tag[iepoch] + ', ' + strmid(spname,8,11)
            rm_checkmratio_1exp,calibname[0],calibdir=calibdir1, topdir=topdir1,$
              title=title, /noanno,nprox=nprox,selfexclud=selfexclud,xyfit=xyfit


         endfor

      endelse
   endfor

   !p.multi = 0
   endplot

   if keyword_Set(oneexp) then begin
      cmd = 'gzip ' + figfile
      spawn, cmd
   endif

end
