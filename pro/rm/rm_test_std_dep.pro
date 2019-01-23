;+
; NAME:
;   rm_check_std_model
;   rm_test_std_dep
;
; PURPOSE:
;   Various tests upon the flux standards
;   Test dependence of std fluxcalib vector on various properties

pro rm_check_std_model, std_indx=std_indx, stdinfo=stdinfo

   platelist=[7338,7338,7338,7339,7339]
   mjdlist=[56660,56664,56669,56683,56686]

   ; default is to return all std
   if n_elements(std_indx) eq 0 then std_indx=lindgen(70)

   ; find the fiberid for all std on all epochs
   std_fiber=rm_findstd(platelist,mjdlist)
   std_fiber=std_fiber[*,std_indx]

   nepoch=n_elements(platelist) & nstd=n_elements(std_indx)
   kmodel=lonarr(nstd,nepoch)
   qgood=lonarr(nstd,nepoch)
   stdinfo=replicate({plate:0L,mjd:0L,fiberid:lonarr(nstd),std_indx:std_indx $
     ,kmodel:lonarr(nstd),qgood:lonarr(nstd)},nepoch)
   stdinfo.plate=platelist & stdinfo.mjd=mjdlist

   for i=0L, nepoch-1 do stdinfo[i].fiberid=std_fiber[i,*]

   ; now loop over each epoch to find the kunz model
   for i=0L, nepoch - 1L do begin

      rm_expinfo,platelist[i],mjdlist[i],expinfo=expinfo,/diet
      ind=where(expinfo.coadded[0,*] eq 1)
      expinfo=expinfo[ind]

      for j=0L, nstd - 1L do begin
         thisfiber = std_fiber[i,j]

         if thisfiber le 500 then $
          calibfile=(expinfo[0].name)[0] else $
           calibfile=(expinfo[0].name)[1]
         calibfile=repstr(calibfile,'spFrame','spFluxcalib')

         calibfile=getenv('BOSS_SPECTRO_REDUX') + '/' + $
            getenv('RUN2D') + '/' + string(platelist[i],format='(i4.4)')+'/' + $
            calibfile + '.gz'
         kindx=mrdfits(calibfile,2,/silent)
         indx = where(kindx.fiberid eq thisfiber)

         kmodel[j,i]=kindx[indx].imodel
         qgood[j,i]=kindx[indx].qgood
      endfor

   endfor
   stdinfo.kmodel = kmodel
   stdinfo.qgood = qgood
end

; Check the effect of fluxcorr
pro rm_check_fluxcorr,plate,mjd,fiberid,calibdir=calibdir,plotkey=plotkey,yrange=yrange,nsmooth=nsmooth

   if not keyword_set(nsmooth) then nsmooth = 1L

   !p.multi=[0,1,2]

   if not keyword_set(yrange) then yrange=[0,50]

   if not keyword_set(calibdir) then calibdir = $
     getenv('BOSS_SPECTRO_REDUX') + '/' + $
      getenv('RUN2D') + '/' + string(plate,format='(i4.4)')+'/'    

   ; Get exposure info
   rm_expinfo,plate,mjd,expinfo=expinfo,/diet
   ind=where(expinfo.coadded[0,*] eq 1,nexp)
   expinfo=expinfo[ind]

   colors = fsc_color(['opposite','red','green','cyan','blue','magenta' $
       ,'gold','dark green','Olive'])

   ; Locate the exposure for the fiber
   nfiber=n_elements(fiberid)
   for ifib=0,nfiber-1 do begin
      thisfiber = fiberid[ifib]
      
      if thisfiber le 500 then begin
        explist = (expinfo.name)[0,*]
        thisindx = thisfiber - 1L
      endif else begin
        explist = (expinfo.name)[1,*]
        thisindx = thisfiber - 501L 
      endelse

      ; loop over the coadded exposures
      ; first do the blue ccd
      plot,[0],[0],xrange=[3600,6000],yrange=yrange,/xsty,/ysty,/nodata
      for iexp=0L,nexp - 1L do begin

         spfile = explist[iexp] + '.gz'
         ; read in the raw spectrum
         ; the wave dispesion has been corrected
         rm_spframe_read,spfile,thisindx,loglam=loglam,objflux=objflux, $
            objivar=objivar,plate=plate
  
         ; apply flux calibration
         calibfile=calibdir+repstr(spfile,'spFrame','spFluxcalib')
         calibfac = (mrdfits(calibfile,0,/silent))[*,thisindx]
         minval = 0.05 * mean(calibfac)
         divideflat, objflux, invvar=objivar, calibfac, minval=minval
         if plotkey[0] eq 1 then $
          oplot,10.D^loglam,smooth(objflux,nsmooth),color=colors[iexp],line=0

         if n_elements(wave) eq 0 then begin
            wave = loglam
            flux = objflux
            ivar = objivar
         endif else begin
            wave = [ [wave], [loglam] ]
            flux = [ [flux], [objflux]]
            ivar = [ [ivar], [objivar]]
         endelse

         ; apply flux-correction factor
         corrfile=calibdir+repstr(spfile,'spFrame','spFluxcorr')
         aterm = (mrdfits(corrfile, 0, corrhdr, /silent))[*,thisindx] > 0.7 < 1.3
         bterm = (mrdfits(corrfile, 1,/silent))[*,thisindx]
         invertcorr = 1. / aterm
         minval = 0.05 / mean(aterm)
         nrownative=(size(objflux,/dimens))[0]
         divideflat, objflux, invvar=objivar, invertcorr[0:nrownative-1,*], minval=minval
         objflux = objflux + bterm
         if plotkey[1] eq 1 then $
          oplot,10.D^loglam,smooth(objflux,nsmooth),color=colors[iexp],line=0

         if plotkey[1] eq 2 then $
          oplot,10.D^loglam,aterm[0:nrownative-1,*],color=colors[iexp],line=0

         if n_elements(fluxcorr) eq 0 then begin
            fluxcorr = objflux
            ivarcorr = objivar
         endif else begin
            fluxcorr = [ [fluxcorr], [objflux]]
            ivarcorr = [ [ivarcorr], [objivar]]
         endelse

         ;message,'stop and diagnose'

      endfor
      items = strmid(explist,8,11)
      legend,items, box=0,pos=[0.7,0.9],/norm,color=colors,textcolor=colors

      ; Next do the red ccd
      if thisfiber le 500 then begin
        explist = (expinfo.name)[2,*]
        thisindx = thisfiber - 1L
      endif else begin
        explist = (expinfo.name)[3,*]
        thisindx = thisfiber - 501L
      endelse
      plot,[0],[0],xrange=[6000,10000.],yrange=yrange,/xsty,/ysty,/nodata
      for iexp=0L,nexp - 1L do begin

         spfile = explist[iexp] + '.gz'
         ; read in the raw spectrum
         ; the wave dispesion has been corrected
         rm_spframe_read,spfile,thisindx,loglam=loglam,objflux=objflux, $
            objivar=objivar,plate=plate
  
         ; apply flux calibration
         calibfile=calibdir+repstr(spfile,'spFrame','spFluxcalib')
         calibfac = (mrdfits(calibfile,0,/silent))[*,thisindx]
         minval = 0.05 * mean(calibfac)
         divideflat, objflux, invvar=objivar, calibfac, minval=minval
         if plotkey[0] eq 1 then $
          oplot,10.D^loglam,smooth(objflux,nsmooth),color=colors[iexp],line=0
         wave = [ [wave], [loglam[0:4111]] ]
         flux = [ [flux], [objflux[0:4111]]]
         ivar = [ [ivar], [objivar[0:4111]]]

         ; apply flux-correction factor
         corrfile=calibdir+repstr(spfile,'spFrame','spFluxcorr')
         aterm = (mrdfits(corrfile, 0, corrhdr, /silent))[*,thisindx] > 0.7 < 1.3
         bterm = (mrdfits(corrfile, 1,/silent))[*,thisindx]
         invertcorr = 1. / aterm
         minval = 0.05 / mean(aterm)
         nrownative=(size(objflux,/dimens))[0]
         divideflat, objflux, invvar=objivar, invertcorr[0:nrownative-1,*], minval=minval
         objflux = objflux + bterm
         if plotkey[1] eq 1 then $
          oplot,10.D^loglam,smooth(objflux,nsmooth),color=colors[iexp],line=0
         if plotkey[1] eq 2 then $
          oplot,10.D^loglam,aterm[0:nrownative-1,*],color=colors[iexp],line=0

         fluxcorr = [ [fluxcorr], [objflux[0:4111]]]
         ivarcorr = [ [ivarcorr], [objivar[0:4111]]]

      endfor

   endfor
  
   !p.multi = 0
 
   if plotkey[0] eq 3 and plotkey[1] eq 3 then begin
 
       finalwave = 3.6 + findgen(4001)*1d-4

       combine1fiber, wave, flux, ivar, newloglam=finalwave, newflux=finalflux
       combine1fiber, wave, fluxcorr, ivarcorr, newloglam=finalwave, newflux=finalfluxcorr

       plot, 10.^finalwave, finalflux, xrange=[4000, 1d4]
       oplot, 10.^finalwave, finalfluxcorr,color=fsc_color('red')

       ; message, 'stop and diagnose'

   endif

end

pro rm_test_std_dep, plate, mjd, calibdir=calibdir

   if n_elements(plate) eq 0 then plate = 7338
   if n_elements(mjd) eq 0 then mjd = 56660

   topdir = getenv('BOSS_SPECTRO_REDUX') + '/' + $
      getenv('RUN2D') + '/7338/'
   if not keyword_set(calibdir) then calibdir=topdir+'recalib/test4/'   

   ; get exp info
   rm_expinfo, plate, mjd, expinfo=expinfo
   ind = where(expinfo.coadded[0,*] eq 1, nexp)
   expinfo = expinfo[ind]

   plotfile = 'stdmratio_dep'  + '.ps'
   set_plot, 'ps'
   dfpsplot, djs_filepath(plotfile, root_dir=getenv('IDLRM_DIR')+'/misc'),$
          /color ;, /landscape


   for i=0L, nexp - 1L do begin

      ; for each camera
      for j=0L, 3L do begin

         spframefile = expinfo[i].name[j] + '.gz'
         spframe_read,topdir+spframefile,plugmap=plugmap,hdr=hdr1

         calibfile = repstr(spframefile,'spFrame','spFluxcalib')
         kindx=mrdfits(calibdir+calibfile,2)
         struct_out = mrdfits(calibdir+calibfile,4)

         ; keep only the good std
         ind = where(kindx.qgood eq 1, ngood)
         kindx = kindx[ind]

         ; get airmass info
         airmass = expinfo[i].airmass
         if strmatch(calibfile,'*calib-b1*') or strmatch(calibfile,'*calib-r1*') then $
           ind_std = kindx.fiberid - 1L else $
           ind_std = kindx.fiberid - 1L - 500L 
         airmass = airmass[ind_std]

         ; get focal plane info
         focal_dist = sqrt(plugmap.xfocal^2 + plugmap.yfocal^2)       
         focal_dist = focal_dist[ind_std]

         minairmass=min(airmass, max=maxairmass)

         thisloglam=struct_out.thisloglam
         thismratio=struct_out.thismratio
         thismrativar=struct_out.thismrativar
         thisflatarr=struct_out.thisflatarr
   
         if strmatch(calibfile,'*calib-b*') then begin
           xrange=[3600,6300]
           yrange=[0,80]
         endif else begin
           xrange=[5800, 1d4]
           yrange=[0,100]
         endelse

         cgloadct, 13, /silent
         !p.multi = [0,1,2]

         ; plot dependence on airmass
         min=minairmass & max=maxairmass
         plot,[0],[0],xrange=xrange,/xsty,yrange=yrange,/nodata,title=calibfile,xtitle='Wavelength' $
            , ytitle='Flux Calibration Vector'

         for ii=0,ngood-1L do begin
            color_use = floor( (airmass[ii]-min)/(max-min)*255L ) < 255L
            yarr = smooth(THISmratio[*,0,ii]*thisflatarr[*,0,ii],10)
            oplot,10.^thisloglam[*,0,ii],yarr,color=color_use
         endfor

         pos1 = [0.14, 0.92, 0.92, 0.94]
         cgcolorbar,/norm, pos = pos1, color=cgcolor('firebrick',255), range=[min, max] $
           , minor=5,ncolors=255
         xyouts, 0.14, 0.87,'airmass',/norm

         ; plot dependence on focal distance to the origin [0,0]
         min=min(focal_dist, max=max)
         plot,[0],[0],xrange=xrange,/xsty,yrange=yrange,/nodata,title=calibfile,xtitle='Wavelength' $
            , ytitle='Flux Calibration Vector'

         for ii=0,ngood-1L do begin
            color_use = floor( (focal_dist[ii]-min)/(max-min)*255L ) < 255L
            yarr = smooth(THISmratio[*,0,ii]*thisflatarr[*,0,ii],10)
            oplot,10.^thisloglam[*,0,ii],yarr,color=color_use
         endfor

         pos1 = [0.14, 0.42, 0.92, 0.44]
         cgcolorbar,/norm, pos = pos1, color=cgcolor('firebrick',255), range=[min, max] $
           , minor=5,ncolors=255
         xyouts, 0.14, 0.37,'focal distance',/norm

         cgloadct, 0
         !p.multi = 0
      endfor

   endfor

   dfpsclose

end


