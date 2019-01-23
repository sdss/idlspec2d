;+
; NAME:
;   rm_findstd
;
; PURPOSE:
;   find std stars of each epoch of the RM spectroscopy
;
; CALLING SEQUENCE:
;   rm_findstd, [platelist, mjdlist, calibdir='recalib/test20/',tag='xyfit']
;
; INPUTS:
;   platelist  - list of plates of the epochs [NEPOCH]
;   mjdlist    - list of mjd of the epochs [NEPOCH]
;
; OPTIONAL INPUTS:
;   
;  
;  
;
; OUTPUTS:
;   fiberid    - list of fiberid of each epoch [NEPOCH, NSTD]
; OPTIONAL OUTPUTS:
;
; COMMENTS:
; 
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;----------------------------------------------------------------------------

function rm_findstd, platelist, mjdlist, calibdir=calibdir, synflux=synflux, $
           calibflux=calibflux, plot_diag=plot_diag, expframe=expframe,tag=tag

    if not keyword_set(topdir) then topdir = getenv('BOSS_SPECTRO_REDUX') $
           + '/' + getenv('RUN2D') + '/'

    nepoch = n_elements(mjdlist)
    nplate = n_elements(platelist)
    if nplate eq 1 and nplate lt nepoch then $
        platelist = replicate(platelist, nepoch)

    ; read in the master file
    target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
    fibermap = mrdfits(target_file,1,/silent)

    ; get all the available epochs
    if nepoch eq 0 then begin
      indd = where(fibermap[0].plate ne 0, nepoch)
      platelist = (fibermap[0].plate)[indd]
      mjdlist = (fibermap[0].mjd)[indd]
    endif

    istd = where(strtrim(fibermap.sourcetype) eq 'STD', nstd)
    fiberid = lonarr(nepoch, nstd)

    for i=0L, nepoch - 1L do begin

       indd = where(fibermap[0].plate eq platelist[i] and $
           fibermap[0].mjd eq mjdlist[i])
       fiberid[i,*] = (fibermap.fiberid)[indd[0], istd]    
   
    endfor

    ; get the synthsized 5-band flux and calibration flux
    synflux = fltarr(nepoch, nstd, 5L)
    calibflux = fltarr(nepoch, nstd, 5L) ; check if calibflux is constant
    for i=0L, nepoch - 1L do begin

      if n_elements(calibdir) ne 0 then path = topdir+string(platelist[i],format='(i4.4)')+'/'+calibdir

      ; read the coadded spectrum
      if not keyword_set(expframe) then begin 
         readspec, platelist[i], (fiberid[i,*]), mjd=mjdlist[i], wave=wave, $
           flux=flux,invvar=invvar,andmask=mask,path=path,plugmap=plugmap,/silent
      endif else begin
         if n_elements(expframe) ne nepoch then message, '# of expframe must equal Nepoch'
         spframe_read, topdir + string(platelist[i],format='(i4.4)') $
              + '/spCFrame-b1-' + expframe[i] + '.fits' $
              , objflux=flux1, loglam=loglam1, plugmap=plugmap1
         ;print, 
         spframe_read, topdir + string(platelist[i],format='(i4.4)') $
              + '/spCFrame-b2-' + expframe[i] + '.fits' $
              , objflux=flux2, loglam=loglam2, plugmap=plugmap2
         flux = [[flux1], [flux2]] & loglam = [[loglam1], [loglam2]]
         wave = 10.^loglam
         ; keep only std
         flux = flux[*, fiberid[i,*]-1L] 
         wave = wave[*, fiberid[i,*]-1L]
         plugmap = ([plugmap1, plugmap2])[reform(fiberid[i,*]-1L)]
      endelse
         
      flambda2fnu = wave^2 / 2.99792e18
      fthru = filter_thru(flux * flambda2fnu, waveimg=wave,/toair)
      thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)
      synflux[i,*,*] = 10.^((22.5-thismag)/2.5)
      calibflux[i,*,*] = transpose(plugmap.CALIBFLUX)

    endfor

    ; plot diagnositic plots if required
    if keyword_set(plot_diag) then begin

       plotfile = 'rmDiagstd.ps'
       if (not keyword_set(tag)) or (not keyword_set(calibdir)) then tag = 'pipe'
       plotfile = repstr(plotfile,'.ps','_'+tag)+'.ps'
       ;set_plot, 'ps'
       ;dfpsplot, djs_filepath(plotfile, root_dir=getenv('IDLRM_DIR')+'/misc'), /color
       figname = djs_filepath(plotfile, root_dir=getenv('IDLRM_DIR')+'/misc')
       begplot, name=figname,/color
       xrange = [0, nepoch+1] & yrange = [0, 8]      
       csize=1.5
     
       platestr = string(platelist,format='(i0)' )
       mjdstr = string(mjdlist, format='(i0)')
       filtstr = 'filter ' + ['u', 'g', 'r', 'i', 'z']
       ;for iband=0L, 4L do begin
       ;   !p.multi = [0,2,5]
       ;   
       ;   for ipanel=0L, 9L do begin
       ;
       ;      yoffset = 0
 
       ;      istd = 0L + ipanel*7
       ;      title = filtstr[iband] + ', std=' $
       ;          + string(istd+1,format='(i0)')+'-'+string(istd+7,format='(i0)')
       ;      ratio = synflux[*,istd,iband]/calibflux[*,istd,iband]
       ;      plot, (1+indgen(nepoch)), ratio,title=title,$
       ;          xrange=xrange, yrange=yrange,/xsty, /ysty, psym=2, /nodata,charsize=csize
       ;      ind = where(ratio ge 0.95 and ratio le 1.05, ngood)
       ;      if ind[0] ne -1 then oplot, (1+indgen(nepoch))[ind], ratio[ind] $
       ;         , psym=2
 
       ;      ind = where( (ratio gt 1.05 or ratio lt 0.95) and ratio ne 0, nbad)
       ;      if ind[0] ne -1 then oplot, (1+indgen(nepoch))[ind], ratio[ind] $
       ;         , psym=2,color=fsc_color('red')
       ;      ind = where(ratio eq 0, nworst)
       ;       if ind[0] ne -1 then oplot, (1+indgen(nepoch))[ind], ratio[ind] $
       ;          , psym=5,color=fsc_color('blue'),symsize=1.5
       ;
       ;      for jj=0L, nepoch-1L do xyouts,1+jj,1.05 + yoffset, $
       ;          string(fiberid[jj,istd],format='(i0)')
       ;      oplot, xrange, replicate(1.05 + yoffset,2),line=1
       ;      oplot, xrange, replicate(0.95 + yoffset,2),line=1

       ;      for inew=1L, 6L do begin

       ;         istd = inew + ipanel*7
       ;         yoffset = yoffset + 1.

       ;         ratio = synflux[*,istd,iband]/calibflux[*,istd,iband]
       ;         ind = where(ratio ge 0.95 and ratio le 1.05, ngood)
       ;         if ind[0] ne -1 then oplot, (1+indgen(nepoch))[ind], ratio[ind] + yoffset, psym=2
       ;         ind = where( (ratio gt 1.05 or ratio lt 0.95) and ratio ne 0, nbad)
       ;         if ind[0] ne -1 then oplot, (1+indgen(nepoch))[ind], ratio[ind] + yoffset $
       ;           , psym=2,color=fsc_color('red')
       ;         ind = where(ratio eq 0, nworst)
       ;         if ind[0] ne -1 then oplot, (1+indgen(nepoch))[ind], ratio[ind] $
       ;           , psym=5,color=fsc_color('blue'),symsize=1.5

       ;         for jj=0L, nepoch-1L do xyouts,1+jj,1.05 + yoffset, $
       ;          string(fiberid[jj,istd],format='(i0)')

       ;         oplot, xrange, replicate(1.05 + yoffset,2),line=1
       ;         oplot, xrange, replicate(0.95 + yoffset,2),line=1
       ;      endfor
       ;   endfor

       ;   !p.multi = 0
         
      ; endfor

       ; NOW, make plots of the distribution of synflux/calibflux
       !p.multi = [0,1,3]
       for iband=1L, 3L do begin
          ; plot the distribution for each epoch
          colors = cgcolor(['black', 'red', 'green', 'magenta', 'cyan', 'blue', 'gold'])
          ncolor = n_elements(colors)
          title = filtstr[iband] + ', ' + tag ;'filter'+string(iband+1,format='(i0)')
          xrange = [-0.2,0.2] & yrange = [0,25]
          plot, [0], [0], xrange=xrange,yrange=yrange,xtitle=textoidl('log synflux/calibflux'), $
              title = title,/xsty, /ysty, charsize=csize
          xpos = xrange[0] + 0.05*(xrange[1]-xrange[0])
          ypos = yrange[1] - 0.1*(yrange[1]-yrange[0])
          dypos = 0.05*(yrange[1]-yrange[0])
          for iepoch=0L, nepoch - 1L do begin
             ratio = synflux[iepoch,*,iband]/calibflux[iepoch,*,iband]
             ; remove unplgged fibers (synflux=0)
             ind = where(ratio gt 0)
             ;if platelist[iepoch] eq 7338 and mjdlist[iepoch] eq 56660 then $
               ; ind = where(ratio gt 0 and fiberid[iepoch,*] ne 998)
 
             if ind[0] ne -1 then begin
               plothist, alog10(ratio[ind]), bin = 0.01, color=colors[iepoch mod ncolor],/over
                
               ;rms_logratio = stddev(alog10(ratio[ind]))
               ;instead of using the stddev, try to use 0.5*(per_84-per_16) 
               arr=alog10(ratio[ind])
               rms_logratio = 0.5*(quantile_2d(0.84,arr,dim=1) - quantile_2d(0.16,arr,dim=1))
               rms_std = stddev(alog10(ratio[ind]))

               xyouts, xpos+(iepoch/16)*0.6*(xrange[1]-xrange[0]), ypos, $
                       platestr[iepoch]+'-'+mjdstr[iepoch] + ':' + $
                       string(rms_logratio,format='(f0.4)')+'   ' + $
                       string(rms_std,format='(f0.4)'), charsize=1., $
                        color=colors[iepoch mod ncolor]
             endif
             ypos = ypos - dypos
             if iepoch eq 15 then ypos=yrange[1] - 0.1*(yrange[1]-yrange[0])
          endfor
       endfor
       !p.multi = 0

       ;dfpsclose
       endplot
    endif

    return, fiberid

end





