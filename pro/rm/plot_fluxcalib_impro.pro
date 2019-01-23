; plot a figure to show the improvement of the custom flux calibration over pipeline

pro plot_fluxcalib_impro

  ; find all std
  std=rm_findstd(synflux=synflux, calibflux=calibflux, plot_diag=0)
  std=rm_findstd(synflux=synflux1, calibflux=calibflux1, plot_diag=0,calibdir='recalib/')  ; 'recalib/'

  ; get the synflux/calibflux ratio
  ratio = synflux/calibflux  ; [nepoch, nstd, nband]
  logratio = alog10(ratio) > (-99.)
  ratio1 = synflux1/calibflux1  ; [nepoch, nstd, nband]
  logratio1 = alog10(ratio1) > (-99.)
 
  figfile=getenv('IDLRM_DIR')+'/misc/fluxcalib_comp.eps'
  filtstr = ['u', 'g', 'r', 'i', 'z']+' band'
  begplot, name=figfile,/color, /cmyk
  !p.multi = [0,1,3]
  bsize=0.01 & charsize=3
  yrange=[0,500.]
  xtitle='log (synflux/calibflux)'
  pos=[0.12, 0.75, 0.95, 0.98]
  for iband=1L, 3L do begin
     
     plot, [0],[0], xrange=[-0.2,0.2] ,yrange=yrange,xtitle=xtitle,ytitle=textoidl('N_{STD}'), $
        charsize=charsize, /nodata, pos=pos, xticklen=0.04
     plothist, logratio[*,*,iband], bin=bsize, /over, xhist,yhist
     plothist, logratio1[*,*,iband], bin=bsize,/over,color=cgcolor('red'), xhist1,yhist1
     ; fit simple Gaussian to the dist
     fit=mpfitfun('gauss1', xhist,yhist,yhist,weights=1.D, [0., 0.1, 100.])
     fit1=mpfitfun('gauss1', xhist1,yhist1,yhist1,weights=1.D, [0., 0.1, 100.])
     xyouts, 0.13, 400, filtstr[iband]
     xyouts, 0.1, 300, textoidl('mean  \sigma')
     xyouts, 0.1, 250, string(fit[0],format='(f0.2)')+'  ' +string(fit[1],format='(f0.3)')
     xyouts, 0.1, 200, string(fit1[0],format='(f0.2)')+'  ' +string(fit1[1],format='(f0.3)'),$
      color=cgcolor('red')

     pos[1]=pos[1]-0.33
     pos[3]=pos[3]-0.33
     if iband eq 1 then begin
       xyouts, 0.18, 0.94, 'BOSS pipeline', /norm
       xyouts, 0.18, 0.91, 'custom', /norm, color=cgcolor('red')
     endif
  endfor
  !p.multi = 0
  endplot

end
