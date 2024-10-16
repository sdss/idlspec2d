; test the SN dependence of line center measurements
; using artifically downgraded spectra

pro test_lineshift_sn, err_scal_arr=err_scal_arr, range=range, emparfile=emparfile, outdir = outdir, $
        linename=linename
; this is to test the global fitting for most broad lines

;err_scal_arr=[2.,4.,6.,8.,10.]
;err_scal_arr=[6.,8.,10.]
if ~keyword_set(err_scal_arr) then err_scal_arr=[8.]

nerr=n_elements(err_scal_arr)

if ~keyword_set(outdir) then outdir = '/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/'
if ~keyword_set(emparfile) then emparfile=getenv('IDLRM_DIR')+'/etc/qsoline_more.par'

for i=0L, nerr - 1 do begin

   outdir1= outdir + 'err_scal_' $
     + string(err_scal_arr[i],format='(f0.1)') + '/'

   if ~keyword_set(linename) then begin
       linename=['SII6718', 'HALPHA_BR', 'OIII5007', 'OIII5007C', 'HBETA_BR', 'HEII4687_BR', $
                 'MGII', 'CIII', 'CIII_BR', 'NIII1750', 'CIV_BR', 'HEII1640', 'SIIV_OIV', 'OI1304', $
                 'LYA_BR', 'NV1240']
   endif

   rm_fitplate,0,56837L,/silent,emparfile=emparfile,outdir=outdir1,linename=linename, $
       err_scal=err_scal_arr[i], add_bhmass=0, range=range

endfor

end

;---------------------------
pro test_lineshift_sn_local
; this is to test the local fitting for CaII, [OII] and [NeV]

;err_scal_arr=[2.,4.,6.,8.,10.]
err_scal_arr=[4.,6.,8.,10.]

nerr=n_elements(err_scal_arr)

for i=0L, nerr - 1 do begin

   outdir='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/local_fit/err_scal_' $
     + string(err_scal_arr[i],format='(f0.1)') + '/'

   linename=['CAII3934', 'OII3728', 'NEV3426']
   emparfile=getenv('IDLRM_DIR')+'/etc/qsoline_few.par'

   rm_fitplate_local,0,56837L,/silent,emparfile=emparfile,outdir=outdir,linename=linename, $
       err_scal=err_scal_arr[i],add_bhmass=0

endfor

end
;----------------------------
pro plot_lineshift_sntest, figfile=figfile, name=name, flag_norm=flag_norm, origfits=origfits, $
   errdir=errdir

; name='peak', 'fwhm' or 'ew'
if n_elements(name) eq 0 then name='peak'

; flag_norm=1 (normalized by measurement error); flag_norm=2 (normalized by original measurements)
if n_elements(flag_norm) eq 0 then flag_norm=1

if ~keyword_set(figfile) then begin
  if name eq 'peak' then begin
    figfile='/data3/yshen/work/lineshifts/vel_shift_sntest.eps'
    indname=0  ; index in the best-fit line properties
  endif
  if name eq 'fwhm' then begin
    figfile='/data3/yshen/work/lineshifts/fwhm_sntest_sig_norm.eps'
    indname=1
  endif
endif

begplot,name=figfile,/color,/landscape

cspeed=2.9979246d5

; this is the original fitting results
if ~keyword_set(origfits) then file='/data3/yshen/work/lineshifts/lineshift.fits' $
  else file=origfits
result0=mrdfits(file,1,/silent)
tags0=tag_names(result0)

if ~keyword_set(errdir) then errdir = '/data3/yshen/work/lineshifts/lineshift_err_scal_'

charsize=0.8
; if n_elements(hist) eq 0 then hist=1
histpeak=0.3 & histbin=0.1 & len=0.04 & thick=2

line=['Hbeta_br', 'OIII5007', 'CaII3934', 'OII3728', 'NeV3426', 'MgII', 'CIII', 'HeII1640', 'CIV_br', 'SIIV_OIV']
line=strupcase(line)
linetag=textoidl(['H\beta_{br}', '[OIII]_a', 'CaII', '[OII]', '[NeV]', 'MgII', 'CIII]_a', 'HeII', 'CIV', 'SiIV'])
wave=[4862.68, 5008.24, 3934.78, 3728.48, 3426.84, 2798.75, 1908.73,1640.42, 1549.06, (1396.76 + 1402.06)*0.5]

err_scal_arr=[2.,4.,6.,8.,10.]
;medsn=27.2434/err_scal_arr

; get the medsn for all 849 objects
rm_readspec, 0, mjd=56837, flux=fluxall, invvar=ivar_all
medsn_all =  median(fluxall*sqrt(ivar_all), dim=1)
medsn_all=medsn_all[0:848]

nerr=n_elements(err_scal_arr)
; set up plot layout
npanel = n_elements(line)
plot_layout, npanel, xypos=xypos, omargin=[0.06, 0.01, 0.98, 0.96], pmargin=[0.05,0.1]
colors=cgcolor(['black', 'blue', 'cyan', 'magenta', 'red'])
if flag_norm eq 1 then begin
   xrange=[-5,5] & bsize=0.2
endif
if flag_norm eq 2 then begin
   xrange=[-1, 0.5] & bsize=0.05
endif
xgrid=-5. + findgen(101)*0.1
for i=0L, npanel - 1 do begin
   if i eq 0 then noerase=0 else noerase=1
   title=linetag[i]
   pos = xypos[*, i]
   for j=0, nerr-1 do begin
      file=errdir+ string(err_scal_arr[j],format='(f0.1)') + '.fits'
      result1=mrdfits(file,1,/silent)
      tags1=tag_names(result1)      

      ind0=where(tags0 eq line[i],nn0) & inderr_0=where(tags0 eq line[i]+'_ERR')
      ind1=where(tags1 eq line[i],nn1) & inderr_1=where(tags1 eq line[i]+'_ERR')
      if nn0 eq 1 and nn1 eq 1 then begin

           wav0=(result0.(ind0))[indname,*] & err0=(result0.(inderr_0))[indname,*]
           wav1=(result1.(ind1))[indname,*] & err1=(result1.(inderr_1))[indname,*]
           indd = where(err0 gt 0 and err1 gt 0, ntot)
           if indd[0] ne -1 then begin

              err_tot=sqrt( err0[indd]^2 + err1[indd]^2  )
              if flag_norm eq 1 then dif = (wav1[indd] - wav0[indd])/err_tot
              if flag_norm eq 2 then dif = (wav1[indd] - wav0[indd])/wav0[indd]
              if j eq 0 then plothist, dif, bin=bsize, pos=pos, noerase=noerase, xticklen=len, yticklen=len, $
                 title=title+'  '+textoidl('N_{tot}=')+string(ntot,format='(i0)'), $
                 color=colors[j], xrange=xrange, charsize=charsize,xhist,yhist,thick=thick,/xsty $
              else plothist, dif, bin=bsize,/over,color=colors[j], thick=thick, xhist, yhist
              oplot, [median(dif), median(dif)], [0,ntot], color=colors[j], thick=thick, line=2

              area=int_tabulated(xhist,double(yhist),/double)
              ;oplot, xgrid, gauss1(xgrid, [0, 1, area]), color=cgcolor('dark gray')
               
              
              if j eq 0 then begin
                 xyouts, pos[0]+0.01, pos[3]-0.03, textoidl('<SNR_{med}>'),charsize=charsize,/norm
                 xyouts, pos[2]-0.09, pos[3]-0.03, textoidl('<\sigma_{mea}> [kms^{-1}]'),charsize=charsize,/norm
              endif
              medsn = median(medsn_all[indd])/err_scal_arr[j]
              xyouts, pos[0]+0.01, pos[3]-0.02*j-0.05, string(medsn,format='(f0.1)'),/norm, color=colors[j],charsize=charsize
              if name eq 'peak' then begin
                 sig_mea_tot=median(err_tot)/wave[i]*cspeed
                 sig_mea0=median(err0[indd])/wave[i]*cspeed
                 sig_mea1=median(err1[indd])/wave[i]*cspeed
              endif
              if name eq 'fwhm' then begin
                 sig_mea_tot=median(err_tot)
                 sig_mea0=median(err0[indd])
                 sig_mea1=median(err1[indd])
              endif

              print, title, ' ', string(medsn,format='(f0.1)'), median(dif), ' ntot=',ntot
              print, sig_mea_tot, sig_mea0, sig_mea1
              xyouts, pos[2]-0.05, pos[3]-0.02*j-0.05, string(round(sig_mea1),format='(i0)'),/norm, color=colors[j],charsize=charsize

              if j eq 0 and flag_norm eq 1 then $
      oplot, xgrid, gauss1(xgrid, [0, 1, ntot*bsize]), color=cgcolor('dark gray'),thick=6 ; this should be almost identical to the dark gray lines
           endif
      endif
   endfor
    
endfor
if name eq 'peak' then xname=textoidl('Normalized Peak Velocity Difference [\sigma]')
if name eq 'fwhm' then begin
  if flag_norm eq 1 then xname=textoidl('Normalized FWHM Difference [\sigma]')
  if flag_norm eq 2 then xname=textoidl('\DeltaFWHM/FWHM')
endif
xyouts, 0.5, 0.01, xname, /norm, align=0.5
xyouts, 0.02, 0.5, textoidl('N_{obj}'), /norm, align=0.5, orient=90

endplot
cgfixps, figfile


end


