; test the SN dependence of spectra measurements
; using artifically downgraded spectra

;----------------------------
pro plot_specfit_sntest, figfile=figfile, name=name, flag_norm=flag_norm, $
   civ3_errdir=civ3_errdir

; name='peak', 'fwhm' or 'ew'
if n_elements(name) eq 0 then name='peak'

; flag_norm=1 (normalized by measurement error); flag_norm=2 (normalized by original measurements)
if n_elements(flag_norm) eq 0 then flag_norm=1

if ~keyword_set(figfile) then begin
  if name eq 'peak' then begin
    indname=0  ; index in the best-fit line properties
  endif
  if name eq 'fwhm' then begin
    indname=1
  endif
  if name eq 'ew' then indname=3
endif
figfile='/data3/yshen/work/sdssrm_sample_char/figs/sntest_' + name + '.eps'

begplot,name=figfile,/color,/landscape
ang = string(197B)
cspeed=2.9979246d5

; this is the original fitting results
file='/data3/yshen/work/lineshifts/lineshift.fits'
result0=mrdfits(file,1,/silent)
; replace CIV measurements with the 3-gaussian fit
file_civ3 = '/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/civ_3gauss/qso_prop-0000-56837.fits'
tmp = mrdfits(file_civ3, 1, /silent)
result0.CIV = tmp.civ & result0.civ_err = tmp.civ_err
result0.civ_br = tmp.civ_br & result0.civ_br_err = tmp.civ_br_err
result0.HEII1640 = tmp.HEII1640 & result0.HEII1640_ERR = tmp.HEII1640_ERR
result0.heII1640_br = tmp.heII1640_br & result0.heII1640_br_err = tmp.heII1640_br_err

tags0=tag_names(result0)

if ~keyword_set(errdir) then errdir = '/data3/yshen/work/lineshifts/lineshift_err_scal_'

charsize=0.8
; if n_elements(hist) eq 0 then hist=1
histpeak=0.3 & histbin=0.1 & len=0.04 & thick=2

;line=['Hbeta_br', 'OIII5007', 'CaII3934', 'OII3728', 'NeV3426', 'MgII', 'CIII', 'HeII1640', 'CIV_br', 'SIIV_OIV']
;line=strupcase(line)
;linetag=textoidl(['H\beta_{br}', '[OIII]_a', 'CaII', '[OII]', '[NeV]', 'MgII', 'CIII]_a', 'HeII1640', 'CIV', 'SiIV'])
;wave=[4862.68, 5008.24, 3934.78, 3728.48, 3426.84, 2798.75, 1908.73,1640.42, 1549.06, (1396.76 + 1402.06)*0.5]

; my EW measurements do not include the 3 locally fitted lines 'CaII3934', 'OII3728', 'NeV3426'
line=['Hbeta_br', 'OIII5007', 'MgII', 'CIII', 'HeII1640', 'CIV_br', 'SIIV_OIV']
line=strupcase(line)
linetag=textoidl(['H\beta_{br}', '[OIII]_a', 'MgII', 'CIII]_a', 'HeII1640', 'CIV', 'SiIV'])
wave=[4862.68, 5008.24, 2798.75, 1908.73,1640.42, 1549.06, (1396.76 + 1402.06)*0.5]

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
      ; replace CIV measurements with the 3gauss fit
      file_civ3 = '/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/civ_3gauss/err_scal_'+ string(err_scal_arr[j],format='(f0.1)') + '/qso_prop-0000-56837.fits'
      tmp = mrdfits(file_civ3,1,/silent)
      ;result1.CIV = tmp.civ & result1.civ_err = tmp.civ_err
      result1.civ_br = tmp.civ & result1.civ_br_err = tmp.civ_err
      result1.HEII1640 = tmp.HEII1640 & result1.HEII1640_ERR = tmp.HEII1640_ERR
      ;result1.heII1640_br = tmp.heII1640_br & result1.heII1640_br_err = tmp.heII1640_br_err
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
              ind_use = where(abs(dif) lt 1d3)
              dif = dif[ind_use]
              if j eq 0 then plothist, dif, bin=bsize, pos=pos, noerase=noerase, xticklen=len, yticklen=len, $
                 title=title+'  '+textoidl('N_{tot}=')+string(ntot,format='(i0)'), $
                 color=colors[j], xrange=xrange, charsize=charsize,xhist,yhist,thick=thick,/xsty $
              else plothist, dif, bin=bsize,/over,color=colors[j], thick=thick, xhist, yhist
              oplot, [median(dif), median(dif)], [0,ntot], color=colors[j], thick=thick, line=2

              area=int_tabulated(xhist,double(yhist),/double)
              ;oplot, xgrid, gauss1(xgrid, [0, 1, area]), color=cgcolor('dark gray')
               
              
              if j eq 0 then begin
                 xyouts, pos[0]+0.01, pos[3]-0.03, textoidl('<SNR_{med}>'),charsize=charsize,/norm
                 if name ne 'ew' then $
                  xyouts, pos[2]-0.09, pos[3]-0.03, textoidl('<\sigma_{mea}> [kms^{-1}]'),charsize=charsize,/norm $
                 else $
                  xyouts, pos[2]-0.09, pos[3]-0.03, textoidl('<\sigma_{mea}> [')+ang+']',charsize=charsize,/norm
              endif
              medsn = median(medsn_all[indd])/err_scal_arr[j]
              xyouts, pos[0]+0.01, pos[3]-0.02*j-0.05, string(medsn,format='(f0.1)'),/norm, color=colors[j],charsize=charsize
              if name eq 'peak' then begin
                 sig_mea_tot=median(err_tot)/wave[i]*cspeed
                 sig_mea0=median(err0[indd])/wave[i]*cspeed
                 sig_mea1=median(err1[indd])/wave[i]*cspeed
              endif
              if name eq 'fwhm' or name eq 'ew' then begin
                 sig_mea_tot=median(err_tot)
                 sig_mea0=median(err0[indd])
                 sig_mea1=median(err1[indd])
              endif

              print, title, ' ', string(medsn,format='(f0.1)'), median(dif), ' ntot=',ntot
              print, sig_mea_tot, sig_mea0, sig_mea1
              if name ne 'ew' then $
     xyouts, pos[2]-0.05, pos[3]-0.02*j-0.05, string(round(sig_mea1),format='(i0)'),/norm, color=colors[j],charsize=charsize $
              else xyouts, pos[2]-0.05, pos[3]-0.02*j-0.05, string(sig_mea1,format='(f0.1)'),/norm, color=colors[j],charsize=charsize

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
if name eq 'ew' then xname=textoidl('Normalized REW Difference [\sigma]')
xyouts, 0.5, 0.01, xname, /norm, align=0.5
xyouts, 0.02, 0.5, textoidl('N_{obj}'), /norm, align=0.5, orient=90

endplot
cgfixps, figfile


end


