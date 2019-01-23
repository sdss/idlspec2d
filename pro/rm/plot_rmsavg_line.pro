; plot the MEDIAN avg/rms and rms/avg profiles

pro plot_rmsavg_line, line=line

file = '/data3/yshen/work/sdssrm_sample_char/sample_char_final.fits'
tt = mrdfits(file,1)

minsnr2 = 10.
if line eq 'Hb' then ind = where(tt.rms_ml_frac_hb gt 0 and tt.snr2_hb gt minsnr2)
if line eq 'MgII' then ind=where(tt.rms_ml_frac_mgii gt 0 and tt.snr2_mgii gt minsnr2)
if line eq 'CIV' then ind=where(tt.rms_ml_frac_civ gt 0 and tt.snr2_civ gt minsnr2)

figfile='/data3/yshen/ftp/sdssrm/collab/prepspec/2014b/' + line + '_AVGRMS.ps'
begplot, name=figfile,/color

make_rmsavg_line, ind, line=line, result=result
wave = (result[0].wave/result[0].lam0 - 1.)*3d5

max = max(result.rms_avg, dim=1)
indd = where(max gt 0 and max lt 5, nobj)

plot, wave, median(result[indd].avg,dim=2), pos=[0.12, 0.70,0.95, 0.95], title=line+', N='+string(nobj,format='(i0)'), /norm, $
   xtickname = replicate('', 10L)
oplot, [0,0],[0,10],line=1
xyouts, 0.65, 0.92, 'median AVG', /norm, charsize=charsize
plot, wave, median(result[indd].rms,dim=2), pos=[0.12, 0.40,0.95, 0.65], /norm, /noerase,xtickname = replicate('', 10L)
oplot, [0,0],[0,10],line=1
xyouts, 0.65, 0.62, 'median RMS', /norm, charsize=charsize
plot, wave, median(result[indd].rms_avg,dim=2), pos=[0.12, 0.1,0.95, 0.35], /norm, /noerase, xtitle='Velocity [km/s]'
oplot, [0,0],[0,10],line=1
xyouts, 0.65, 0.32, 'median RMS/AVG', /norm, charsize=charsize


endplot

end
