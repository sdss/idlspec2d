; plot the peak detections against line SNR output by PrepSpec
; using some threshold for the statistical signifcance of the peak
; this is to test if there is a preference for positive lags

pro plot_lag_snr_test, min_peak_sig=min_peak_sig, figfile=figfile, ccf_file=ccf_file,plot_det=plot_det

if n_elements(plot_det) eq 0 then plot_det=1

if ~keyword_set(min_peak_sig) then min_peak_sig=0.999

file='/data3/quasar/yshen/work/lags/prepspec/ACBFJ/snrplot.dat'
readcol,file,format='l,x,x,x,d,a', rmid_snr, logsnr_line, line_snr
line_snr=strupcase(line_snr)

if not keyword_set(ccf_file) then ccf_file='/data3/quasar/yshen/work/lags/prepspec/ACBFJ/ccf_sym_range_rej_ep7'
readcol,ccf_file,format='l,a,x,d,d,x,d', rmid, line, lag, lag_low, peak_sig
line=strupcase(line)
; keep only hb and mg2 measurements
ind=where(line eq 'MG2' or line eq 'HB', nobj)
rmid=rmid[ind] & line=line[ind] & lag=lag[ind] & peak_sig=peak_sig[ind]
logsnr=dblarr(nobj)-99
for i=0L, nobj-1 do begin
  indd=where(rmid_snr eq rmid[i] and line_snr eq line[i],nnn)
  if nnn eq 1 then logsnr[i]=logsnr_line[indd]
endfor


if not keyword_set(figfile) then figfile='/data3/quasar/yshen/work/lags/prepspec/ACBFJ/peak_snr.eps'
;'/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/peak_snr.eps'
begplot, name=figfile, /color,/encap, /cmyk,ysize=5
pos=[0.12,0.13, 0.95,0.92]
plot, [0], [0], xrange=[-100,100], yrange=[0, 2.1], xtitle='CCF Peak [days]', /ysty, $
 ytitle='log SNR (line variability)', /nodata, pos=pos,title='based on spectroscopic LCs only'
oplot, [-100,100],[1.,1.],line=2
ind=where(line eq 'HB')
oplot, lag[ind], logsnr[ind],psym=1, color=cgcolor('dark gray')
ind=where(line eq 'HB' and peak_sig gt min_peak_sig and logsnr gt 0, nn)
if ind[0] ne -1 then oplot, lag[ind], logsnr[ind],psym=1, color=cgcolor('red')
; decide how many are positive and negative lags with r>0.999
ind1=where(line eq 'HB' and peak_sig gt min_peak_sig and lag ge 0 and logsnr gt 0,nn1)
ind2=where(line eq 'HB' and peak_sig gt min_peak_sig and lag lt 0 and logsnr gt 0,nn2)
print, 'Hbeta(tot/pos/neg lags):', nn, nn1, nn2

ind=where(line eq 'MG2')
oplot, lag[ind], logsnr[ind],psym=2, color=cgcolor('dark gray')
ind=where(line eq 'MG2' and peak_sig gt min_peak_sig and logsnr gt 0, nn)
if ind[0] ne -1 then oplot, lag[ind], logsnr[ind],psym=2, color=cgcolor('red')
; decide how many are positive and negative lags with r>0.999
ind1=where(line eq 'MG2' and peak_sig gt min_peak_sig and lag ge 0 and logsnr gt 0,nn1)
ind2=where(line eq 'MG2' and peak_sig gt min_peak_sig and lag lt 0 and logsnr gt 0,nn2)
print, 'MgII(tot/pos/neg lags):', nn, nn1, nn2

xyouts, 0.16, 0.86, 'All', color=cgcolor('dark gray'), /norm
xyouts, 0.16, 0.82, 'Peak sig>' + string(min_peak_sig,format='(f0.3)'), color=cgcolor('red'),/norm
;xyouts, 0.16, 0.78, 'Reported detections', color=cgcolor('green'),/norm
;xyouts, 0.55, 0.86, textoidl('\circ')+': Detections', color=cgcolor('black'),/norm
legend, pos=[0.5,0.9], ': Detections', box=0,psym=symcat('9'), /norm

xyouts, 0.75, 0.86, textoidl('+: H\beta'), /norm
xyouts, 0.75, 0.82, textoidl('*: MgII'), /norm

; now overplot the reported detections
if keyword_set(plot_det) then begin

; Here I am using my measurements of lags instead of Kate Grier's values for the 15
; detection reported in Table 1 of the first-lag paper; this is because Kate's values
; do not align with the points from the full 100 objects, which were used in all the red
; points 
file='/data3/quasar/yshen/work/lags/prepspec/ACBFJ/output/keep/tau_all'
readcol,file,format='l,a,d,d',rmid,line,zz,tau_rest
line=strupcase(line)
nobj=n_elements(rmid)
tau=tau_rest*(1. + zz)
logsnr=dblarr(nobj)-99
for i=0L, nobj-1 do begin
  indd=where(rmid_snr eq rmid[i] and line_snr eq line[i],nnn)
  if nnn eq 1 then logsnr[i]=logsnr_line[indd]
endfor

ind=where(line eq 'HB')
oplot, tau[ind], logsnr[ind],psym=1, color=cgcolor('red')
oplot, tau[ind], logsnr[ind],psym=symcat(9),symsize=1.4, color=cgcolor('red')
ind=where(line eq 'MG2')
oplot, tau[ind], logsnr[ind],psym=2, color=cgcolor('red')
oplot, tau[ind], logsnr[ind],psym=symcat(9), symsize=1.4, color=cgcolor('red')

endif


endplot

end

