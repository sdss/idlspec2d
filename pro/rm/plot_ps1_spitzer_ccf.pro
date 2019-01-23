; make CCF plots of the PS1xSpitzer data

pro plot_ps1_spitzer_ccf, channel=channel, sp_dir=sp_dir, fig_suffix=fig_suffix, ccf_outdir=ccf_outdir

file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
target=mrdfits(file,1)

file = '/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/spitzer/spitzer_c11-12_raw_lc.fits'
ps1_lc = mrdfits(file,1)

spherematch, target.ra, target.dec, ps1_lc.ra_1, ps1_lc.dec_1, 0.1/3600.D, match1, match2, dist1
rmid = match1
nnn = n_elements(rmid)

; re-arrange the orders
ps1_lc = ps1_lc[match2]
str = replicate({RMID:-1}, nnn)
str.rmid = rmid
ps1_lc = struct_addtags(ps1_lc, str)
outfile = '/data3/quasar/yshen/work/spitzer/sdss-rm_ps1_lc.fits'
if file_test(outfile) ne 1 then mwrfits, ps1_lc, outfile, /create

if ~keyword_set(sp_dir) then sp_dir = '/data3/quasar/yshen/work/spitzer/spitzer-sdss-rm2/'  ; 2 means 2-year spitzer
if ~keyword_set(channel) then channel=2
if ~keyword_set(fig_suffix) then fig_suffix = '_ccf2'
figfile='/data3/quasar/yshen/work/spitzer/ps1_spitzer_ch' + string(channel, format='(i0)') + fig_suffix + '.ps'
begplot, name=figfile

if ~keyword_set(ccf_outdir) then ccf_outdir = '/data3/quasar/yshen/work/spitzer/ccf2/'

for i=0, nnn-1 do begin
   !P.multi = [0,1,3]
   charsize=2.5 & ticklen=0.02

   sp_file = sp_dir + 'sdss-rm-' + string(rmid[i],format='(i3.3)')+'-ch' + string(channel, format='(i0)') + '-161202.txt'
   if file_test(sp_file) ne 1 then sp_file = sp_dir + 'sdss-rm-' + string(rmid[i],format='(i3.3)')+'-ch' + string(channel, format='(i0)') + '-161205.txt'   

   ccf_outfile = ccf_outdir + 'sdss-rm-' + string(rmid[i],format='(i3.3)') + '_ccf_ch' + string(channel, format='(i0)') 
   fmt = '(f6.1, " ", f6.2, " ", i3)'

   ;print, sp_file
   if file_test(sp_file) eq 1 then begin
       openw, lun, ccf_outfile, /get_lun
       printf, lun, '# tau  ccf  ndata'

       mjd_ps1 = ps1_lc[i].lc_mjd_g
       mag_ps1 = ps1_lc[i].lc_mag_g
       err_ps1 = ps1_lc[i].lc_err_g
       ind0 = where(err_ps1 gt 0)
       mjd_ps1 = mjd_ps1[ind0]
       mag_ps1 = mag_ps1[ind0]
       err_ps1 = err_ps1[ind0]
       ind0 = sort(mjd_ps1)  ; sort the LC (which was not originally sorted)
       mjd_ps1 = mjd_ps1[ind0]
       mag_ps1 = mag_ps1[ind0]
       err_ps1 = err_ps1[ind0]

       ; use the PS1 data as continuum to be interpolated
       t_conti = mjd_ps1
       f_conti = 10.D^(-0.4*mag_ps1)*3631d6 ; microJy
       e_conti = err_ps1*alog(10.D)/2.5*f_conti 
       ; reject apparent outliers
       tmp = f_conti[sort(f_conti)]
       nep_tmp = n_elements(tmp)
       tmp = tmp[0.05*nep_tmp:0.95*nep_tmp]
       ind_good = where(abs(f_conti - median(f_conti)) le 5.*stddev(tmp) )
       t_conti = t_conti[ind_good]
       f_conti = f_conti[ind_good]
       e_conti = e_conti[ind_good]

       readcol, sp_file, format='d,d,d', mjd_sp, mag_sp, err_sp
       ind0 = sort(mjd_sp)
       mjd_sp = mjd_sp[ind0] & mag_sp = mag_sp[ind0] & err_sp = err_sp[ind0]
       t_line = mjd_sp
       f_line = 10.D^(-0.4*mag_sp)*3631d6 ; microJy
       e_line = err_sp*alog(10.D)/2.5*f_line

       ; create a time array
       tau_min = min(mjd_sp) - max(mjd_ps1)
       tau_max = max(mjd_sp) - min(mjd_ps1)
       dt = 10.
       ccf_tau = tau_min + findgen((tau_max - tau_min)/dt)*dt
       ccf_xrange = [min(ccf_tau), max(ccf_tau)]

       ccf0 = xcorr_interp(t_line,f_line,t_conti,f_conti,tau=ccf_tau,ndata=ndata0 $
            , nogap=0, cent_tau=cent_tau0, peak_tau=peak_tau0 $
            , peak_sig=peak_sig)
       print, rmid[i], peak_sig

       title = 'RMID'+string(rmid[i],format='(i3.3)') + ', PS1-g, Spitzer Ch' + string(channel, format='(i0)')
       pos = [0.11,0.65,0.95,0.95]
       plot,[0],[0],xrange=ccf_xrange,yrange=yrange,xtickname=replicate(' ',10), $
       ytitle='CCF/DCF', $
       pos=pos,charsize=charsize,xticklen=ticklen, title=title
       oplot, ccf_tau, ccf0 ; ,color=colors
       oplot, [peak_tau0, peak_tau0],[0,1],linestyle=2
       xyouts, 2000, 0.88, textoidl('peak_{sig}=')+string(peak_sig,format='(f0.3)')
 
       dt4 = dt*2. ; this is the DCF grid
       tau_grid = tau_min + findgen((tau_max - tau_min)/dt4)*dt4
       rm_dcf,t_conti,f_conti,t_line,f_line,err1=median(e_conti),err2=median(e_line),$
              DCF=DCF,eDCF=eDCF,tau_grid=tau_grid
       ind1 = where(eDCF gt 0)
       tmp = where(dcf gt 1)
       if tmp[0] ne -1 then print, tau_grid[tmp], dcf[tmp]
       if ~keyword_set(noplot) then begin
         if ind1[0] ne -1 then begin
           oploterror, tau_grid[ind1],dcf[ind1],edcf[ind1],psym=symcat(6),color=cgcolor('red'),$
            errcolor=cgcolor('red')
         endif
       endif

       nnn=n_elements(ccf_tau)
       for ii=0L, nnn-1 do printf, lun, format=fmt, ccf_tau[ii], ccf0[ii], ndata0[ii]
       close, lun
       free_lun, lun

       pos = [0.11,0.34,0.95,0.64]
       plot, ccf_tau, ndata0, xrange=ccf_xrange, ytitle = textoidl('N_{data}'), pos=pos, xtitle='Time Delay (Days)', yrange=[0,50],charsize=charsize

       pos = [0.11, 0.07, 0.95, 0.27]
       xrange = [min(t_conti),max(t_line)]
       ploterror, t_conti, f_conti, e_conti, psym=5, xtitle='MJD', ytitle='LC (microJy)', xrange=xrange,pos=pos, XTICKFORMAT='(i0)',charsize=charsize, symsize=0.5
       norm1 = median(f_line)/median(f_conti)
       oploterror, t_line, f_line/norm1, e_line/norm1, psym=5,color=cgcolor('red'),errcolor=cgcolor('red'),symsize=0.5
      

       !P.multi = 0
   endif

endfor
endplot

end
