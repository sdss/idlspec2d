; check the pt correction from PrepSpec, plot the
; difference in stddev in corrected and raw spec LCs
; as a function of emission line S/N measured from prepspec

pro check_pt, prepspecdir=prepspecdir, outfile=outfile, synfile=synfile

  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  epoch_info=mrdfits(file,2)
  if ~keyword_set(synfile) then file='/data3/yshen/ftp/sdssrm/collab/photo_lc/spec_syn.fits' $
    else file=synfile
  spec_syn=mrdfits(file,1)
  mjd_all=epoch_info[0:31].mean_mjd

  nnn=849 ; the RM targets
  result={rmid:0L, z:0.D,rms_raw:0.d, rms_pt:0.d, mad_raw:0.D, mad_pt:0.D, $
     flag_pt:0L, maxsn_nlr:-1.d, wave_maxsn:0.D}
  result=replicate(result, nnn)
  result.rmid=lindgen(nnn)
  result.z=target[0:nnn-1].zfinal

  if ~keyword_set(prepspecdir) then prepspecdir='/data3/yshen/ftp/sdssrm/collab/prepspec/ACBFJ/'
  for i=0, nnn-1 do begin
    ; find the difference in r band LC stddev
    rmag_raw=(spec_syn[i].gri_mag_ori)[1,*]
    indd=where(rmag_raw gt 0)
    rmag_raw=rmag_raw[indd]
    rmag_pt=(spec_syn[i].gri_mag)[1,*]
    indd=where(rmag_pt gt 0)
    rmag_pt=rmag_pt[indd]
    result[i].rms_raw=stddev(rmag_raw)
    result[i].mad_raw=median( abs(rmag_raw - median(rmag_raw) )  )
    result[i].rms_pt=stddev(rmag_pt)
    result[i].mad_pt=median( abs(rmag_pt - median(rmag_pt) ) )
 
    ; Is PrepSpec used for the object?
    pt=rm_get_pt(i,err=pt_err,mjd_all=mjd_all,topdir=prepspecdir)
    if (where(pt_err gt 1d-5))[0] ne -1 then result[i].flag_pt=1L

    ; find the NLR S/N measured by PrepSpec
    rmtag='rm'+string(i,format='(i3.3)')
    nlrfile=prepspecdir+rmtag+'/'+rmtag+'_vnlr.dat' 
    if file_test(nlrfile) eq 1 then begin
      readcol,nlrfile,format='d,d,d',wave_line,Jint,Jerr,/silent
      ; keep the avg and err from MC trials
      indd=where(Jerr gt -1d-5)
      NLR_sn=Jint[indd]/Jerr[indd]
      result[i].maxsn_nlr=max(nlr_sn, ind_max)
      result[i].wave_maxsn=wave_line[indd[ind_max]]
    endif

  endfor

  if ~keyword_set(outfile) then outfile='/data3/yshen/work/lags/prepspec/ACBFJ/rmag_scat.fits'
  mwrfits, result, outfile, /create

end

; plot the rms in r-band (synthetic) from pt corrected and uncorrected LCs
pro plot_rms_rband, infile=infile, outdir=outdir

  if ~keyword_set(infile) then file='/data3/yshen/work/lags/prepspec/ACBFJ/rmag_scat.fits' $
    else file=infile
  rms=mrdfits(file,1)
  ind=where(rms.flag_pt gt 0,nnn)
  indd=where(rms.flag_pt gt 0 and rms.z lt 0.7855)  

  if ~keyword_set(outdir) then outdir='/data3/yshen/work/lags/prepspec/ACBFJ/'

  figfile=outdir+'rmag_scat.ps'
  begplot, name=figfile,/landscape
  diff=rms[ind].rms_pt - rms[ind].rms_raw
  nlrsn=rms[ind].maxsn_nlr
  plot, nlrsn, diff, psym=2,/xlog, yrange=[-0.5, 1],xrange=[3,1d3],/xsty, $
     xtitle='Maximum Narrow Line S/N (coadd)', ytitle=textoidl('RMS_{pt}-RMS_{orig} [mag]'), $
     title='r-band (synthetic) LC stddev'
  oplot, rms[indd].maxsn_nlr,rms[indd].rms_pt - rms[indd].rms_raw,psym=2,color=cgcolor('red')
  oplot, [3,1d3], [0,0],line=2
  xyouts, 0.7, 0.9, textoidl('N_{obj}=')+string(nnn,format='(i0)'),/norm
  xyouts, 0.7, 0.82, textoidl('first 100 obj'), /norm, color=cgcolor('red')
  endplot
  cgfixps, figfile

  figfile=outdir+'rmag_mad.ps'
  begplot, name=figfile,/landscape
  diff=rms[ind].mad_pt - rms[ind].mad_raw
  nlrsn=rms[ind].maxsn_nlr
  plot, nlrsn, diff, psym=2,/xlog, yrange=[-0.3, 0.6],xrange=[3,1d3],/xsty, /ysty, $
     xtitle='Maximum Narrow Line S/N (coadd)', ytitle=textoidl('MAD_{pt}-MAD_{orig} [mag]'), $
     title='r-band (synthetic) LC MAD'
  oplot, rms[indd].maxsn_nlr,rms[indd].rms_pt - rms[indd].rms_raw,psym=2,color=cgcolor('red')
  oplot, [3,1d3], [0,0], line=2
  xyouts, 0.7, 0.9, textoidl('N_{obj}=')+string(nnn,format='(i0)'),/norm
  xyouts, 0.7, 0.82, textoidl('first 100 obj'), /norm, color=cgcolor('red')
  endplot
  cgfixps, figfile

end
