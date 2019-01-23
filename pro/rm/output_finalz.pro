; Make a final list of redshifts for the 849 RM targets

pro output_finalz

  coadd_MJD=56697
  rm_readspec,0000,mjd=coadd_mjd,zans=zans

  ; manually fix some redshifts based on visual inspection
  VI_file= getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
       + '/0000/qsofit/visual_z'
  readcol, vi_file, format='l,d', fiber_vi, z_vi
  zans[fiber_vi-1].z = z_vi

  outfile=getenv('IDLRM_DIR')+'/etc/zfinal.list'
  openw, lun, outfile,/get_lun
  printf, lun, '# Plate=0000 (coadd)'
  printf, lun, '# MJD='+string(coadd_mjd,format='(i5.5)')
  printf, lun, '# Fiber  zfinal'  

  fmt = '(i4.4, " ", f6.4)'
  for i=0L, 848L do begin
    
    printf, lun, format=fmt, i+1, zans[i].z

  endfor

 
  close, lun
  free_lun, lun

end
