; output an ascii table of the observed peak wave for a list of lines
; from my qsofit for Keith to optimize the PrepSpec line-fitting windows

pro output_lineprop_keith

  file = '/data3/yshen/work/sdssrm_sample_char/sample_char_final.fits'
  target = mrdfits(file,1)
  tags = tag_names(target)

  line=['SII6718','HALPHA','HBETA','HEII4687','OIII5007','OII3728', $
        'NEV3426', 'MGII', 'CIII_ALL', 'NIII1750', 'CIV', 'HEII1640', $
        'SIIV_OIV', 'OI1304', 'LYA', 'NV1240']
  header = [['# RMID', 'ZPIP', 'ZSYS', 'ZSYS_ERR'], line] 
  zpip = target.zpip & zsys = target.zsys & zsys_err = target.zsys_err
  fmt = '(i3.3, " ", f7.5, " ", f7.5, " ", f7.5, 16(" ", f9.3, " ", f9.3) )'

  outfile = '/data3/yshen/work/sdssrm_sample_char/peak_wave.txt'
  openw, lun, outfile, /get_lun
  printf, lun, header
  nnn = n_elements(target)
  nline = n_elements(line)
  for i=0L, nnn - 1 do begin
      for j=0L, nline - 1 do begin
          ind1 = where(tags eq line[j]) & ind2 = where(tags eq line[j]+'_ERR')
          if j eq 0 then arr = [ (target[i].(ind1))[0], (target[i].(ind2))[0] > 0 ] * (1. + target[i].zpip) $
          else arr = [arr, [ (target[i].(ind1))[0], (target[i].(ind2))[0] > 0 ] * (1. + target[i].zpip)] 
      endfor        
       
      printf, lun, format=fmt, target[i].rmid, target[i].zpip, $
        target[i].zsys, target[i].zsys_err, arr
  endfor
  close, lun
  free_lun, lun

  ; now output FWHM
  line=['HALPHA_BR','HBETA_BR','HEII4687_BR','OIII5007', $
        'MGII_BR', 'CIII_ALL', 'CIV', 'HEII1640_BR', $
        'SIIV_OIV', 'LYA']
  header = [['# RMID', 'ZPIP', 'ZSYS', 'ZSYS_ERR'], line]
  fmt = '(i3.3, " ", f7.5, " ", f7.5, " ", f7.5, 10(" ", f9.3, " ", f9.3) )'

  outfile = '/data3/yshen/work/sdssrm_sample_char/bl_fwhm.txt'
  openw, lun, outfile, /get_lun
  printf, lun, header
  nnn = n_elements(target)
  nline = n_elements(line)
  for i=0L, nnn - 1 do begin
      for j=0L, nline - 1 do begin
          ind1 = where(tags eq line[j]) & ind2 = where(tags eq line[j]+'_ERR')
          if j eq 0 then arr = [ (target[i].(ind1))[1], (target[i].(ind2))[1] > 0 ] $
          else arr = [arr, [ (target[i].(ind1))[1], (target[i].(ind2))[1] > 0 ] ]
      endfor

      printf, lun, format=fmt, target[i].rmid, target[i].zpip, $
        target[i].zsys, target[i].zsys_err, arr
  endfor

  close, lun
  free_lun, lun


end

