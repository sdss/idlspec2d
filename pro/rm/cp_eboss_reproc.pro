; all plates taken during 2015- (eBOSS) had an extraction problem with earlier version
; of eBOSS pipeline (the BOSS v5_7_1 pipeline works fine but was not run on most eBOSS plates)
; the new v5_10_10 reduction seems to fix the problem, which I reprocessed with custom 
; flux calibration and sky subtraction
; USE this routine to copy the reprocessed spPlate files to the ftp

pro cp_eboss_reproc, plate, mjd

fromdir = '/data3/yshen/spectro/bossredux/v5_10_10/'
todir = '/data3/yshen/ftp/sdssrm/collab/bossredux/v5_10_10/wh_skysub/'

nnn = n_elements(plate)

for i=0, nnn - 1 do begin
  pstr = string(plate[i],format='(i4.4)')
  mstr = string(mjd[i], format='(i5.5)')

  file = 'spPlate-' + pstr + '-' + mstr + '.fits'
  cmd = 'cp -p ' + fromdir + pstr + '/wh_skysub/' + file + ' ' + todir + pstr + '/'

  spawn, cmd
endfor


end
