; delete the mismatched-fiber ascii spectra in /data3/yshen/ftp/sdssrm/collab/bossredux/v5_7_1/ascii_spec/wh_skysub/

pro delete_misfiber_epoch

file = '/data3/yshen/ftp/sdssrm/collab/bossredux/fiber_mismatch_eboss/corr_fiber_eboss_all'
readcol,file,format='l,l,l,l',rmid, plate, mjd, fiber

nnn = n_elements(rmid)
dir = '/data3/yshen/ftp/sdssrm/collab/bossredux/v5_7_1/ascii_spec/wh_skysub/'
for i=0, nnn - 1 do begin

  tag = string(plate[i],format='(i4.4)')+'-'+string(mjd[i],format='(i5.5)')+ $
   '-' + string(fiber[i], format='(i4.4)')

  filename = dir + 'RMID_' + string(rmid[i],format='(i3.3)')+'/' + $
     '*' + tag + '.dat.gz'

  spawn, 'rm ' + filename

endfor

end
