; make the LC table for the first-lag paper

pro make_firstlag_table

file='/data3/quasar/yshen/work/lags/prepspec/ACBFJ/output/lc_all'
readcol, file, format='l,x,x,d,d,d,d,d,l',rmid, mjd, f_conti, e_conti, f_line, e_line, mask
mjd=mjd+50000.

outfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/first_lags/rv1/lc_table.txt'
openw, lun, outfile, /get_lun
fmt='(i3.3, " ", f9.3, " ", 4(e12.6, " "), i0)'

rmid_list=[101,191,229,267,272,320,457,589,645,694,767,769,775,789,840]

for i=0L, 14 do begin
  ind=where(rmid eq rmid_list[i], nnn)
  for j=0, nnn - 1 do $
    printf, lun, rmid[ind[j]], mjd[ind[j]], f_conti[ind[j]], e_conti[ind[j]], $
      f_line[ind[j]], e_line[ind[j]], mask[ind[j]], format=fmt
endfor

close, lun
free_lun, lun

end
