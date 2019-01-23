; add columns to the sample characterization table

pro add_sample_char_host

file = '/data3/yshen/work/sdssrm_sample_char/sample_char_final.fits'
tt = mrdfits(file,1)
nnn = n_elements(tt)
column = {F_H_5100:0., sigma:0.D, sigma_err:-1.D, sigma_err_warning:0L}
tt = struct_addtags(tt, replicate(column, nnn))

file2 = '/data3/yshen/ftp/anon/public/sdssrm/paper_data/Shen_2015_ApJ_msigma/table1_apj.fits'
tt2 = mrdfits(file2,1)

nobj=n_elements(tt2)

for i=0, nobj-1 do begin
  ind = where(tt.rmid eq tt2[i].rmid)
  tt[ind].f_H_5100 = tt2[i].f_H_5100
  tt[ind].sigma = tt2[i].sigma
  tt[ind].sigma_err = tt2[i].sigma_err
  tt[ind].sigma_err_warning = tt2[i].sigma_err_warning
endfor

mwrfits, tt, file, /create

end

pro add_sample_char_col


file = '/data3/yshen/work/sdssrm_sample_char/sample_char_refined.fits'
tt = mrdfits(file,1)
nnn = n_elements(tt)
column = {Mi:0.D, Mi_z2:0.D}

tt = struct_addtags(tt, replicate(column, nnn))

mi = reform((tt.PSFmag)[3,*])
zz = reform(tt.zsys)
tt.mi_z2 = get_abs_mag(mi, zz)
tt.mi = tt.mi_z2 + 0.596


mwrfits, tt, file, /create

end
