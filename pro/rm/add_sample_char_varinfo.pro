; add PrepSpec variability info to the sample char table

pro add_sample_char_varinfo


file='/data3/yshen/work/sdssrm_sample_char/sample_char_refined.fits'
tt=mrdfits(file,1)

; first create a struct that contains prepspec columns

prefix = ['RMS_ML_', 'RMS_ML_', 'SNR_RMS_ML_', 'RMS_ML_FRAC_', 'RMS_ML_FRAC_', 'N_RMS_GOOD_', 'SNR2_']
suffix = ['', '_ERR', '', '', '_ERR', '', '']
ncol = n_elements(prefix)
tag = ['C1700','C3000','C5100','HA','HB','HEII4687','MGII', 'CIII', 'CIV', 'LYA']
prepspec_name = ['c1700', 'c3000', 'c5100', 'ha_t', 'hb_t', 'he2_4686_t', 'mg2_t', 'c3_t', 'c4_t', 'lya_t'] + '_stats.dat'
ntag = n_elements(tag)
; now make the struct
str = create_struct('null', 0.)
for i=0, ntag - 1 do begin
  for j=0, ncol - 1 do begin
     name = prefix[j] + tag[i] + suffix[j]
     str = struct_addtags(str, create_struct(name, -1.D)  )
  endfor
endfor
remove_tags, str, 'null', new_str
alltags = tag_names(tt)
noldtag = n_elements(alltags)
str = new_str

nobj=n_elements(tt)
str = replicate(str, nobj)
tt = struct_addtags(tt, str)

outdir = '/data3/yshen/ftp/sdssrm/collab/prepspec/2014b/'
for i=0, nobj - 1 do begin
  rmtag = 'rm' + string(i,format='(i3.3)')
  subdir = outdir + rmtag + '/'

  for j=0, ntag - 1 do begin
    file = subdir + rmtag + '_' + prepspec_name[j]
    if file_test(file) eq 1 then begin
       data = 0
       readcol, file, format='d,a', skipline=28L, data,note, /silent
       ; rearrange order of data to match the tag order
       if n_elements(data) eq 7 then begin
          data1 = [data[1], data[2], data[3], data[4], data[5], data[0], data[6]]
          for kk=0,6 do tt[i].(noldtag + j*7 + kk) = data1[kk]
       endif
    endif
  endfor
  splog, 'Finished: ', i+1
endfor

outfile = '/data3/yshen/work/sdssrm_sample_char/sample_char_refined_w_var.fits'
mwrfits, tt, outfile, /create

end
