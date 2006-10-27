pro test_design_plate

filename='plPlugMapT-0583.par.orig'
hamin=15
hamax=15
airtemp=4

plugs=yanny_readone(filename, hdr=hdr)
for i=0L, n_elements(hdr)-1L do begin
    words=strsplit(hdr[i], /extr)
    if(strupcase(strtrim(string(words[0]),2)) eq 'RACEN') then $
      racen=double(words[1])
    if(strupcase(strtrim(string(words[0]),2)) eq 'DECCEN') then $
      deccen=double(words[1])
    if(strupcase(strtrim(string(words[0]),2)) eq 'TEMP') then $
      airtemp=double(words[1])
    if(strupcase(strtrim(string(words[0]),2)) eq 'HAMIN') then $
      hamin=double(words[1])
    if(strupcase(strtrim(string(words[0]),2)) eq 'HAMAX') then $
      hamax=double(words[1])
endfor
lst=racen+(0.5*(hamin+hamax))

ii=where(plugs.ra ne 0., nii)
plugs=plugs[ii]

simple_plate, plugs, racen=racen, deccen=deccen, airtemp=airtemp

;design_plate, plugs, racen=racen, deccen=deccen, airtemp=airtemp, $
;  nminsky=32, nstd=8

end
