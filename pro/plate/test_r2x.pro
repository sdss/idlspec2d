pro test_r2x, filename

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
help,/st,plugs

radec_to_xyfocal, plugs.ra, plugs.dec, xfocal, yfocal, racen=racen, $
  deccen=deccen, airtemp=airtemp, lst=lst

radec_to_xyfocal, plugs.ra, plugs.dec, xfocal2, yfocal2, racen=racen, $
  deccen=deccen, airtemp=airtemp, lst=lst, /norefrac

splot, xfocal-plugs.xfocal, yfocal-plugs.yfocal, psym=4
splot, xfocal2-plugs.xfocal, yfocal2-plugs.yfocal, psym=4
rfocal=sqrt(plugs.xfocal^2+plugs.yfocal^2)
splot, rfocal, rfocal/plugs.yfocal, psym=4
splot, rfocal, xfocal/plugs.xfocal, psym=4
splot, xfocal, xfocal-plugs.xfocal, psym=4
splot, yfocal, yfocal-plugs.yfocal, psym=4
splot_vec, xfocal, yfocal, xfocal-plugs.xfocal, yfocal-plugs.yfocal, $
  scale=250.
save, filename='blah.sav'

end
