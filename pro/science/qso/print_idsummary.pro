
pro print_idsummary, s, filename=filename

  nin = n_elements(s)

  if keyword_set(filename) then begin 
    openw, ilun, filename, /get_lun
    printf, ilun, 'INDEX Plate-MJD-FIB    RA(2000)      DEC    z_qso  <z_abs>  ew2796 ew2803 ew2852 ew2600 ew2382'
  endif else $
    print, 'INDEX Plate-MJD-FIB    RA(2000)      DEC    z_qso  <z_abs>  ew2796 ew2803 ew2852 ew2600 ew2382'

  for i=0,nin-1 do begin

     if keyword_set(filename) then $
     printf, ilun, s[i].index , s[i].plate, s[i].mjd, s[i].fiberid, s[i].ra, $
               s[i].dec, s[i].zqso, s[i].zabs, s[i].ew2796, s[i].ew2803, $
               s[i].ew2852, s[i].ew2600, s[i].ew2382, $
            format = '(i5.5, i5.4, "-",i5,"-",i3.3, 2d11.6, 2f8.5, 5f7.3)' $
     else $
     print, s[i].index , s[i].plate, s[i].mjd, s[i].fiberid, s[i].ra, s[i].dec,$
               s[i].zqso, s[i].zabs, s[i].ew2796, s[i].ew2803, s[i].ew2852, $
               s[i].ew2600, s[i].ew2382, $
            format = '(i5.5, i5.4, "-",i5,"-",i3.3, 2d11.6, 2f8.5, 5f7.3)'
   endfor

   if keyword_set(filename) then free_lun, ilun

return
end
     

