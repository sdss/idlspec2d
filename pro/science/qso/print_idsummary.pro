
pro print_idsummary, s, filename=filename

  nin = n_elements(s)

  
  if keyword_set(filename) then begin 
    openw, ilun, filename, /get_lun
    printf, ilun, 'INDEX Plate-MJD-FIB    RA(2000)      DEC   iband  z_qso  <z_abs>  RANK  2796 2803 2852 2600 2382'
  endif else $
    print,        'INDEX Plate-MJD-FIB    RA(2000)      DEC   iband  z_qso  <z_abs>  RANK  2796 2803 2852 2600 2382'

  for i=0,nin-1 do begin

     if keyword_set(filename) then $
     printf, ilun, s[i].index , s[i].plate, s[i].mjd, s[i].fiberid, s[i].ra, $
               s[i].dec, s[i].imag, s[i].zqso, s[i].zabs, s[i].rank, $
               s[i].ew2796*1000, s[i].ew2803*1000, s[i].ew2852*1000, $
               s[i].ew2600*1000, s[i].ew2382*1000, $
          format = '(i5.5, i5.4, "-",i5,"-",i3.3, 2d11.6, f6.2, 2f8.5, f7.3, 5i5)' $
     else $
     print, s[i].index , s[i].plate, s[i].mjd, s[i].fiberid, s[i].ra, s[i].dec,$
               s[i].imag, s[i].zqso, s[i].zabs, s[i].rank, $
               s[i].ew2796*1000, s[i].ew2803*1000, $
               s[i].ew2852*1000, s[i].ew2600*1000, s[i].ew2382*1000, $
            format = '(i5.5, i5.4, "-",i5,"-",i3.3, 2d11.6, f6.2, 2f8.5, f7.3, 5i5)'
   endfor

   if keyword_set(filename) then free_lun, ilun

return
end
     

