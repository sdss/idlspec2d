
pro id_summarize, id_list, spec, summary=s, filename=filename

  mg1 = where(fix(id_list.restwave) EQ 2796,nin)

  if nin EQ 0 then return

  template = { index: 0L, plate: 0L, mjd : 0L, fiberid : 0L, $
               ra : 0.0d, dec : 0.0d, imag : 0.0, zqso : 0.0, zabs : 0.0, $
               ew2796: 0.0, ew2803: 0.0, ew2852: 0.0, ew2600 : 0.0, $
               ew2382 : 0.0 , rank: 0.0 }

  s = replicate(template,nin)
  index = id_list[mg1].index
  s.index   = index
  s.ra      = spec[index].zans.plug_ra
  s.dec     = spec[index].zans.plug_dec
  s.plate   = spec[index].zans.plate
  s.mjd     = spec[index].zans.mjd
  s.fiberid = spec[index].zans.fiberid
  s.zqso    = spec[index].zans.z
;
;  Guess  i' band magnitude
;
  s.imag    = spec[index].plug.mag[3] - 0.2

;
;  Now fill in absorption information and print out
;
  tag_wave = strmid(tag_names(template),2)

  for i=0,nin-1 do begin

     lines = where(id_list.index EQ id_list[mg1[i]].index AND $
                   abs(id_list[mg1[i]].zabs - id_list.zabs) LT 0.001)

     zabs = 0.0
     weight = 0.0

     for j=9,13 do begin
       exist = where(tag_wave[j] EQ fix(id_list[lines].restwave))
       if exist[0] NE -1 then begin
          s[i].(j) = id_list[lines[exist]].ew_pix
          zabs = zabs + id_list[lines[exist]].zabs
          weight = weight + 1.0
       endif
     endfor

     s[i].zabs = zabs/weight

   endfor

   print_idsummary, s, filename=filename

return
end
     

