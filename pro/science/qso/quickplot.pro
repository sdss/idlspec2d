pro quickplot, tt, index, list, z=z

    if NOT keyword_set(z) then pz = 1.0 else pz = 1.0 + z

    flux = (1.0 - tt[index].flux)*tt[index].cont 

    splot, 10^tt[index].loglam/pz, flux, ps=10
    soplot, 10^tt[index].loglam/pz, tt[index].cont, color='blue', thick=3

    bad = where(tt[index].mask EQ 0)
    if bad[0] NE -1 then $
      soplot, 10^tt[index].loglam[bad]/pz, flux[bad], color='red', ps=4

    if size(list,/tname) NE 'STRUCT' then return
    lines = where(list.index EQ index)

    if lines[0] NE -1 then begin

      pix = long(list[lines].pixel)
      soplot, list[lines].wave/pz, tt[index].cont[pix]*1.1+1.0, ps=1, $
         symsize=3.0, color='green'
      sxyouts, list[lines].wave/pz, tt[index].cont[pix]*1.1+1.0, $
         list[lines].label, orient=90., charsize=1.8

    endif
    return
end

