PRO spgo

; Read template
  listfile = filepath('brg.08.gk.01', $
                      root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
  
  readcol, listfile[0], wave, template, /silent

  reduce_plate, 303, wave, template, result
  reduce_plate, 305, wave, template, result
  reduce_plate, 306, wave, template, result
  reduce_plate, 307, wave, template, result
  reduce_plate, 308, wave, template, result
  reduce_plate, 309, wave, template, result


  return
END


; Explore Eisenstein templates
PRO sp5

  listfile = filepath('template5gal.dat', $
                      root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
  
  readcol, listfile[0], wave, a, b, c, d, e

  ztemplate = [[a], [b], [c], [d], [e]]
  n_ztemp = (size(ztemplate))[2]
  ztemplate_wave = wave#(fltarr(n_ztemp))

  reduce_plate, 308, wave, ztemplate, result
stop

  return
END
