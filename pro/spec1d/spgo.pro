PRO spgo

; Read template
  listfile = filepath('brg.08.gk.01', $
                      root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
  
  readcol, listfile[0], wave, template, /silent

  reduce_plate, 303, wave, template, result


  return
END
