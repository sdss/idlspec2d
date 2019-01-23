; output ascii files from fits table

pro fits2ascii, result, fitsfile=fitsfile, column=column, fmt=fmt, outfile=outfile

if ~keyword_Set(outfile) then outfile='ascii.out'

if n_elements(result) eq 0 then $
  result=mrdfits(fitsfile, 1,column=column)
nnn=n_elements(result)

openw, lun, outfile, /get_lun

printf, lun, '# ', column
for i=0L, nnn-1 do begin

  printf, lun, result[i], format=fmt

endfor


close, lun
free_lun, lun

end
