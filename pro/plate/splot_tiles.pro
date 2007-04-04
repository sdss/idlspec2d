;+
; NAME:
;   splot_tiles
; PURPOSE:
;   plot circular tiles from a given file
; CALLING SEQUENCE:
;   splot_tiles, filename
; INPUTS:
;   filename - name of Yanny .par file
; REVISION HISTORY:
;   26-Oct-2006  Written by MRB, NYU
;-
;------------------------------------------------------------------------------
pro splot_tiles, filename

radius=1.49
num=100L

tiles=yanny_readone(filename)

soplot, tiles.racen, tiles.deccen, psym=2, color='red'
for i=0L, n_elements(tiles)-1L do begin
    cosd=cos(!DPI/180.*tiles[i].deccen)
    th=findgen(num)/float(num-1L)*!DPI*2.
    rr=tiles[i].racen+radius*cos(th)/cosd
    dd=tiles[i].deccen+radius*sin(th)
    soplot, rr, dd, th=2, color='red'
endfor

end
