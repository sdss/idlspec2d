;+
; NAME:
;   dr1list_best
;
; PURPOSE:
;   Determine which plates *should* have been in DR1 that match the tiles
;   of those actually chosen.
;
; CALLING SEQUENCE:
;   dr1list_best
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Outputs are written to the file 'jill.out'.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   platelist
;   struct_print
;
; REVISION HISTORY:
;   22-Aug-2002  Written by David Schlegel, Princeton (not checked in then)
;-
;------------------------------------------------------------------------------
pro dr1list_best

   ; Read all plates and sort by TILE number
   platelist, plist=plist
   plist = plist[sort(plist.tile)]

   ; Find the list of tiles that should be in DR1
   indx = where(strmatch(plist.public,'*DR1*'))
   tiles = plist[indx].tile
   tiles = tiles[sort(tiles)]
   tiles = tiles[uniq(tiles)]
   ntile = n_elements(tiles)

   ; Make the logic of which plates correspond to these tile numbers
   nplate = n_elements(plist)
   qkeep = bytarr(nplate)
   for i=0, nplate-1 do $
    if (total(plist[i].tile EQ tiles) GT 0) then qkeep[i] = 1B

   ; For each tile in DR1, pick the best plate+MJD for that tile.
   qbest = bytarr(nplate)
   for itile=0, ntile-1 do begin
      jj = where(plist.tile EQ tiles[itile] $
       AND strmatch(plist.platequality,'good*'), nj)
      if (nj GT 0) then begin
         junk = max(plist[jj].platesn2, jbest)
         qbest[jj[jbest]] = 1B
      endif else begin
         jj = where(plist.tile EQ tiles[itile] $
          AND strmatch(plist.platequality,'marginal*'), nj)
         if (nj GT 0) then begin
            junk = max(plist[jj].platesn2, jbest)
            qbest[jj[jbest]] = 1B
         endif else begin
            print, 'No good/marginal plate for tile = ', tiles[itile]
         endelse
      endelse
   endfor

   pout = replicate(create_struct( $
    'plate', 0L, $
    'tile' , 0L, $
    'mjd'  , 0L, $
    'progname', '', $
    'platequality', '', $
    'platesn2', 0.0, $
    'public', '', $
    'qbest',  0B ), nplate)
   copy_struct, plist, pout
   pout.qbest = qbest
   struct_print, pout, filename='jill.out'

   return
end
;------------------------------------------------------------------------------
