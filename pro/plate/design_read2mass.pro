; Read a list of 2MASS objects from an ASCII file as generated by Jill Knapp.
;
; The 2MASS magnitudes are roughly related to the Tycho B+R magnitudes
; as follows:
;   B = J + 4.51 * (J-K) + 0.22
;   R = J + 2.26 * (J-K) + 0.08
;   (B-R) = 2.25 * (J-K) + 0.14
; with a standard deviation of 0.65 mag in B, and 0.40 mag in R. This
; extrapolation would be much, much worse for late-type stars not in Tycho.
;------------------------------------------------------------------------------
function design_read2mass, filename

   readcol, filepath(filename, $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory=['pro','plate']), $
    ra, dec, bmag, rmag, jmag, hmag, kmag, flag, $
    format='(D,D,F,F,F,F,F,A)'

   bad_bmag = (bmag LT -9 OR bmag GT 90)
   bad_rmag = (rmag LT -9 OR rmag GT 90)

   ibad = where(bad_bmag, nbad)
   if (nbad GT 0) then begin
      bmag[ibad] = rmag[ibad] + 1.3 ; Approximate as B-R = 1.3
   endif

   ibad = where(bad_rmag, nbad)
   if (nbad GT 0) then begin
      rmag[ibad] = bmag[ibad] - 1.3 ; Approximate as B-R = 1.3
   endif

   ibad = where(bad_bmag AND bad_rmag, nbad)
   if (nbad GT 0) then begin
      bmag[ibad] = jmag[ibad] + 4.51 * (jmag[ibad] - kmag[ibad]) + 0.22
      rmag[ibad] = jmag[ibad] + 2.26 * (jmag[ibad] - kmag[ibad]) + 0.08
   endif

   result = design_starstruct(n_elements(ra))
   result.ra = ra
   result.dec = dec
;   result.mag = transpose( [[bmag], [rmag], [jmag], [hmag], [kmag]] )
   result.mag = transpose( [[bmag], [bmag], [rmag], [rmag], [rmag]] )
   result.objtype = 'SERENDIPITY_MANUAL'

   return, result
end
;------------------------------------------------------------------------------
