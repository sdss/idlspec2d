; D. Finkbeiner
; 27 Jun 2000
; Create structure to contain results of mveldisp
FUNCTION mveldisp_struc, N, ntemplate

  IF keyword_set(ntemplate) THEN arr = fltarr(ntemplate) ELSE arr = 0.

  struc = { $
            plate               : 0,   $
            mjd                 : 0L,  $
            fiber               : 0,   $
            zchic               : 0.0, $
            class               : '',  $
            primtarget          : 0L,  $
            zbest               : 0.0, $
            zbest_err           : 0.0, $
            sigma2best          : 0.0, $
            sigma2best_err      : 0.0, $
            z                   : arr, $
            z_err               : arr, $
            zconf               : arr, $
            sigma2_cc           : arr, $
            sigma2_cc_err       : arr, $
            sigma2_quotient     : arr, $
            sigma2_quotient_err : arr, $
            sigma2_diff         : arr, $
            sigma2_diff_err     : arr, $
            run                 : 0L,  $
            rerun               : 0L,  $
            camcol              : 0L,  $
            field               : 0L,  $
            id                  : 0L   $
}
  
  IF keyword_set(N) THEN $
    arr = replicate(struc, N) $
  ELSE $
    arr = struc
  
  return, arr
END 
