; D. Finkbeiner
; 27 Jun 2000
; Create structure to contain results of veldisp
FUNCTION veldisp_struc, N


  struc = { $
            plate              : 0,   $
            mjd                : 0L,  $
            fiber              : 0,   $
            zchic              : 0.0, $
            class              : '',  $
            primtarget         : 0L,  $
            z                  : 0.0, $
            z_err              : 0.0, $
            sigma_cc           : 0.0, $
            sigma_cc_err       : 0.0, $
            sigma_quotient     : 0.0, $
            sigma_quotient_err : 0.0, $
            sigma_diff         : 0.0, $
            sigma_diff_err     : 0.0, $
            run                : 0L,  $
            rerun              : 0L,  $
            camcol             : 0L,  $
            field              : 0L,  $
            id                 : 0L,  $
}
  
  IF keyword_set(N) THEN $
    arr = replicate(struc, N) $
  ELSE $
    arr = struc
  
  return, arr
END 
