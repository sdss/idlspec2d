;------------------------------------------------------------------------------
function psolve_obj_struct, nobj

   obj1 = create_struct( $
    'XCEN', 0., $
    'XCEN_ERR', 0., $
    'YCEN', 0., $
    'YCEN_ERR', 0., $
    'FLUX', 0., $
    'FLUX_ERR', 0., $
    'QSATUR', 0b, $
    'GOODMASK', 0b, $
    'BESTMASK', 0b, $
    'FITSTATUS', 0, $
    'FITITER', 0, $
    'CHI2', 0., $
    'DOF', 0., $
    'RMOMENT', 0., $
    'NEFF', 0., $
    'NERR', 0., $
    'PVALUE', fltarr(16) ) ; This should really be dimensioned by nmap !!!???
   if (keyword_set(nobj)) then return, replicate(obj1, nobj) $
    else return, obj1

end
;------------------------------------------------------------------------------
