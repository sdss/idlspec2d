
;------------------------------------------------------------------------------
function design_starstruct, num
   
   result = create_struct( $
    name = 'STARDATA', $
    'ra'        ,  0.d, $
    'dec'       ,  0.d, $
    'mag'       , fltarr(5), $
    'holetype'  ,   '', $
    'objtype'   ,   '', $
    'priority'  ,   0L, $
    'primtarget',   0L, $
    'sectarget' ,   0L, $
    'tilenum'   ,   0L, $
    'comments'  ,   '' )
   if (keyword_set(num)) then result = replicate(result, num)

   return, result
end
;------------------------------------------------------------------------------

