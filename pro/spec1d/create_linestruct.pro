;------------------------------------------------------------------------------
function create_linestruct, nline

   ;----------
   ; Initialize the output structure, and set default return values

   linestruct = create_struct( name='EMLINEFIT',   $
    'linename'          ,   ' ', $
    'linewave'          ,  0.0d, $
    'linez'             ,   0.0, $
    'linez_err'         ,   0.0, $
    'linesigma'         ,   0.0, $
    'linesigma_err'     ,   0.0, $
    'linearea'          ,   0.0, $
    'linearea_err'      ,   0.0, $
    'lineew'            ,   0.0, $
    'lineew_err'        ,  -1.0, $
    'linecontlevel'     ,   0.0, $
    'linecontlevel_err' ,  -1.0, $
    'linenpixleft'      ,   0L,  $
    'linenpixright'     ,   0L,  $
    'linedof'           ,   0.0, $
    'linechi2'          ,  -1.0  )
   if (keyword_set(nline)) then $
    linestruct = replicate(linestruct, nline)

   return, linestruct
end
;------------------------------------------------------------------------------
