
function sdss_flags, flagprefix

   temp = {STRUCT_FLAG, name: '', num:  0 }

   if (flagprefix EQ 'TARGET') then begin
     ; Flags for primary target
     retval = [ $
      { STRUCT_FLAG, 'QSO_HIZ'             ,  0 }, $
      { STRUCT_FLAG, 'QSO_CAP'             ,  1 }, $
      { STRUCT_FLAG, 'QSO_SKIRT'           ,  2 }, $
      { STRUCT_FLAG, 'QSO_FIRST_CAP'       ,  3 }, $
      { STRUCT_FLAG, 'QSO_FIRST_SKIRT'     ,  4 }, $
      { STRUCT_FLAG, 'GALAXY_RED'          ,  5 }, $
      { STRUCT_FLAG, 'GALAXY'              ,  6 }, $
      { STRUCT_FLAG, 'GALAXY_BIG'          ,  7 }, $
      { STRUCT_FLAG, 'GALAXY_BRIGHT_CORE'  ,  8 }, $
      { STRUCT_FLAG, 'ROSAT_A'             ,  9 }, $
      { STRUCT_FLAG, 'ROSAT_B'             , 10 }, $
      { STRUCT_FLAG, 'ROSAT_C'             , 11 }, $
      { STRUCT_FLAG, 'ROSAT_D'             , 12 }, $
      { STRUCT_FLAG, 'STBHB'               , 13 }, $
      { STRUCT_FLAG, 'STCARBON'            , 14 }, $
      { STRUCT_FLAG, 'STBROWN_DWARF'       , 15 }, $
      { STRUCT_FLAG, 'STSUB_DWARF'         , 16 }, $
      { STRUCT_FLAG, 'STCATY_VAR'          , 17 }, $
      { STRUCT_FLAG, 'STRED_DWARF'         , 18 }, $
      { STRUCT_FLAG, 'STWHITE_DWARF'       , 19 }, $
      { STRUCT_FLAG, 'SERENDIP_BLUE'       , 20 }, $
      { STRUCT_FLAG, 'SERENDIP_FIRST'      , 21 }, $
      { STRUCT_FLAG, 'SERENDIP_RED'        , 22 }, $
      { STRUCT_FLAG, 'SERENDIP_DISTANT'    , 23 }, $
      { STRUCT_FLAG, 'SERENDIP_MANUAL'     , 24 }, $
      { STRUCT_FLAG, 'QSO_FAINT'           , 25 } ]

   endif else if (flagprefix EQ 'TTARGET') then begin
     ; Flags for secondary target
     retval = [ $
      { STRUCT_FLAG, 'LIGHT_TRAP'          ,  0 }, $
      { STRUCT_FLAG, 'REDDEN_STD'          ,  1 }, $
      { STRUCT_FLAG, 'TEST_TARGET'         ,  2 }, $
      { STRUCT_FLAG, 'QA_TARGET'           ,  3 }, $
      { STRUCT_FLAG, 'SKY'                 ,  4 }, $
      { STRUCT_FLAG, 'SPECTROPHOTO_STD'    ,  5 }, $
      { STRUCT_FLAG, 'GUIDE_STAR'          ,  6 }, $
      { STRUCT_FLAG, 'BUNDLE_HOLE'         ,  7 }, $
      { STRUCT_FLAG, 'QUALITY_HOLE'        ,  8 }, $
      { STRUCT_FLAG, 'HOT_STD'             ,  9 } ]

   endif else if (flagprefix EQ 'OBJ1') then begin
     ; Photo's flags
     retval = [ $
      { STRUCT_FLAG, 'CANONICAL_CENTER'    ,  0 }, $
      { STRUCT_FLAG, 'BRIGHT'              ,  1 }, $
      { STRUCT_FLAG, 'EDGE'                ,  2 }, $
      { STRUCT_FLAG, 'BLENDED'             ,  3 }, $
      { STRUCT_FLAG, 'CHILD'               ,  4 }, $
      { STRUCT_FLAG, 'PEAKCENTER'          ,  5 }, $
      { STRUCT_FLAG, 'NODEBLEND'           ,  6 }, $
      { STRUCT_FLAG, 'NOPROFILE'           ,  7 }, $
      { STRUCT_FLAG, 'NOPETRO'             ,  8 }, $
      { STRUCT_FLAG, 'MANYPETRO'           ,  9 }, $
      { STRUCT_FLAG, 'NOPETRO_BIG'         , 10 }, $
      { STRUCT_FLAG, 'DEBLEND_TOO_MANY_PEAKS', 11 }, $
      { STRUCT_FLAG, 'CR'                  , 12 }, $
      { STRUCT_FLAG, 'MANYR50'             , 13 }, $
      { STRUCT_FLAG, 'MANYR90'             , 14 }, $
      { STRUCT_FLAG, 'BAD_RADIAL'          , 15 }, $
      { STRUCT_FLAG, 'INCOMPLETE_PROFILE'  , 16 }, $
      { STRUCT_FLAG, 'INTERP'              , 17 }, $
      { STRUCT_FLAG, 'SATUR'               , 18 }, $
      { STRUCT_FLAG, 'NOTCHECKED'          , 19 }, $
      { STRUCT_FLAG, 'SUBTRACTED'          , 20 }, $
      { STRUCT_FLAG, 'NOSTOKES'            , 21 }, $
      { STRUCT_FLAG, 'BADSKY'              , 22 }, $
      { STRUCT_FLAG, 'PETROFAINT'          , 23 }, $
      { STRUCT_FLAG, 'TOO_LARGE'           , 24 }, $
      { STRUCT_FLAG, 'DEBLENDED_AS_PSF'    , 25 }, $
      { STRUCT_FLAG, 'DEBLEND_PRUNED'      , 26 }, $
      { STRUCT_FLAG, 'ELLIPFAINT'          , 27 }, $
      { STRUCT_FLAG, 'BINNED1'             , 28 }, $
      { STRUCT_FLAG, 'BINNED2'             , 29 }, $
      { STRUCT_FLAG, 'BINNED4'             , 30 }, $
      { STRUCT_FLAG, 'MOVED'               , 31 } ]
;      { STRUCT_FLAG, 'MOVED'               , 31 }, $
;      { STRUCT_FLAG, 'DETECTED'            , 32 } ]

   endif else if (flagprefix EQ 'OBJ2') then begin
     ; Photo's flags
     retval = [ $
      { STRUCT_FLAG, 'DEBLENDED_AS_MOVING' ,  0 }, $
      { STRUCT_FLAG, 'NODEBLEND_MOVING'    ,  1 }, $
      { STRUCT_FLAG, 'TOO_FEW_DETECTIONS'  ,  2 }, $
      { STRUCT_FLAG, 'BAD_MOVING_FIT'      ,  3 }, $
      { STRUCT_FLAG, 'STATIONARY'          ,  4 }, $
      { STRUCT_FLAG, 'PEAKS_TOO_CLOSE'     ,  5 }, $
      { STRUCT_FLAG, 'MEDIAN_CENTRE'       ,  6 }, $
      { STRUCT_FLAG, 'LOCAL_EDGE'          ,  7 }, $
      { STRUCT_FLAG, 'BAD_COUNTS_ERROR'    ,  8 }, $
      { STRUCT_FLAG, 'BAD_MOVING_FIT_CHILD',  9 }, $
      { STRUCT_FLAG, 'DEBLEND_UNASSIGNED_FLUX', 10 } ]

   endif else if (flagprefix EQ 'AR_OBJECT_STATUS') then begin
     ; Flags for object status
     retval = [ $
      { STRUCT_FLAG, 'SET'                 ,  0 }, $
      { STRUCT_FLAG, 'GOOD'                ,  1 }, $
      { STRUCT_FLAG, 'DUPLICATE'           ,  2 }, $
      { STRUCT_FLAG, 'OK_RUN'              ,  3 }, $
      { STRUCT_FLAG, 'RESOLVED'            ,  4 }, $
      { STRUCT_FLAG, 'PSEGMENT'            ,  5 }, $
      { STRUCT_FLAG, 'FIRST_FIELD'         ,  6 }, $
      { STRUCT_FLAG, 'OK_SCANLINE'         ,  7 }, $
      { STRUCT_FLAG, 'OK_STRIPE'           ,  8 }, $
      { STRUCT_FLAG, 'SECONDARY'           ,  9 }, $
      { STRUCT_FLAG, 'PRIMARY'             , 10 }, $
      { STRUCT_FLAG, 'TARGET'              , 11 } ]

   endif else if (flagprefix EQ 'TYPE') then begin
     ; Flags for object type
     retval = [ $
      { STRUCT_FLAG, 'UNK'                 ,  0 }, $
      { STRUCT_FLAG, 'CR'                  ,  1 }, $
      { STRUCT_FLAG, 'DEFECT'              ,  2 }, $
      { STRUCT_FLAG, 'GALAXY'              ,  3 }, $
      { STRUCT_FLAG, 'GHOST'               ,  4 }, $
      { STRUCT_FLAG, 'KNOWNOBJ'            ,  5 }, $
      { STRUCT_FLAG, 'STAR'                ,  6 }, $
      { STRUCT_FLAG, 'TRAIL'               ,  7 }, $
      { STRUCT_FLAG, 'SKY'                 ,  8 } ]

   endif

   return, retval
end
;------------------------------------------------------------------------------
function sdss_flagname, flagprefix, flagvalue, concat=concat

   flagnames = sdss_flags(flagprefix)
   indx = where(djs_int2bin(flagvalue))
   if (indx[0] EQ -1) then begin
      retval = ''
   endif else begin
      retval = flagnames[indx].name
   endelse

   if (keyword_set(concat)) then begin
      for i=1, n_elements(retval)-1 do $
       retval[0] = retval[0] + ' ' + retval[i]
      retval = retval[0]
   endif

   return, retval
end
;------------------------------------------------------------------------------
