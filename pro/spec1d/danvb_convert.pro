; Convert the file "plateEDR.summ" from Dan van den Berk from ASCII
; format to a FITS binary table, "plateEDR.summ.fits".
pro danvb_convert

   filename = 'plateEDR.summ'
   outfile = 'plateEDR.summ.fits'

   readcol, filename, plate, fiber, mjd, id_manual, flag_manual, ra, dec, $
    run, camcol, field, u_fib, g_fib, r_fib, i_fib, z_fib, z_manual, $
    z_1d, zerr_1d, z_confidence, z_status, class_1d, primtarget, $
    sectarget, $
    format='(L,L,L,A,L,D,D,L,L,L, F,F,F,F,F, F,F,F,L,L,L,L)', skipline=2
   readfmt, filename, '99X,70X,A100', comments, skipline=2
   comments = strtrim(comments,2)

   danstruct = create_struct( $
    'plate'       ,  0L, $
    'fiber'       ,  0L, $
    'mjd'         ,  0L, $
    'id_manual'   , ' ', $
    'flag_manual' ,  0L, $
    'ra'          ,  0D, $
    'dec'         ,  0D, $
    'run'         ,  0L, $
    'camcol'      ,  0L, $
    'field'       ,  0L, $
    'u_fib'       ,  0., $
    'g_fib'       ,  0., $
    'r_fib'       ,  0., $
    'i_fib'       ,  0., $
    'z_fib'       ,  0., $
    'z_manual'    ,  0., $
    'z_1d'        ,  0., $
    'zerr_1d'     ,  0., $
    'z_confidence',  0., $
    'z_status'    ,  0L, $
    'class_1d'    ,  0L, $
    'primtarget'  ,  0L, $
    'sectarget'   ,  0L, $
    'comments'    ,  ' ')
   danstruct = replicate(danstruct, n_elements(plate))
   danstruct.plate        = plate
   danstruct.fiber        = fiber
   danstruct.mjd          = mjd
   danstruct.id_manual    = id_manual
   danstruct.flag_manual  = flag_manual
   danstruct.ra           = ra
   danstruct.dec          = dec
   danstruct.run          = run
   danstruct.camcol       = camcol
   danstruct.field        = field
   danstruct.u_fib        = u_fib
   danstruct.g_fib        = g_fib
   danstruct.r_fib        = r_fib
   danstruct.i_fib        = i_fib
   danstruct.z_fib        = z_fib
   danstruct.z_manual     = z_manual
   danstruct.z_1d         = z_1d
   danstruct.zerr_1d      = zerr_1d
   danstruct.z_confidence = z_confidence
   danstruct.z_status     = z_status
   danstruct.class_1d     = class_1d
   danstruct.primtarget   = primtarget
   danstruct.sectarget    = sectarget
   danstruct.comments     = comments

   mwrfits, danstruct, outfile, /create
end
