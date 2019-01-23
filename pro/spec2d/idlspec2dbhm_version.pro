;+
; NAME:
;   idlspec2dbhm_version
; PURPOSE:
;   Return the version name for the black hole mapper idlspec2d product 
; CALLING SEQUENCE:
;   vers = idlspec2dbhm_version()
; OUTPUTS:
;   vers       - Version name for the product idlspec2d 
; COMMENTS:
;   Depends on shell script in $IDLSPEC2D_DIR/bin
;-
;------------------------------------------------------------------------------
function idlspec2dbhm_version
   spawn, 'idlspec2d_version', stdout, /noshell
   return, stdout[0]
end
;------------------------------------------------------------------------------
