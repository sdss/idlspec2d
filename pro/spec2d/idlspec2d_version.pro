;+
; NAME:
;   idlspec2d_version
; PURPOSE:
;   Return the version name for the idlspec2d product
; CALLING SEQUENCE:
;   vers = idlspec2d_version()
; OUTPUTS:
;   vers       - Version name for the product idlspec2d
; COMMENTS:
;   Requires that the IDLSPEC2D_DIR environment variable be set
; VERSION:
;   $Id$
;-
;------------------------------------------------------------------------------
FUNCTION idlspec2d_version
    RETURN, FILE_BASENAME(GETENV('IDLSPEC2D_DIR'))
END
;------------------------------------------------------------------------------
