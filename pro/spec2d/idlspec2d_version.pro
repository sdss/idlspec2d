;+
; NAME:
;   idlspec2d_version
;
; PURPOSE:
;   Return the version name for the product idlspec2d
;
; CALLING SEQUENCE:
;   vers = idlspec2d_version()
;
; INPUTS:
;
; OUTPUTS:
;   vers       - Version name for the product idlspec2d
;
; COMMENTS:
;   If this version is not tagged by CVS, then we return 'NOCVS: TOPLEVEL'
;   where TOPLEVEL is the last directory in the environment variable
;   $IDLSPEC2D_DIR.  For example, if you are using a version of the code
;   in the directory '/u/schlegel/idlspec2d/v0_0', then this returns
;   'NOCVS:v0_0'.
;
; BUGS:
;
; PROCEDURES CALLED:
;
;; REVISION HISTORY:
;   01-Dec-1999  Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------

function idlspec2d_version

   ; The following expression in dollar signs is expanded by CVS
   ; and replaced by the tag name for this version.
   name = '$Name$'

   words = str_sep(strcompress(name), ' ')

   if (words[0] EQ '$Name:' AND N_elements(words) EQ 3) then begin
      vers = words[1]
   endif else begin
      dirname = getenv('IDLSPEC2D_DIR')
      if (dirname NE '') then begin
         words = str_sep(dirname,'/')
         nword = N_elements(words)
         vers = 'NOCVS:' + words[nword-1]
      endif else begin
         vers = 'NOCVS:Unknown'
      endelse
   endelse

   return, vers
end
;------------------------------------------------------------------------------
