;------------------------------------------------------------------------------
; This is a routine to rebin Eisenstein's stellar templates for use
; as velocity dispersion standards.  This is basically a throw-away
; routine, until we properly generate such templates.
pro vdisprebin

   readcol, 'sp4', lam, s1, s2, s3, s4, format='(D,D,D,D)'
   readcol, 'sp321', lam, s5, s6, s7, format='(D,D,D,D)'
   newflux  = [[s1],[s2],[s3],[s4],[s5],[s6],[s7]]

   ;----------
   ; Write output file

   loglam = alog10(lam)
   dloglam = 1.d-4

   sxaddpar, hdr, 'COEFF0', loglam[0]
   sxaddpar, hdr, 'COEFF1', dloglam
   mwrfits, float(newflux), 'spEigenVdisp.fits', hdr, /create

   return
end
;------------------------------------------------------------------------------
