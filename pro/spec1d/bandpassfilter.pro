;+
; NAME:
;   bandpassfilter
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   newdata = bandpassfilter( datafft, [ klo_cut=, khi_cut= ] )
;
; INPUTS:
;   datafft    - Vector of Fourier-transformed data
;
; OPTIONAL KEYWORDS:
;   klo_cut    - Low-frequency cutoff
;   khi_cut    - High-frequency cutoff
;
; OUTPUTS:
;   newdata    - Filtered version of DATAFFT
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Units???
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   30-Mar-2000  Written by D. Schlegel, APO
;-
;------------------------------------------------------------------------------
function bandpassfilter, datafft, klo_cut=klo_cut, khi_cut=khi_cut

   if (size(datafft, /n_dimen) NE 1) then $
    message, 'DATAFFT is not 1-dimensional'

   if (NOT keyword_set(klo_cut) AND NOT keyword_set(khi_cut)) then $
    return, datafft

   if (size(datafft, /tname) EQ 'DOUBLE') then PI = !dpi $
    else PI = !pi

   ndata = N_elements(datafft)
   knums = fft_wavenums(ndata)
   hipass = fltarr(ndata) + 1.0
   lopass = fltarr(ndata) + 1.0

   if (keyword_set(klo_cut)) then begin
      ii = where(abs(knums) LT klo_cut)
      if (ii[0] NE -1) then begin
         hipass[ii] = 0.5 * (1.0 - cos(PI*abs(knums[ii])/klo_cut))
      endif
   endif

   if (keyword_set(khi_cut)) then begin
      ii = where(abs(knums) GT khi_cut)
      if (ii[0] NE -1) then begin
         sep = max(knums) - khi_cut
         lopass[ii] = $
          0.5 * (1.0 - cos(PI*(max(knums) - abs(knums[ii]))/sep))
      endif
   endif

   return, datafft * hipass * lopass
end
;------------------------------------------------------------------------------
