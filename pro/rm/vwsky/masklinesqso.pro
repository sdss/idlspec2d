;
;+
; NAME:
;	masklinesqso.pro

;
; PURPOSE:
;	Identify and mask emission/absorption features in a QSO
;	spectrum. Differs from masklines.pro in having maximum
;	linemask size depending on filter size used in median filter.
;
;
; CALLING SEQUENCE:
;
;	pixel_mask = MASKLINES_QSO(Flux, Wave, Z, Linefile, Nfilt, /SILENT)
;
;
; INPUTS:
;	Flux:       spectrum
;       Wave:       observed wavelength array (Angstroms)
;       Z:          redshift of QSO
;       Linefile:   Path and name of file containing array of
;                   wavelengths to be masked.
;       Nfilt:      Size of median filter used on spectrum before masking.
;
; KEYWORD PARAMETERS: 
;       SILENT:     Stop program writing warnings
;
; OUTPUTS:
;	This function returns an array the same length of flux where 1
;	identifies a pixel with a line, and 0 a pixel without a line.
;
;               
; MODIFICATION HISTORY: Vivienne Wild, vw@ast.cam.ac.uk, 05/01/05
;
;****************************************************************************************
FUNCTION MASKLINESQSO, Flux, Wave, Z, Linefile, Nfilt, silent=silent

if N_ELEMENTS(wave) ne N_ELEMENTS(flux) then message, 'flux and wave must be same size'

;read in line file: 
;Note the broad lines have a single wavelength, the width of the mask
;is then set by the median filter size. The narrow lines have an upper
;and lower wavelength

OPENR, lun, linefile, /get_lun
nlines = 0
READF, lun, nlines  
emlines = FLTARR(2,nlines)
READF, lun, emlines
FREE_LUN, lun

;create pixel mask
nbin = N_ELEMENTS(wave)
pixel_mask = FLTARR(nbin)

for i=0, nlines-1 do begin
    ind = WHERE(wave/(1+z) ge emlines[0,i], count)

    ;the approx. size of the median filter in angstroms, observed frame
    if count gt 2 then a = nfilt*(wave[ind[1]]-wave[ind[0]]) else continue 

    if emlines[0,i] eq emlines[1,i] then begin ;the broad lines
        emlines[0,i] = emlines[0,i]*(1+z)-a/3.
        emlines[1,i] = emlines[1,i]*(1+z)+a/3.
    endif else begin            ;the narrow lines
        emlines[0,i] = emlines[0,i]*(1+z)
        emlines[1,i] = emlines[1,i]*(1+z)
    endelse
    
    index = where(wave gt emlines[0,i] and wave lt emlines[1,i],count) 
    if count gt 0 then pixel_mask[index] = 1.

endfor


;warnings
junk = WHERE(pixel_mask eq 1,n)
if n eq 0 and NOT(KEYWORD_SET(silent)) then print, 'masklines_qso.pro: no qso lines masked'

return, pixel_mask

end
