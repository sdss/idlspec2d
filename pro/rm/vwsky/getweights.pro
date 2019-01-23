;
;+
; NAME:
;	getweights.pro

;
; PURPOSE:
;	Collect median error array for given plate from file and
;	weight to account for SDSS scaling of Poisson errors.
;
;
; CALLING SEQUENCE:
;
;	weights = GETWEIGHTS(Plate, Dir)
;
;
; INPUTS:
;	Plate:      SDSS plate number (integer), can be array of
;	            plate numbers.
;       Dir:        String giving directory of weight file
;
; OUTPUTS:
;	This function returns an array containing the appropriate
;	weighting for this plate number(s).
;
;               
; MODIFICATION HISTORY: Vivienne Wild, vw@ast.cam.ac.uk, 05/01/05
;
;****************************************************************************************

FUNCTION GETWEIGHTS, Plate, DIR

;*** read in plate noise file
;contains parameters plateerror, continuum, plateid
RESTORE, DIR+'plate_weights.sav'
nbin = (SIZE(plateerror,/dim))[0]

nplates = N_ELEMENTS(plate)

;*** set parameters by which the SDSS noise arrays are altered
maxscale = 0.7
alpha = 1.

weights = FLTARR(nbin,nplates)

for i = 0, nplates -1 do begin

;*** find correct plate in array
    ind_plate = WHERE(plateid eq plate[i], count) 
    if NOT(count) then MESSAGE, 'Plate not found - please supply own weights'

    noise = plateerror[*,ind_plate] 
    cont = continuum[*,ind_plate]

;*** scale error
    maxnoise = MAX(noise-cont)
    index = WHERE(abs(noise-cont) gt 0.02 and noise-cont gt 0) 
    scale = FLTARR(nbin)+1
    scale[index] = (1-( ((noise[index]-cont[index])/maxnoise)^alpha * (1-maxscale)) ) < 1 ;max value=1

;*** calculate weights

    posvar = WHERE(noise ne 0)
    weights[posvar,i] = 1./(noise[posvar]*scale[posvar])

endfor

RETURN, weights

END
