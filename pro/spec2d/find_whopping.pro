;+
; NAME:
;   find_whopping
;
; PURPOSE:
;   Use simple medians to detect whopping fibers
;   Take care to not choose adjacent fibers
;
; CALLING SEQUENCE:
;   whopping = find_whopping(flux, thresh, whopct)
;
; INPUTS:
;   flux       - extracted spectra [npix, nfiber]
;   thresh     - minimum threshold for median counts in whopping fiber
;                 Default 10000.0
;
; OUTPUTS:
;   whopping   - indices of whopping fiber(s) in flux
;
; OPTIONAL OUTPUTS:
;   whopct     - Number of whopping fibers detected
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_median()
;
; REVISION HISTORY:
;   24-Jan-2000  Written by S. Burles, Chicago
;    4-Apr-2001  Moved to its own routine for use by SoS
;-
;------------------------------------------------------------------------------

function find_whopping, flux, thresh, whopct

    if NOT keyword_set(thresh) then thresh=10000.0

    scrunch = djs_median(flux, 1) ; Find median counts/row in each fiber
    scrunch_sort = sort(scrunch)
    i5 = n_elements(scrunch)/20
    i95 = i5 * 19
    boxcar = scrunch - djs_median(flux) ; Find median counts/row in all fibers

    candidates = where(boxcar GT thresh, whopct)

    if (whopct LT 2) then return, candidates

    testc = [-20, candidates, n_elements(boxcar)+20]
    diff = testc[1:whopct+1] - testc[0:whopct]
    bottom = where(diff[0:whopct-1] NE 1, bottomct)
    top = where(diff[1:whopct] NE 1,topct)

    if (topct NE bottomct) then begin
      message, 'Bug introduced by Scott, look in find_whopping in extract_object.pro' ; ???
    endif

    whopping = lonarr(topct)
    for i=0, topct -1 do begin
       mmax = max(boxcar[bottom[i]:top[i]], place)
       whopping[i] = candidates[place + bottom[i]]
    endfor

    whopct = topct
    return, whopping
end

