;+
; NAME:
;   skyline_dispersion
;
; PURPOSE:
;   Find median skyline width relative to arcline dispset for each line
;
; CALLING SEQUENCE:
;   dispset = skyline_dispersion(flux, fluxivar, xcen, iskies, dispset)
;
; INPUTS:
;   flux        - Object extracted flux
;   fluxivar    - corresponding inverse variance
;   xcen        - xpeaks of sky lines [ntrace,nlines]
;   iskies      - fibers corresponding to good sky fibers
;   dispset     - Original arcline solution for resolution vs pixel
;
; OPTIONAL KEYWORDS:
;
; OUTPUTS:
;   dispset   - a traceset structure containing adjusted coefficients
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_iterstat
;   extract_image
;   splog
;   traceset2xy
;   xy2traceset
;
; REVISION HISTORY:
;   20-Feb-2002  Written by S. Burles, MIT (Adopted from fitdispersion)
;-
;------------------------------------------------------------------------------
function skyline_dispersion, flux, fluxivar, xcen, iskies, dispset

        if (NOT keyword_set(xmin)) then xmin = 0.0
        if (NOT keyword_set(xmax)) then xmax = 2047.0

        nline = (size(xcen,/dimen))[1]
        ntrace = (size(xcen,/dimen))[0]
        if n_elements(iskies) EQ 0 then iskies = lindgen(ntrace)
  

        ;  sort positions for extract_image

        hmm = sort(xcen[0,*])
        xsky = xcen[*,hmm]

        ;  and make mask for extraction
        ysky = lindgen(ntrace) # replicate(1,nline) 
        npix = (size(fluxivar,/dimen))[0]
        mask = long(fluxivar) * 0L

        for offset = -12,12 do begin
          xtemp = (long(xsky + offset) > 0) < (npix - 1L)
          mask[xtemp,ysky] = 1
        endfor

        ; Extract sky lines

        traceset2xy, dispset, transpose(xsky), sigma
        sigma=transpose(sigma)

        extract_image, flux, fluxivar*mask, xsky, sigma, yrow=iskies, $
          skylineflux, skylineivar, ansimage=ansimage, wfixed=[1,1], $
          highrej=10, lowrej=10, relative=1, npoly=5, proftype=1, $
          ymodel=ymodel
         
        ;
        ;  Prepare width terms
        ;  

       splog, '  xpos  ngood flux   arc    min    max  median med_diff stddev'

        finaldiff = fltarr(nline)

        for iline=0, nline-1 do begin
          mask2 = transpose((skylineivar[*,iline] GT 0) $
           * (skylineflux[*,iline] GT 0))
          igood = where(mask2, ngood)

          if (ngood GT 0) then begin
             widthterm = ansimage[iline*2+1,igood]
             width = (1 + widthterm / skylineflux[igood,iline]) $
              * sigma[igood,iline]
             quad_diff = width^2 - sigma[igood,iline]^2
             djs_iterstat, quad_diff, median=md, sigma=ss

             splog, median(xsky[igood,iline]), ngood, $
              median(skylineflux[igood,iline]), $
               median(sigma[igood,iline]), $
               min(width), max(width), median(width), $
               md, ss, format='(f8.3, i4, i6, 6f7.3)'

             finaldiff[iline] = md
           endif
        endfor
    
        diff25 = finaldiff[(sort(finaldiff))[nline/4]] > 0
        splog, 'Adjusting dispset with 25% difference: ', sqrt(diff25)
 
        traceset2xy, dispset, x, sigma
        sigma = sqrt(sigma^2 + diff25)
        ncoeff = (size(dispset.coeff))[1]

        xy2traceset, x, sigma, skydispset, ncoeff=ncoeff, $
           xmin=dispset.xmin, xmax=dispset.xmax, yfit=yfit

        return, skydispset
end
