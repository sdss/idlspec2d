;+
; NAME:
;   fitdispersion
;
; PURPOSE:
;   Fit polynomials to the line width (dispersion axis) for each fiber bundle
;
; CALLING SEQUENCE:
;   dispset = fitdispersion(arc_flux, arc_fluxivar, xcen, $
;             sigma=sigma, ncoeff=ncoeff, xmin=xmin, xmax=xmax)
;
; INPUTS:
;   arc_flux     - arc image extracted flux
;   arc_fluxivar - corresponding inverse variance
;   xcen         - xpeaks of arc lines [ntrace,nlines]
;
; OPTIONAL KEYWORDS:
;   ncoeff     - order of legendre polynomial to apply to width vs. row
;   sigma      - The initial sigma (in pixels, 1.0 is a good guess)
;   xmin, xmax - limits in traceset
;
; OUTPUTS:
;   dispset   - a traceset structure containing fitted coefficients
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;    Used to fill arcstruct.dispset, which can then be applied
;     to psf-corrected sky subtraction 
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   1-Mar-2000  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function fitdispersion, arc_flux, arc_fluxivar, xcen, $
             sigma=sigma, ncoeff=ncoeff, xmin=xmin, xmax=xmax

        if (NOT keyword_set(sigma)) then sigma = 1.0
        if (NOT keyword_set(ncoeff)) then ncoeff = 4
        if (NOT keyword_set(xmin)) then xmin = 0.0
        if (NOT keyword_set(xmax)) then xmax = 2047.0

        nline = (size(xcen,/dimen))[1]
        ntrace = (size(xcen,/dimen))[0]

        if (ntrace NE 320 OR nline LT 10) then begin
            splog, 'WARNING: Could not find dispersion traceset'
            return, -1
        endif

        ;
        ;	Put blue arcs in correct order for extraction
        ;
        if (xcen[0,0] GT xcen[0,1]) then xcentemp = reverse(xcen,2) $
         else xcentemp = xcen

	;
        ; Extract arc lines  	
        ;
        extract_image, arc_flux, arc_fluxivar, xcentemp, sigma, $
          arclineflux, arclineivar, ansimage=ansimage, wfixed=[1,1], $
          highrej=10, lowrej=10, relative=1

         
        ;
        ;  Prepare width term
        ;  


        mask = transpose((arclineivar GT 0) * (arclineflux GT 0))
	
        good = where(mask)
	width = transpose(arclineflux * 0.0)

        widthterm = ansimage[lindgen(nline)*2+1,*]
	width[good] = widthterm[good]/(transpose(arclineflux))[good]

	width = reform(width,nline,20,16)
	mask = reform(mask,nline,20,16)

        width_bundle = fltarr(nline,16)

;
;	Perform median across bundles on good arclines only
;       somewhat tedious, but it works
;
        for i = 0, nline - 1 do begin
          for j = 0, 15 do begin
             ss = where(mask[i,*,j])  
             if (ss[0] NE -1) then $
               width_bundle[i,j] = djs_median(width[i,ss,j]) 
          endfor
        endfor

        expand_width = rebin(width_bundle,nline,320,/sample) + sigma

        xy2traceset, transpose(xcentemp), expand_width, dispset, $
             ncoeff=ncoeff, xmin=xmin, xmax=xmax, yfit=yfit

        return, dispset
end
