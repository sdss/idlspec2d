;+
;
; NAME: gausspix
;
; PURPOSE: pixel-integrated 1D Gaussian.
;
; USAGE:
;  profile = gausspix(x, par [, /skew, /kurt])
;
; ARGUMENTS:
;  x: pixel-centered pixel coordinate over which to evaluate
;  par: parameter vector:
;   par[0] = centroid, in pixels, same frame as x
;   par[1] = sigma value of non-integrated profile,
;            in pixels
;  *note*: passing "par" of length greater than 2 will
;          generate a basis array, using each successive
;          pair of par entries to specify each basis row.
;   skew: set this keyword to return not the pixel-integrated
;         Gaussian, but rather a pixel-integrated skewing term:
;         Gaussian times (u^3 - 3*u), where u = (x - xcen) / sigma.
;   kurt: set this keyword to return not the pixel-integrated
;         Gaussian, but rather a pixel-integrated kurtosis term:
;         Gaussian times (u^4-6*u^2+3), where u = (x - xcen) / sigma.
;
; COMMENTS:
;   Amplitude is always (integrated) unity.  If you want to dial
;   in fixed amplitudes, wrap this function in something else.
;
;   Skewness and kurtosis bases are pixel-integrated Hermite
;   functions, and have no net flux associated with them.
;
; WRITTEN:
;   Adam S. Bolton, U. of Utah, 2010 May
;
;-

function gausspix, x, par, skew=skew, kurt=kurt
; Output is float or double depending upon inputs:
  xtype = size(x, /type)
  ptype = size(par, /type)
  otype = ((xtype > ptype) > 4) < 5
  nx = n_elements(x)
  ngauss = n_elements(par) / 2L
  outarray = make_array(nx, ngauss, type=otype)
  for i = 0L, ngauss-1 do begin
     xcen = par[2L*i]
     sigma = abs(par[2L*i+1L])
     uhi = (x - xcen + 0.5d0) / sigma
     ulo = (x - xcen - 0.5d0) / sigma
     if ((not keyword_set(skew)) and (not keyword_set(kurt))) then begin
        outarray[*,i] = gaussint(uhi) - gaussint(ulo)
     endif else begin
        if keyword_set(skew) then outarray[*,i] = (1.d0/sqrt(2.d0*!dpi)) * $
         ( (1.d0 - uhi^2) * exp(-0.5d0 * uhi^2) - (1.d0 - ulo^2) * exp(-0.5d0 * ulo^2))
        if keyword_set(kurt) then outarray[*,i] = (1.d0/sqrt(2.d0*!dpi)) * $
         ( (3.d0*uhi - uhi^3) * exp(-0.5d0 * uhi^2) - (3.d0*ulo - ulo^3) * exp(-0.5d0 * ulo^2))
     endelse
;     print,sigma, xcen
  endfor
  return, outarray
end
