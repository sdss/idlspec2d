;+
; NAME:
;   gen_scatter
;
; PURPOSE:
;   Generate a scattering surface for a given image. This function is
;   used both as an MPFIT function and to generate the surface.
;
; CALLING SEQUENCE:
;   gen_scatter(img, lambdas, params, $
;               evalfunc=, $
;               bandrows=, maxkernr=, dosmall=)
;
; INPUTS:
;   img        - an image to generate a scattering surface
;                from. [NCOLS, NROWS]
;   lambdas    - the wavelengths of the rows, in _microns_. There can
;                be more wavelengths than rows, in which case the
;                surface will be evaluated for rows "above" the image,
;                using the flux in the top row of the image. [NCOLS, NLAMBDAS]
;   params     - fit parameters, passed to evalfunc. These are the
;                coefficients which MPFITFUN might be
;                modifying/fitting. [NPARAMS]
;   evalfunc   - a function which generates a model surface given an
;                array of distances, a wavelength, and the parameters
;
; OPTIONAL KEYWORDS:
;   bandrows   - how many rows of the image to evaluate at once. The
;                parameters for the middle row are used. This gives a
;                direct linear speedup. 
;   maxkernr   - The maximum allowed kernel radius. Defaults to 256.
;   dosmall    - Only evaluate a region around the middle of the
;                image. The region size is scaled to the size of the
;                scattering model; this lets you avoid evaluating a
;                large region, and is only really useful for fitting
;                to an extremely sparse flat. And is only needed
;                because the IDL FFT is wretchedly slow, even if you are
;                careful about input sizes.
; OUTPUT:
;   the evaluated scattering surface.
;
; COMMENTS:
;   The slightly odd arrangement of function arguments is because the
;   function has to be callable by mpfitfun.
;
;   The peculiar micro-optimizations (and traces of older ones) are
;   because calling this with a large kernel per row from mpfitun does
;   actually take time.
;
;   When banding, it would be better to linearly interpolate between
;   the surfaces generated with the endpoint parameters.
;
;   The wavelengths are in microns mostly for historical reasons.
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;   redhalo     - generates red side backscattering halo
;   bluescat    - generates blue side wing + gaussian scattering
;     Both take the following:
;      lam     - the wavelength, in microns, of the center row.
;      par     - an array of parameters.
;      r       - an array of distances.
;
;   nextsize_for_IDL_fft   - currently unused.
;   nextsize_for_IDL_2dfft - currently unused.
;

function redhalo, lam, par, r
                                ; See
                                ; $IDLSPEC2D_DIR/doc/IR_scattering.txt
                                ; for the form of the function. We are
                                ; using different coefficients.
    r0 = par[2] * lam^2
    f = exp(par[3] * (lam-1.05))
    halo = (f/(2*!PI*r0)) * exp(-r/r0)/(r>1)

                                ; Notch out flux in the real core
    if par[1] GT 0 then begin
       z = 0.5*(r/par[1])^2
       icore = max(halo) * exp(-(z < 50))
       return, halo - icore
    end else begin
       return, halo
    end
 end

function bluescat, lam, par, r
    sigma = par[2]
    b = par[3]
    p = par[4]
                                ; The wings
    z = 0.5*(r/sigma)^2
    wing = b / (1. + z/p)^p

                                ; The wide second core gaussian.
    sigma2 = par[5]
    a2 = par[6]
    z2 = 0.5*(r/sigma2)^2
    core2 = a2 * exp(-(z2 < 50))

                                ; Notch out flux in the real core
    if par[1] GT 0 then begin
       model = wing + core2
       z3 = 0.5*(r/par[1])^2
       icore = max(model) * exp(-(z3 < 50))
       return, model - icore
    end else begin
       return, wing+core2
    end
 end


; Return the next even number with no prime factors other than 2,3,5
;
function nextsize_for_IDL_fft, n
    repeat begin
       faclist = factorize(n)
       flist = where(faclist GT 5, fcnt) 
       if (fcnt EQ 0) AND (n mod 2 EQ 0) then $
          return, n 
       n += 1
    endrep until 0
end

; Return the next even number N with no prime factors other than 2,3,5,
; where (N + ax0) also meets the rules.
;
function nextsize_for_IDL_2dfft, n, ax0
    if ax0 MOD 2 eq 1 then message, "ax0 must be even"
    repeat begin
       if (n mod 2 eq 1) then begin
          n += 1
          continue
       end

       faclist1 = factorize(n)
       flist1 = where(faclist1 GT 5, fcnt1) 
       if (fcnt1 EQ 0) then begin
          faclist2 = factorize(ax0 + n)
          flist2 = where(faclist2 GT 5, fcnt2) 
          if (fcnt2 EQ 0) then $
             return, n 
       end

       n += 1
    endrep until 0
 end

function gen_scatter, img, lambdas, params, evalfunc=evalfunc, $
                      bandrows=bandrows, maxkernr=maxkernr, dosmall=dosmall, $
                      keeprows=keeprows, doplot=doplot
                      
    imgsize = size(img,/dim)
    nr = imgsize[1]
    nc = imgsize[0]

    if not keyword_set(bandrows) then bandrows = 1

    if nc mod 2 EQ 1 then message,'number of columns must be even'
    if nr mod 2 EQ 1 then message,'number of rows must be even'

    if size(params, /n_dim) NE 1 and n_elements(params[0,*]) NE nr then $
       message,'params must be a vector or have the same number of rows as the image'

    splog, "params: ", params

                                ; How big can the kernel possibly get?
                                ; If banding is used, that should be
                                ; allowed for by setting maxkernr big
                                ; enough.
    maxneedr = keyword_set(maxkernr) ? maxkernr : 256
    maxneedw = 2 * maxneedr
                                ; Make a result image large enough not
                                ; to need to clip, etc.  Build a
                                ; corresponding scratch band to apply
                                ; the kernel to.
    maxw = keyword_set(dosmall) ? maxneedw : nc + maxneedw
    res = fltarr(nc+maxneedw, nr+maxneedw) 
    res += params[0]            ; DC term
    band = fltarr(maxw, maxneedw)

                                ; Make bookkeeping easier by
                                ; constructing an input image large
                                ; enough to hold the any extra
                                ; rows. Beyond the detector we have no
                                ; idea what the actual incoming flux
                                ; is, so use what we saw on the last
                                ; row.
    timg = fltarr(nc, nr)
    timg[0, 0] = img

                                ; calculate the distance map once
    xs = rebin(findgen(maxw) - maxw/2, [maxw,maxneedw])
    ys = rebin(transpose(findgen(maxneedw) - maxneedr), [maxw, maxneedw])
    r = sqrt(xs^2 + ys^2)
    rcenter = [maxw/2, maxneedw/2]
    if r[rcenter[0], rcenter[1]] ne 0 then message, "rcenter oops"

    nbands = nr / bandrows
    for b_i=0,nbands-1 do begin
       t0 = systime(1)
                                ; the parameters for this band/row
       r_i = b_i * bandrows + bandrows/2
       bparams = size(params, /n_dim) eq 1 ? params : params[*,r_i]

       band_r0 = b_i * bandrows
       band_r1 = (b_i+1)*bandrows - 1
       bandh = band_r1-band_r0+1
                                ; How big a kernel do we need for
                                ; lam[r_i]? This is the biggest
                                ; controllable factor for the running
                                ; time.
       lam = lambdas[r_i]
       rx = findgen(maxneedw)
       thishalo = call_function(evalfunc, lam, bparams, rx)
       needr = max(where(thishalo GT 1e-7)) > 16; Do this better - CPL
       needr = needr > bandrows
                                ; generate the model surface for this
                                ; row and extract the image row into
                                ; an otherwise empty padded buffer.

                                ; The IDL FFT can do non-powers of 2,
                                ; but the running time is proportional
                                ; to the sum of the prime factors of
                                ; each (?) array dimension. Crudely
                                ; avoid punitive sizes. [Discovery:
                                ; the IDL FFT is poor. ]
;if keyword_set(dosmall) then $
;   needw = (nextsize_for_IDL_fft(needr*2) < 2*maxneedr) $
;else $
;   needw = (nextsize_for_IDL_2dfft(needr*2, nc) < 2*maxneedr)
;needr = needw/2
       needr = (32 * ((needr + 31)/32)) < maxneedr
       needw = 2*needr

       inset = maxneedr - needr
       bandrow0 = needr - bandh/2
       resrow0 = inset + band_r0 + bandh/2

       if keyword_set(dosmall) then begin
                                ; Solve for the halo region only,
                                ; centered on the middle
                                ; column. Sometimes necessary for
                                ; fitting. With the FFTs, needr better
                                ; be plenty large.
          imgcols = indgen(needw) + nc/2 - needr
          bandcol0 = 0
          rescol0 = nc/2 + inset 
          thisr = r[rcenter[0]-needr:rcenter[0]+needr-1, rcenter[1]-needr:rcenter[1]+needr-1]
       end else begin
                                ; Solve for the whole row. Used for
                                ; generating the scattering surfaces.
          imgcols = indgen(nc)
          bandcol0 = needr
          rescol0 = inset
          thisr = r[inset:inset+nc+needw-1, inset:inset+needw-1]
       end 
       rowhalo = call_function(evalfunc, lam, bparams, thisr)
       
                                ; Insert the image rows into the
                                ; scratch band.
       band = thisr * 0
       band[bandcol0, bandrow0] = timg[imgcols, band_r0:band_r1]
       bandsize = size(band, /dim)
       
                                ; Convolve the row band with the
                                ; kernel, sum into the full result.
                                ; Basically: 
                                ; rowres = convol(band, rowhalo, /edge_xxxx)
       t1 = systime(1)
       halofft = fft(rowhalo)
       bandfft = fft(band)
       xx = real_part(fft(halofft*bandfft,/inverse))
                                ; The IDL FFT doesn't scale the
                                ; inverse. Do that here.
       rowres = shift(xx, bandsize/2) * bandsize[0] * bandsize[1]
       res[rescol0:rescol0+bandsize[0]-1, resrow0:resrow0+bandsize[1]-1] += rowres

                                ; Depressingly useful counter of step number...
       print, format='("Step ",i5," of ",i6,f5.2,i4,f5.2,f5.2, a1,$)', $
              r_i, nr, lam, needr, t1-t0, systime(1)-t0, string(13b)
       
    end

    keepr0 = maxneedr
    keepr1 = keepr0+nr-1
    if keyword_set(keeprows) then begin
       keepr0 += keeprows[0]
       keepr1 = keepr0 + (keeprows[1]-keeprows[0])
    end

    rres = res[maxneedr:nc+maxneedr-1, keepr0:keepr1]
    if keyword_set(doplot) then $
       pres,img,rres,rows=[(keepr1-keepr0)/2]
    return, rres
 end
