; Fit a scattered light image to the 
;+
; NAME:
;   fitscatter
;
; PURPOSE:
;   Simple scattered light fitting on BOSS spectro images
;
; CALLING SEQUENCE:
;   scatimg = fitscatter( image, invvar, xcen)
;
; INPUTS:
;   image      - Image [NX,NY]
;   innvar     - Inverse variance image corresponding to IMAGE [NX,NY]
;   xcen       - Trace centers for traces running vertical [NY,NTRACE]
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   scatimg    - Scattered light image [NX,NY]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This is a very naive scattered light fitting routine that fits to
;   the light levels between the bundles of 20 fibers.  The assumption is
;   that one should force the model to fit the light level between fibers
;   20 and 21, 40 and 41, etc.  For BOSS data, this light is actually
;   dominated by wings in the PSF and not large-angle scattering.
;   The Lorentzian wings from the VPH grating are in the vertical (wavelength)
;   direction, so do not show up in the inter-fiber-bundle region.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   splog
;
; REVISION HISTORY:
;   22-Dec-2011  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function fitscatter, image, invvar, xcen
   nx = (size(image,/dimens))[0]
   dims = size(xcen,/dimens)
   ny = dims[0]
   nfiber = dims[1]
   ycen = djs_laxisgen(dims, iaxis=0)
   nbundle = nfiber / 20
   xt = fltarr(ny,nbundle+1)
   yt = djs_laxisgen([ny,nbundle+1], iaxis=0)
   xt[*,0] = xcen[*,0] - 1.2 * (xcen[*,1] - xcen[*,0])
   for ib=1, nbundle-1 do $
    xt[*,ib] = 0.5 * (xcen[*,ib*20-1] + xcen[*,ib*20])
   xt[*,nbundle] = xcen[*,nfiber-1] + 1.2 * (xcen[*,nfiber-1] - xcen[*,nfiber-2])

   ; Nearest-pixel value (or use extract_asymbox2)
   imgt = image[round(xt),round(yt)]
   ivart = invvar[round(xt),round(yt)]

   ; Replace these values with a b-spline fit
   for ib=0, nbundle do begin
      sset1 = bspline_iterfit(yt[*,ib], imgt[*,ib], invvar=ivart[*,ib], $
       nord=3, maxiter=5, bkspace=30)
      imgt[*,ib] = bspline_valu(yt[*,ib], sset1)
   endfor

   ; Interpolate across a full image
   scatimg = fltarr(nx,ny)
   xvec = findgen(nx)
   for iy=0, ny-1 do $
    scatimg[*,iy] = interpol(reform(imgt[iy,*]), reform(xt[iy,*]), xvec)

   ; Use constant values beyond the used regions of the CCD on top + bottom
   scatimg = djs_maskinterp(scatimg, (invvar EQ 0), iaxis=1, /const)
   ; Use constant values beyond the used regions of the CCD on left + right
   for iy=0, ny-1 do $
    scatimg[*,iy] = djs_maskinterp(scatimg[*,iy], xvec LT xt[iy,0] OR $
     xvec GT xt[iy,nbundle], /const)

   return, scatimg
end
