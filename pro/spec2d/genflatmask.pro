;+
; NAME:
;   genflatmask
;
; PURPOSE:
;   Read or generate a mask for pixels to ignore when generating flat-fields.
;
; CALLING SEQUENCE:
;   maskimg = genflatmask( [ flatname, spectrographid=, color=, indir=, ] )
;
; OPTIONAL INPUTS:
;   flatname   - Name of flat-field image
;                Note that many flats from many nights can be combined.
;
; OPTIONAL KEYWORDS:
;   spectrographid - Spectrograph ID (1 or 2); reqired if FLATNAME not set.
;   color      - Spectrograph color ('red' or 'blue');
;                reqired if FLATNAME not set.
;   indir      - Input directory for FLATNAME; default to '.'
;
; OUTPUTS:
;   maskimg    - Mask image with 0=good, 1=bad
;
; COMMENTS:
;   If FLATNAME is specified, then a bad pixel mask is generated from
;   that flat-field image assuming it has a mean value approximately equal
;   to one.  If not, then read the mask from the distribution of IDLSPEC2D
;   in the "etc" subdirectory.
;
;   Mask values =0 for good pixels, =1 for bad.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   mrdfits()
;   readfits()
;
; DATA FILES:
;   $IDLSPEC2D_DIR/etc/flatmask-b1.fits.gz
;   $IDLSPEC2D_DIR/etc/flatmask-b2.fits.gz
;   $IDLSPEC2D_DIR/etc/flatmask-r1.fits.gz
;   $IDLSPEC2D_DIR/etc/flatmask-r2.fits.gz
;
; REVISION HISTORY:
;   16-Dec-1999  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
function genflatmask, flatname, spectrographid=spectrographid, color=color, $
 indir=indir

   if (keyword_set(flatname)) then begin

      if (NOT keyword_set(indir)) then indir = '.'

      if (keyword_set(indir)) then inpath = filepath(flatname, root_dir=indir) $
       else inpath = flatname
      fullname = (findfile(inpath, count=ct))[0]
      if (ct NE 1) then $
       message, 'Cannot find image ' + flatname

      flatimg = readfits(flatname)
      dims = size(flatimg,/dimens)
      nx = dims[0]
      ny = dims[1]

      ; First mask all points less than 0.5 and everything within NGROW pixels

      ngrow = 2
      maskimg1 = flatimg LT 0.5
      ; Remove edge effects...
      maskimg1[0,*] = 0
      maskimg1[nx-1,*] = 0
      maskimg1[*,0] = 0
      maskimg1[*,ny-1] = 0
      maskimg1 = smooth(maskimg1 * (2*ngrow+1)^2, ngrow) GT 0

      ; Now find large regions of points less than 0.9,
      ; and grow them by many pixels
      maskimg2 = flatimg LT 0.9
      ; Remove edge effects...
      maskimg2[0,*] = 0
      maskimg2[nx-1,*] = 0
      maskimg2[*,0] = 0
      maskimg2[*,ny-1] = 0
      maskimg2 = smooth((smooth(maskimg2+0.0,7,/edge) GT 0.5) +0., 35) GT 0.15

      maskimg = maskimg1 OR maskimg2

   endif else begin

      if (NOT keyword_set(spectrographid) OR NOT keyword_set(color)) then $
       message, 'Must specify SPECTROGRAPHID and COLOR if FLATNAME not set'

      ; Generate the name of the mask file in the IDLSPEC2D distribution,
      ; and the name of the temporary gunzip-ed version.

      rootname = 'flatmask-' + string(format='(a1,i1)',color,spectrographid) $
       + '.fits'
      masknamein = filepath(rootname+'.gz', $
       root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='etc')
;      masknameout = filepath(rootname, root_dir=tmpdir)

      ; Make sure that the mask file exists

      junk = findfile(masknamein, count=ct)
      if (ct NE 1) then $
       message, 'Cannot find mask ' + masknamein

      ; Gunzip the mask file, read it, and dispose of the file

;      spawn, 'zcat ' + masknamein + ' > ' + masknameout
;      maskimg = readfits(masknameout)
;      rmfile, masknameout
      maskimg = mrdfits(masknamein)

   endelse

   return, maskimg
end
;------------------------------------------------------------------------------
