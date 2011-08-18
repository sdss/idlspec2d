;+
; NAME:
;   ccdtilt
;
; PURPOSE:
;   Compute the CCD tilts from a series of Hartmann mask exposures.
;
; CALLING SEQUENCE:
;   ccdtilt, expnum1, piston=, [ expnum2, docams=, indir=, nregx=, nregy=, $
;    dy=, maxshift=, /nocheck, /doplot ]
;
; INPUTS:
;   expnum1    - Hartmann-left exposure numbers
;   piston     - Piston position of focus ring in microns; if known in degrees,
;                then multiply by 793./360
;
; OPTIONAL KEYWORDS:
;   expnum2    - Hartmnann-right exposure numbers; default to EXPNUM1+1.
;   docams     - Cameras to analyze; default to ['b1','b2','r1','r2'].
;   indir      - Input directory for files; default to searching for
;                files in $BOSS_SPECTRO_DATA/*.  If $BOSS_SPECTRO_DATA is not set,
;                then it is assumed to be /data/spectro.
;   nregx      - Number of sub-regions in the X dimension; default to 8.
;   nregy      - Number of sub-regions in the Y dimension; default to 8.
;   dy         - Sampling of shift when comparing Hartmann-L and Hartmann-R
;                images; default to 0.05 pix.
;   maxshift   - Maximum pixel shift to search in both X and Y; default to 3.
;   nocheck    - If set, then assume that the 1st exposure is Hartmann-l,
;                and the 2nd exposure is Hartmann-r, rather than looking at
;                the OBSCOMM header keyword.  This correct keywords are
;                added by the SOP goSpecFocus command, but will not be if
;                you simply move the collimator and shutters yourself.
;                (Deprecated to never do this check, since BOSS has
;                no header keywords.)
;   doplot     - If set, then compute the offsets in different regions
;                of the CCD and generate a contour plot.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The focus of the collimator is measured by comparing two Hartmann
;   exposures of arc lamps, and looking for shifts in the arc line positions.
;   The required piston of the collimator mirrors is assumed to be
;   linearly related to the offset in arc lines in the wavelength direction.
;
;   If /DOPLOT is set, then the arc-line shifts are computed in NREGX
;   by NREGY regions on each CCD, and a plot is made.
;
;   The following files are output for each camera:
;     Collimate-$MJD-$CAMERA-$EXPNUM1.log
;     Collimate-$MJD-$CAMERA-$EXPNUM1.ps
;
;   The position of the Hartmann shutters is read from the OBSCOMM header
;   keywords for SDSS-I.  It is expected to be '{focus, hartmann l}'
;   for one exposure and '{focus, hartmann r}' for the other (in either
;   order).  For BOSS, the HARTMANN keyword is read and must be either
;   'Left' or 'Right'.
;   It is assumed that the collimator position is identical for both exposures.
;
;   The sense of the pixel shifts reported is what one would have to shift
;   the Hartmann-r exposure by in Y to agree with the Hartmann-l exposure.
;
; EXAMPLES:
;   Solve for the tilt of the r1 CCD from Hartmann data on MJD 55218:
;    expnum1 = [107910,107912,107908,107916,107918]
;    piston = [-50,-25,0,25,50] * 793./360.
;    ccdtilt, expnum1, piston=piston, docams='r1'
;
; BUGS:
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_filepath()
;   djs_icolor()
;   sdssproc
;   splog
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   collimate_obscomm()
;   collimate_xcorr()
;
; REVISION HISTORY:
;   03-Feb-2010  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
; Note that OBSCOMM should be '{focus, hartmann l}'
; or '{focus, hartmann r}'
function collimate_obscomm, hdr

   obscomm = strtrim(sxpar(hdr, 'OBSCOMM'),2)
   hartmann = strtrim(sxpar(hdr, 'HARTMANN'),2)
   if (obscomm EQ '{focus, hartmann l}' OR hartmann EQ 'Left') then $
    retval = -1 $
   else if (obscomm EQ '{focus, hartmann r}' OR hartmann EQ 'Right') then $
    retval = 1 $
   else retval = 0

   return, retval
end
;------------------------------------------------------------------------------
function collimate_xcorr, img1, img2, mask, ngrow=ngrow, $
 dy=dy, maxshift=maxshift, qgood=qgood

   dims = size(img1, /dimens)
   nx = dims[0]
   ny = dims[1]

   ;----------
   ; Make a mask of pixels that are good in both images, and grow
   ; any bad pixels by NGROW in each dimension.
   ; Also mask any pixels within MAXSHIFT from the edge of the CCD.

   width = ngrow*2 + 1
   newmask = smooth(mask*width^2, width, /edge_truncate) EQ width^2
   newmask[0:maxshift,*] = 0
   newmask[nx-1-maxshift:nx-1,*] = 0
   newmask[*,0:maxshift] = 0
   newmask[*,ny-1-maxshift:ny-1] = 0

   newimg1 = img1 * newmask
   newimg2 = img2 * newmask

   nshift = ceil(2*maxshift/dy)
   yshift = -maxshift + dy * lindgen(nshift)

   ans = fltarr(nshift)
   for j=0, nshift-1 do $
    ans[j] = total( newimg1 * sshift2d(newimg2,[0,yshift[j]]) * newmask )
   junk = max(ans, ibest)
   yoffset = yshift[ibest]

   ; Good solution iff >75% pixels are good in the subregion
   qgood = mean(newmask) GT 0.75

   return, yoffset
end
;------------------------------------------------------------------------------
pro ccdtilt1, expnum1, expnum2, docams=docams, indir=indir, $
 nregx=nregx, nregy=nregy, dy=dy, maxshift=maxshift, $
 nocheck=nocheck, doplot=doplot, $
 yoffset=yoffset, moffset=moffset, xcenter=xcenter, ycenter=ycenter

   if (NOT keyword_set(maxshift)) then maxshift = 3
   if (NOT keyword_set(nregx)) then nregx = 8
   if (NOT keyword_set(nregy)) then nregy = 8
   if (NOT keyword_set(dy)) then dy = 0.05
   ngrow = 2 ; Grow bad pixels by this number of pixels in every dimension
   focustol = 0.25

   ;----------
   ; Locate the input files (either compressed or un-compressed)

   if (NOT keyword_set(indir)) then begin
      indir = getenv('BOSS_SPECTRO_DATA')
      if (NOT keyword_set(indir)) then $
       indir = '/data/spectro'
      indir = indir + '/*'
   endif

   print,''
   print,'ANALYZING FOCUS FOR CAMERA ', docams
   filename1 = 'sdR-' + docams + '-' + string(expnum1, format='(i8.8)') $
    + '.fit*'
   filename2 = 'sdR-' + docams + '-' + string(expnum2, format='(i8.8)') $
    + '.fit*'

   filename1 = (findfile(djs_filepath(filename1, root_dir=indir), count=ct1))[0]
   filename2 = (findfile(djs_filepath(filename2, root_dir=indir), count=ct2))[0]

   if (ct1 EQ 0 OR ct2 EQ 0) then begin
      print, 'Files not found'
      return
   endif

   ;----------
   ; Read the two images

   sdssproc, filename1, bigimg1, bigivar1, hdr=hdr1, /silent
   sdssproc, filename2, bigimg2, bigivar2, hdr=hdr2, /silent

   camname = strtrim(sxpar(hdr1, 'CAMERAS'),2)
   camcolor = strmid(camname,0,1)
   mjd = sxpar(hdr1, 'MJD')
   mjdstr = string(mjd,format='(i5.5)')
   expnum1 = sxpar(hdr1, 'EXPOSURE')
   expnum2 = sxpar(hdr2, 'EXPOSURE')
   expstr = string(expnum1, expnum2, format='(i8.8,"-",i8.8)')

   ; Remember original image size for contour plots
   dims_big = size(bigimg1, /dimens)
   nx_big = dims_big[0]
   ny_big = dims_big[1]
   qboss = nx_big GT 4000 ; True if BOSS image

   ; Explicitly set the subregion in which to compute focus
   if (camcolor EQ 'b' AND qboss) then begin
      xtrim = [250,3820]
      ytrim = [820,3600]
   endif else if (camcolor EQ 'r' AND qboss) then begin
      xtrim = [240,3840]
      ytrim = [500,3560]
   endif else begin
      xtrim = [0,2047]
      ytrim = [0,2047]
   endelse

   if (keyword_set(nocheck)) then begin
      hartpos1 = -1
      hartpos2 = 1
   endif else begin
      hartpos1 = collimate_obscomm(hdr1)
      hartpos2 = collimate_obscomm(hdr2)
      if (hartpos1 EQ 0 OR hartpos2 EQ 0) then begin
         print, 'FITS header do not indicate these are Hartmann exposures'
         return
      endif
      if (hartpos1 EQ hartpos2) then begin
         print,' FITS header indicate both exposures had same Hartmann position'
         return
      endif
   endelse

   plotfile = 'Collimate-' + mjdstr + '-' + camname $
    + string(expnum1, format='("-",i8.8)') + '.ps'
   logfile = 'Collimate-' + mjdstr + '-' + camname $
    + string(expnum1, format='("-",i8.8)') + '.log'
   title = 'Collimation for MJD=' + mjdstr + ' Camera=' + camname $
    + ' Exp=' + expstr

   ;----------
   ; Compute offsets across the image

   moffset = fltarr(nregx,nregy)
   yoffset = fltarr(nregx,nregy)
   xcenter = fltarr(nregx,nregy)
   ycenter = fltarr(nregx,nregy)

   for iregx=0, nregx-1 do begin
      for iregy=0, nregy-1 do begin
         print, 'Region ', iregx*nregy+iregy, nregx*nregy, string(13b), $
          format='(a,i,i,a,$)'
         ix1 = long(xtrim[0] + iregx * float(xtrim[1]-xtrim[0]) / nregx)
         ix2 = long(xtrim[0] + (iregx+1) * float(xtrim[1]-xtrim[0]) / nregx - 1)
         iy1 = long(ytrim[0] + iregy * float(ytrim[1]-ytrim[0]) / nregy)
         iy2 = long(ytrim[0] + (iregy+1) * float(ytrim[1]-ytrim[0]) / nregy - 1)
         xcenter[iregx,iregy]=0.5*(ix1+ix2)
         ycenter[iregx,iregy]=0.5*(iy1+iy2)
         yoffset[iregx,iregy] = collimate_xcorr( $
          bigimg1[ix1:ix2,iy1:iy2], bigimg2[ix1:ix2,iy1:iy2], $
          (bigivar1[ix1:ix2,iy1:iy2] NE 0) $
           AND (bigivar2[ix1:ix2,iy1:iy2] NE 0), $
          ngrow=ngrow, dy=dy, maxshift=maxshift, qgood=qgood1)
         moffset[iregx,iregy] = qgood1
      endfor
   endfor
   print

   ;----------
   ; We have assumed that the sequence is Hartmann-l, then Hartmann-r.
   ; If in the other order, then reverse the results

   if (hartpos1 GT hartpos2) then begin
      yoffset *= -1
   endif

   ;----------

   splog, /close

   ;----------
   ; Make the contour plot of focus if arcs look okay

   if (keyword_set(doplot)) then begin
      dfpsplot, plotfile, /color, /square

      ; Resample+smooth the image by a factor of 4 for the contour plot
      yoffimg = min_curve_surf(yoffset, nx=nregx*4, ny=nregy*4)
      xaxis = nx_big * (findgen(nregx*4) + 0.5) / (nregx*4)
      yaxis = ny_big * (findgen(nregy*4) + 0.5) / (nregy*4)

      nlevel = 10
      level0 = min(yoffset)
      level1 = max(yoffset)
      levels = level0 + findgen(nlevel+1) * (level1 - level0) / nlevel
      c_colors = (levels GE 0) * djs_icolor('blue') $
       + (levels LT 0) * djs_icolor('red')
      contour, yoffimg, xaxis, yaxis, /follow, levels=levels, $
       c_colors=c_colors, c_charsize=2, $
       xrange=[0,nx_big], yrange=[0,ny_big], /xstyle, /ystyle, $
       xtitle='X [pix]', ytitle='Y [pix]', title=title
      meanyoffstr = strtrim(mean(yoffimg),2)
      xyouts, 0.5*nx_big, 0.92*ny_big, align=0.5, $
       'Mean offset = ' + meanyoffstr + ' pix', charsize=1.5, charthick=2

      dfpsclose
   endif

   return
end
;------------------------------------------------------------------------------
function ccdtilt_zero, piston1, yoffset1

   ; Sort everything by piston position
   isort = sort(piston1)
   piston = piston1[isort]
   yoffset = yoffset1[isort]

   np = n_elements(piston)
   i1 = max(where(piston LE 0) > 0)
   i2 = min(where(piston GE 0) < (np-1))
   i2 = (i2 > (i1+1))
   if (i2 EQ np) then begin
      i1 = i1-1
      i2 = i2-1
   endif

   slope = (piston[i2] - piston[i1]) / (yoffset[i2] - yoffset[i1]) 
   pzero = piston[i1] - slope * yoffset[i1]

   return, pzero
end
;------------------------------------------------------------------------------
pro ccdtilt, expnum1, expnum2, piston=piston, docams=docams, _EXTRA=Extra

   if (n_params() LT 1) then begin
      print, 'Syntax - ccdtilt, expnum1, piston=, [ expnum2, docams=, indir=, $'
      print, ' nregx=, nregy=, maxshift=, /nocheck, /doplot ]'
      return
   endif

   if (n_elements(piston) NE n_elements(expnum1)) then $
    message, 'Elements in PISTON must agree with EXPNUM1'

   ;----------
   ; Set defaults

   if (NOT keyword_set(expnum2)) then expnum2 = expnum1 + 1
   if (NOT keyword_set(docams)) then docams = ['b1','b2','r1','r2']

   quiet = !quiet
   !quiet = 1

   ;----------
   ; If DOCAMS is an array, then call this routine recursively

   ncam = n_elements(docams)
   if (ncam GT 1) then begin
      for icam=0, ncam-1 do begin
         ccdtilt, expnum1, expnum2, piston=piston, docams=docams[icam], $
          _EXTRA=Extra
      endfor
      !quiet = quiet
      return
   endif

   ;----------
   ; Loop through all the exposures

   nfile = n_elements(expnum1)
   for ifile=0, nfile-1 do begin
      print, 'Working on file ', ifile+1,' of ',nfile
      ccdtilt1, expnum1[ifile], expnum2[ifile], docams=docams, $
       indir=indir, _EXTRA=Extra, $
       yoffset=yoffset1, moffset=moffset1, xcenter=xcenter1, ycenter=ycenter1
      if (ifile EQ 0) then begin
         yoffset = yoffset1
         moffset = moffset1
      endif else begin
         yoffset = [[[yoffset]], [[yoffset1]]]
         moffset = [[[moffset]], [[moffset1]]]
      endelse
   endfor

   ;----------
   ; At each CCD location, solve for the best piston position

   dims = (size(yoffset, /dimens))[0:1]
   pbest = fltarr(dims) - 999.
   for i=0, dims[0]-1 do begin
   for j=0, dims[1]-1 do begin
      indx = where(moffset[i,j,*] NE 0, ct)
      if (ct GT 1) then $
       pbest[i,j] = ccdtilt_zero(piston[indx], yoffset[i,j,indx])
   endfor
   endfor

   print, 'Focus offset map in microns ' + docams
   print, 'Positive values have focus closer to camera'
   print, 'Fiber number runs left-to-right at X=' $
    +string(minmax(xcenter1),format='(2i8)')
   print, 'Wavelength runs bottom to top at Y=' $
    +string(minmax(ycenter1),format='(2i8)')
   for i=dims[1]-1, 0, -1 do $
    print, pbest[*,i],format='('+string(dims[0])+'i6)'

   print, ''
   !quiet = quiet

   return
end
;------------------------------------------------------------------------------
