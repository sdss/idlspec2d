;+
; NAME:
;   collimate
;
; PURPOSE:
;   Compute the spectrograph collimation focus from Hartmann mask exposures.
;
; CALLING SEQUENCE:
;   collimate, expnum1, [ expnum2, docams=, indir=, nregx=, nregy=, $
;    maxshift=, /nocheck ]
;
; INPUTS:
;   expnum1    - First exposure number of raw sdR file.
;
; OPTIONAL KEYWORDS:
;   expnum2    - Second exposure number; default to EXPNUM1+1.
;   docams     - Cameras to analyze; default to ['b1','b2','r1','r2'].
;   indir      - Input directory for files; default to searching for
;                files in $RAWDATA_DIR/*.  If $RAWDATA_DIR is not set,
;                then it is assumed to be /data/spectro.
;   nregx      - Number of sub-regions in the X dimension; default to 8.
;   nregy      - Number of sub-regions in the Y dimension; default to 8.
;   maxshift   - Maximum pixel shift to search in both X and Y; default to 3.
;   nocheck    - If set, then assume that the 1st exposure is Hartmann-l,
;                and the 2nd exposure is Hartmann-r, rather than looking at
;                the OBSCOMM header keyword.  This correct keywords are
;                added by the SOP goSpecFocus command, but will not be if
;                you simply move the collimator and shutters yourself.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The focus of the collimator is measured by comparing two Hartmann
;   exposures of arc lamps, and looking for shifts in the arc line positions.
;   A linear correlation coefficient is computed independently in NREGX
;   by NREGY regions on each CCD as a function of pixel shifts of the 2nd
;   image in both X and Y.  The best-focus value is found in each region
;   by maximizing the linear correlation in the Y (wavelength) direction.
;   
;   The following files are output for each camera:
;     Collimate-$MJD-$CAMERA-$EXPNUM1.log
;     Collimate-$MJD-$CAMERA-$EXPNUM1.ps
;
;   The position of the Hartmann shutters is read from the OBSCOMM header
;   keywords.  It is expected to be '{focus, hartmann l}' for one exposure
;   and '{focus, hartmann r}' for the other (in either order).  It is assumed
;   that the collimator position is identical for both exposures.
;
;   The sense of the pixel shifts reported is what one would have to shift
;   the Hartmann-r exposure by in Y to agree with the Hartmann-l exposure.
;
; EXAMPLES:
;   Solve for the focus of all 4 CCD's from exposures 10812+10813
;   on MJD 52161 (assuming the files exist in /data/spectro/52161):
;     IDL> collimate, 10812
;
; BUGS:
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_filepath()
;   djs_icolor()
;   djs_svdfit()
;   sdssproc
;   splog
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   collimate_obscomm()
;
; REVISION HISTORY:
;   28-Mar-2002  Written by D. Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
; Note that OBSCOMM should be '{focus, hartmann l}'
; or '{focus, hartmann r}'
function collimate_obscomm, hdr

   obscomm = sxpar(hdr, 'OBSCOMM')
   if (obscomm EQ '{focus, hartmann l}') then retval = -1 $
    else if (obscomm EQ '{focus, hartmann r}') then retval = 1 $
    else retval = 0

   return, retval
end

;------------------------------------------------------------------------------
pro collimate, expnum1, expnum2, docams=docams, indir=indir, $
 nregx=nregx, nregy=nregy, maxshift=maxshift, nocheck=nocheck

   if (n_params() LT 1) then begin
      print, 'Syntax - collimate, expnum1, [ expnum2, docams=, indir=, $'
      print, ' nregx=, nregy=, maxshift=, /nocheck ]'
      return
   endif

   ;----------
   ; Set defaults

   if (NOT keyword_set(expnum2)) then expnum2 = expnum1 + 1
   if (NOT keyword_set(docams)) then docams = ['b1','b2','r1','r2']
   if (NOT keyword_set(maxshift)) then maxshift = 3
   if (NOT keyword_set(nregx)) then nregx = 8
   if (NOT keyword_set(nregy)) then nregy = 8

   ngrow = 2 ; Grow bad pixels by this number of pixels in every dimension
   focustol = 0.15

   quiet = !quiet
   !quiet = 1

   ;----------
   ; If DOCAMS is an array, then call this routine recursively

   ncam = n_elements(docams)
   if (ncam GT 1) then begin
      for icam=0, ncam-1 do begin
         collimate, expnum1, expnum2, docams=docams[icam], indir=indir, $
          nregx=nregx, nregy=nregy, maxshift=maxshift
      endfor
      !quiet = quiet
      return
   endif

   ;----------
   ; Locate the input files (either compressed or un-compressed)

   if (NOT keyword_set(indir)) then begin
      indir = getenv('RAWDATA_DIR')
      if (NOT keyword_set(indir)) then $
       indir = '/data/spectro'
      indir = indir + '/*'
   endif

   filename1 = 'sdR-' + docams[0] + '-' + string(expnum1, format='(i8.8)') $
    + '.fit*'
   filename2 = 'sdR-' + docams[0] + '-' + string(expnum2, format='(i8.8)') $
    + '.fit*'

   filename1 = (findfile(djs_filepath(filename1, root_dir=indir), count=ct1))[0]
   filename2 = (findfile(djs_filepath(filename2, root_dir=indir), count=ct2))[0]
   if (ct1 EQ 0 OR ct2 EQ 0) then begin
      print, 'Files not found'
      return
   endif

   ;----------
   ; Read the two images

   print
   print, 'Reading ' + filename1
   sdssproc, filename1, bigimg1, ivar1, hdr=hdr1, /silent
   print, 'Reading ' + filename2
   sdssproc, filename2, bigimg2, ivar2, hdr=hdr2, /silent
   dims = size(bigimg1, /dimens)
   nx = dims[0]
   ny = dims[1]
   if (n_elements(bigimg1) NE n_elements(bigimg2)) then begin
      print, 'Dimensions of the two images do not agree'
      return
   endif

   if (keyword_set(nocheck)) then begin
      hartpos1 = -1
      hartpos2 = 1
   endif else begin
      hartpos1 = collimate_obscomm(hdr1)
      hartpos2 = collimate_obscomm(hdr2)
      if (NOT keyword_set(hartpos1) OR NOT keyword_set(hartpos2)) then begin
         print, 'OBSCOMM keywords in header do not indicate these are Hartmann exposures'
         return
      endif
      if (hartpos1 EQ hartpos2) then begin
         print,' OBSCOMM keywords in header indicate both exposures had same Hartmann position'
         return
      endif
   endelse

   camname = strtrim(sxpar(hdr1, 'CAMERAS'),2)
   mjd = sxpar(hdr1, 'MJD')
   mjdstr = string(mjd,format='(i5.5)')
   expnum1 = sxpar(hdr1, 'EXPOSURE')
   expnum2 = sxpar(hdr2, 'EXPOSURE')
   expstr = string(expnum1, expnum2, format='(i8.8,"-",i8.8)')

   plotfile = 'Collimate-' + mjdstr + '-' + camname $
    + string(expnum1, format='("-",i8.8)') + '.ps'
   logfile = 'Collimate-' + mjdstr + '-' + camname $
    + string(expnum1, format='("-",i8.8)') + '.log'
   title = 'Collimation for MJD=' + mjdstr + ' Camera=' + camname $
    + ' Exp=' + expstr

   ;----------
   ; Make a mask of pixels that are good in both images, and grow
   ; any bad pixels by NGROW in each dimension.
   ; Also mask any pixels within MAXSHIFT from the edge of the CCD.

   width = ngrow*2 + 1
   gmask = (ivar1 NE 0) AND (ivar2 NE 0)
   gmask = smooth(gmask*width^2, width, /edge_truncate) EQ width^2
   gmask[0:maxshift,*] = 0
   gmask[nx-1-maxshift:nx-1,*] = 0
   gmask[*,0:maxshift] = 0
   gmask[*,ny-1-maxshift:ny-1] = 0

   bigimg1 = bigimg1 * gmask
   bigimg2 = bigimg2 * gmask

   ;----------
   ; Loop over each possible sub-image

   xoffset = fltarr(nregx,nregy)
   yoffset = fltarr(nregx,nregy)

   for iregx=0, nregx-1 do begin
      for iregy=0, nregy-1 do begin
         ix1 = long(iregx * float(nx) / nregx)
         ix2 = long((iregx+1) * float(nx) / nregx - 1)
         iy1 = long(iregy * float(ny) / nregy)
         iy2 = long((iregy+1) * float(ny) / nregy - 1)
         subimg1 = bigimg1[ix1:ix2,iy1:iy2]
         subimg2 = bigimg2[ix1:ix2,iy1:iy2]

         ;----------
         ; Compute the linear correlation coefficients

         nshift = 2*maxshift + 1
         xshift = lindgen(nshift) - maxshift
         yshift = lindgen(nshift) - maxshift
         ans = dblarr(nshift,nshift)
         for i=0, nshift-1 do $
          for j=0, nshift-1 do $
           ans[i,j] = total( subimg1 * shift(subimg2,xshift[i],yshift[j]) )
         ans = ans / max(ans)

         ;----------
         ; Find the peak in the correlation coefficient

         junk = max(ans, ibest)
         ibestx = ibest MOD nshift
         ibesty = ibest / nshift

         if (ibestx GT 0 AND ibestx LT nshift-1) then begin
;            coeff = svdfit(xshift[ibestx-1:ibestx+1], $
;             ans[ibestx-1:ibestx+1,ibesty], 3, /double)

            thisx = xshift[ibestx-1:ibestx+1]
            thisy = ans[ibestx-1:ibestx+1,ibesty]
            mmatrix = [[1d0,1d0,1d0], [thisx], [thisx^2]]
            coeff = invert(mmatrix, /double) # thisy

            xoffset[iregx,iregy] = -0.5 * coeff[1] / coeff[2]
         endif

         if (ibesty GT 0 AND ibesty LT nshift-1) then begin
;            coeff = svdfit(yshift[ibesty-1:ibesty+1], $
;             transpose(ans[ibestx,ibesty-1:ibesty+1]), 3, /double)

            thisx = yshift[ibesty-1:ibesty+1]
            thisy = transpose(ans[ibestx,ibesty-1:ibesty+1])
            mmatrix = [[1d0,1d0,1d0], [thisx], [thisx^2]]
            coeff = invert(mmatrix, /double) # thisy

            yoffset[iregx,iregy] = -0.5 * coeff[1] / coeff[2]
         endif
      endfor
   endfor

   ;----------
   ; We have assumed that the sequence is Hartmann-l, then Hartmann-r.
   ; If in the other order, then reverse the results

   if (hartpos1 GT hartpos2) then begin
      xoffset = -xoffset
      yoffset = -yoffset
   endif

   meanyoff = mean(yoffset)
   meanydev = stddev(yoffset)
   minyoff = min(yoffset)
   maxyoff = max(yoffset)

   meanyoffstr = string(meanyoff, format='(f5.2)')
   meanydevstr = string(meanydev, format='(f5.2)')
   minyoffstr = string(minyoff, format='(f5.2)')
   maxyoffstr = string(maxyoff, format='(f5.2)')

   ;----------
   ; Write info to the log file, and echo only some to the terminal

   splog, file=logfile
   splog, 'FILE1 = ' + filename1, /no_stdout
   splog, 'FILE2 = ' + filename2, /no_stdout
   splog, 'MJD = ', mjdstr, /no_stdout
   splog, 'Camera = ', camname, /no_stdout
   splog, ' ', /no_stdout
   splog, 'IMAGE OF WAVELENGTH-OFFSETS', /no_stdout
   splog, '(in the correct orientation that 0,0 is lower left)', /no_stdout
   splog, ' ', /no_stdout
   for iregy=nregy-1, 0, -1 do begin
      format = '(' + string(nregx) + 'f6.2)'
      splog, yoffset[*,iregy], format=format, /no_stdout
   endfor
   splog, ' ', /no_stdout
   splog, 'Min offset = ' + minyoffstr + ' pix'
   splog, 'Max offset = ' + maxyoffstr + ' pix'
   splog, 'RMS across CCD = ' + meanydevstr + ' pix'
   splog, ' ', /no_stdout
   splog, 'Mean offset = ' + meanyoffstr + ' pix'
   splog, ' ', /no_stdout

   ;----------
   ; Output the predicted movements in order to collimate

   camcolor = strmid(camname,0,1)
   spectroid = strmid(camname,1,1)

   tolstr = strtrim(string(focustol, format='(f6.2)'),2)
   if (abs(meanyoff) LT focustol) then begin
      splog, 'Camera ' + camname $
       + ' appears to be IN-FOCUS (|mean| < ' + tolstr + ' pix)'
   endif else begin
      splog, 'Camera ' + camname $
       + ' appears to be OUT-OF-FOCUS (|mean| > ' + tolstr + ' pix)'
   endelse

   if (camcolor EQ 'r') then begin
      val = long( -9150. * meanyoff )
      splog, 'Predict (red) piston movement of ', val, $
       ' steps for spectro-' + spectroid
      splog, '  Issue SOP command: "sp' + spectroid $
       + '; mechPiston ' + strtrim(string(val),2) + '"'
   endif else if (camcolor EQ 'b') then begin
      val = -35.69 * meanyoff
      splog, 'Predict blue camera ring movement of ' $
       + string(val, format='(f6.1)') + ' degrees for spectro-' + spectroid
      splog, '  (if red is already in focus)'
   endif

   splog, /close

   ;----------
   ; Make the contour plot of focus

   dfpsplot, plotfile, /color, /square

   ; Resample+smooth the image by a factor of 4 for the contour plot
   yoffimg = min_curve_surf(yoffset, nx=nregx*4, ny=nregy*4)
   xaxis = nx * (findgen(nregx*4) + 0.5) / (nregx*4)
   yaxis = ny * (findgen(nregy*4) + 0.5) / (nregy*4)

   lspace = 0.1
   level0 = floor(min(yoffset) / lspace) * lspace
   nlevel = ceil((max(yoffset) - level0) / lspace) + 1
   levels = level0 + lindgen(nlevel) * lspace
   c_colors = (levels GE 0) * djs_icolor('blue') $
    + (levels LT 0) * djs_icolor('red')
   contour, yoffimg, xaxis, yaxis, /follow, levels=levels, $
    c_colors=c_colors, c_charsize=2, $
    xrange=[0,nx], yrange=[0,ny], /xstyle, /ystyle, $
    xtitle='X [pix]', ytitle='Y [pix]', title=title
   xyouts, 0.5*nx, 0.96*ny, align=0.5, $
    'Mean offset = ' + meanyoffstr + ' pix', charsize=1.5

   dfpsclose

   print, 'Log file        = ' + logfile
   print, 'PostScript file = ' + plotfile

   !quiet = quiet
end
;------------------------------------------------------------------------------
