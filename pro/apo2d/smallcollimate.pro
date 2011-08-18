;+
; NAME:
;   smallcollimate
;
; PURPOSE:
;   Compute the spectrograph collimation focus from Hartmann mask exposures.
;
; CALLING SEQUENCE:
;   smallcollimate, expnum1, [ expnum2, docams=, indir=, nregx=, nregy=, $
;    maxshift=, /nocheck, /debug, /simple ]
;
; INPUTS:
;   expnum1    - First exposure number of raw sdR file.
;
; OPTIONAL KEYWORDS:
;   expnum2    - Second exposure number; default to EXPNUM1+1.
;   docams     - Cameras to analyze; default to ['b1','b2','r1','r2'].
;   indir      - Input directory for files; default to searching for
;                files in $BOSS_SPECTRO_DATA/*.  If $BOSS_SPECTRO_DATA is not set,
;                then it is assumed to be /data/spectro.
;   nregx      - Number of sub-regions in the X dimension; default to 8.
;   nregy      - Number of sub-regions in the Y dimension; default to 8.
;   maxshift   - Maximum pixel shift to search in both X and Y; default to 3.
;   nocheck    - If set, then assume that the 1st exposure is Hartmann-l,
;                and the 2nd exposure is Hartmann-r, rather than looking at
;                the OBSCOMM header keyword.  This correct keywords are
;                added by the SOP goSpecFocus command, but will not be if
;                you simply move the collimator and shutters yourself.
;                (Deprecated to never do this check, since BOSS has
;                no header keywords.)
;   debug      - flag to remove hardwired subregions and fits and instead
;                process all of the data in the array.  Much slower than the
;                nominal routine but informative for tests of focus variation
;                across full illuminated area
;   simple      - generate simpler output, for the benefit of external programs.
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
;   sdssproc
;   splog
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   collimate_obscomm()
;   collimate_bad_hdr()
;
; REVISION HISTORY:
;   28-Mar-2002  Written by D. Schlegel, Princeton
;   01-Sep-2009  K. Dawson, Utah
;-
;------------------------------------------------------------------------------
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
; Return 1 if the header indicates there are no lamps on or no flat-field
; petals closed.  If the wrong number, then just print a warning.
function collimate_bad_hdr, hdr

   qbad = 0

   ffs = sxpar(hdr, 'FFS')
   if (keyword_set(ffs)) then begin
      ffs_sum = fix( total( fix( str_sep(ffs,' ') ) ) )
      if (ffs_sum LT 8) then $
       splog, 'WARNING: Only ' + strtrim(ffs_sum,2) + ' of 8 flat-field petals closed'
      if (ffs_sum EQ 0) then qbad = 1
   endif else begin
      splog, 'WARNING: FFS not in file header'
      qbad = 1
   endelse

   lamp_ne = sxpar(hdr, 'NE')
   lamp_hgcd = sxpar(hdr, 'HGCD')
   if (keyword_set(lamp_ne) AND keyword_set(lamp_hgcd)) then begin
      ne_sum = fix( total( fix( str_sep(lamp_ne,' ') ) ) )
      hgcd_sum = fix( total( fix( str_sep(lamp_hgcd,' ') ) ) )
      if (ne_sum LT 4) then $
       splog, 'WARNING: Only ' + strtrim(ne_sum,2) + ' of 4 Ne lamps turned on'
      if (hgcd_sum LT 4) then $
       splog, 'WARNING: Only ' + strtrim(ne_sum,2) + ' of 4 HgCd lamps turned on'
      if (ne_sum + hgcd_sum EQ 0) then qbad = 1
   endif else begin
      splog, 'WARNING: NE or HGCD not in file header'
      qbad = 1
   endelse

   return, qbad
end
;------------------------------------------------------------------------------
pro smallcollimate, expnum1, expnum2=expnum2, docams=docams, indir=indir, $
 nregx=nregx, nregy=nregy, maxshift=maxshift, nocheck=nocheck,debug=debug, $
 simple=simple
 

   if (n_params() LT 0) then begin
      print, 'Syntax - collimate, expnum1, [ expnum2, docams=, indir=, $'
      print, ' nregx=, nregy=, maxshift=, /nocheck, /debug, /simple ]'
      return
   endif

   ;----------
   ; Set defaults

   if (NOT keyword_set(expnum2)) then expnum2 = expnum1 + 1
   if (NOT keyword_set(docams)) then docams = ['b1','b2','r1','r2']
   if (NOT keyword_set(maxshift)) then maxshift = 2
   if (NOT keyword_set(nregx)) then nregx = 8
   if (NOT keyword_set(nregy)) then nregy = 8

   ngrow = 2 ; Grow bad pixels by this number of pixels in every dimension
   focustol = 0.20 ;if yoffset less than this, in focus

   ;----------
   ; If DOCAMS is an array, then call this routine recursively

   ncam = n_elements(docams)
   if (ncam GT 1) then begin
      for icam=0, ncam-1 do begin
         smallcollimate, expnum1, expnum2=expnum2, docams=docams[icam], indir=indir, $
          nregx=nregx, nregy=nregy, maxshift=maxshift, nocheck=nocheck,debug=debug,simple=simple
      endfor
      return
   endif

   ;----------
   ; Locate the input files (either compressed or un-compressed)

   if (NOT keyword_set(indir)) then begin
      indir = getenv('BOSS_SPECTRO_DATA')
      if (NOT keyword_set(indir)) then $
       indir = '/data/spectro'
      indir = indir + '/*'
   endif
   
   filename1pat = djs_filepath('sdR-' + docams[0] + '-' + string(expnum1, format='(i8.8)') $
                               + '.fit*', root_dir=indir)
   filename2pat = djs_filepath('sdR-' + docams[0] + '-' + string(expnum2, format='(i8.8)') $
                               + '.fit*', root_dir=indir)

   filename1 = (findfile(filename1pat, count=ct1))[0]
   filename2 = (findfile(filename2pat, count=ct2))[0]

   if (ct1 EQ 0 OR ct2 EQ 0) then begin
      print, 'All files not found: ' + filename1pat + ", " + filename2pat
      return
   endif

   ;----------
   ; Read the two images

   sdssproc, filename1, bigimg1, ivar1, hdr=hdr1, /silent
   sdssproc, filename2, bigimg2, ivar2, hdr=hdr2, /silent

   if (collimate_bad_hdr(hdr1) OR collimate_bad_hdr(hdr2)) then begin
      return
   endif
   
   ;region to look at
   yl=855
   yh=1300
   xl=1500
   xh=2500

   im=bigimg1[xl:xh,yl:yh]
   djs_iterstat,im,sigma=sigma,mean=meanim
   var=sigma*sigma
   maxim=max(im)
   diff=maxim-meanim  ;condition that there is light
   ; Remember original image size for contour plots
   dims_orig = size(bigimg1, /dimens)
   nx_orig = dims_orig[0]
   ny_orig = dims_orig[1]
   camname = strtrim(sxpar(hdr1, 'CAMERAS'),2)
   camcolor = strmid(camname,0,1)

;first check that the camera is capturing light by requiring
;that a minimum flux falls in the first three apertures
;flux in first three apertures must be greater than 1000 counts and
;also must be monotonically increasing

   cam_flag=0

   if (var lt 100 and diff lt 3000) then begin
      if keyword_set(simple) then begin
         print,"THERE DOES NOT APPEAR TO BE ANY LIGHT FROM THE ARCS IN ", camname, "!!!"
         return
      endif else begin
         print,""
         print,""
         print,""
         print,"ERROR FOUND!!"
         print,""
         print,""
         wait,5
         print,"THERE DOES NOT APPEAR TO BE ANY LIGHT FROM THE ARCS IN THIS CAMERA!!!"
         print,""
         print,""
         print,""
         wait,5
         cam_flag=1
      endelse
   endif

   dims = size(bigimg1, /dimens)
   nx = dims[0]
   ny = dims[1]


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

   moffset = fltarr(nregx,nregy)
   yoffset = fltarr(nregx,nregy)
   xcenter = fltarr(nregx,nregy)
   ycenter = fltarr(nregx,nregy)

   ix1=xl
   ix2=xh
   iy1=yl
   iy2=yh
   
   submask = gmask[ix1:ix2,iy1:iy2]
   submask[0:maxshift,*] = 0
   submask[*,0:maxshift] = 0
   submask[ix2-ix1-maxshift:ix2-ix1,*] = 0
   submask[*,iy2-iy1-maxshift:iy2-iy1] = 0
   subimg1 = bigimg1[ix1:ix2,iy1:iy2]
   subimg2 = bigimg2[ix1:ix2,iy1:iy2]
   xcenter=0.5*(ix1+ix2)
   ycenter=0.5*(iy1+iy2)


         ;----------
         ; Compute the linear correlation coefficients
;maxshift=3
   dy = 0.05
   nshift = ceil(2*maxshift/dy)
   yshift = -maxshift + dy * lindgen(nshift)
;
   ans = fltarr(nshift)
   for j=0, nshift-1 do $
      ans[j] = total( subimg1 * sshift2d(subimg2,[0,yshift[j]]) * submask )
   junk = max(ans, ibest)
   yoffset = yshift[ibest]
;         ; Good solution iff >75% pixels are good in the subregion
   moffset = mean(submask)      ;GT 0.75




   ;----------
   ; We have assumed that the sequence is Hartmann-l, then Hartmann-r.
   ; If in the other order, then reverse the results

   if (hartpos1 GT hartpos2) then begin
      yoffset = -yoffset
   endif

;   Determine region for optimal focus for R and B on BOSS spectrograph
;  for BOSS, choose region of red CCD that lies in the ~6000-8500 /AA range to avoid bias from the far red absorption depth
;   deep_red=2000
;   red_x=[min(xcenter),max(xcenter)]
;   red_y=[min(ycenter),ny/2-1]
;  for BOSS, choose region of blue CCD near center of image, excluding extremes of saddle that sit near edges
;   blue_x=[700,2700]            ;[1996,2196];[700,2700]
;   blue_y=[400,2000]            ;[1232,1285];[400,2000] ;how it was


;   if (nx_orig LT 4000) then igood = where(moffset)
;   if (nx_orig GE 4000 AND camcolor eq 'r') then igood = where(moffset AND ycenter gt red_y[0] AND ycenter lt red_y[1] AND xcenter gt red_x[0] AND xcenter lt red_x[1])
;   if (nx_orig GE 4000 AND camcolor eq 'b') then igood = where(moffset AND ycenter ge blue_y[0] AND ycenter le blue_y[1] AND xcenter ge blue_x[0] AND xcenter le blue_x[1])

;only assign offsets if camera is capturing light
   if cam_flag eq 1 then begin
      igood=-1
      meanyoffstr = 'N/A'
      meanydevstr = 'N/A'
      minyoffstr = 'N/A'
      maxyoffstr = 'N/A'
      yoffset[*,*]=0.0
      meanyoff = 0.0
      meanydev = 0.0
      minyoff = 0.0
      maxyoff = 0.0
   endif

   if cam_flag eq 0 then begin
;offsets for regions picked

      if docams eq 'b1' then begin
         m=1.01252
         b=.173888
      endif

      if docams eq 'b2' then begin
         m=1.04467
         b=-0.100018
      endif

      if docams eq 'r1' then begin
         m=1.00308
         b=0.0833203
      endif

      if docams eq 'r2' then begin
         m=0.999042
         b=0.0178034
      endif

;m=1   ;for testing uncomment these
;b=0

      yoffset=yoffset*m + b

      igood=0

      meanyoff = mean(yoffset[igood])
      meanydev = 0;stddev(yoffset[igood])
      minyoff = min(yoffset[igood])
      maxyoff = max(yoffset[igood])
      meanyoffstr = string(meanyoff, format='(f5.2)')
      meanydevstr = string(meanydev, format='(f5.2)')
      minyoffstr = string(minyoff, format='(f5.2)')
      maxyoffstr = string(maxyoff, format='(f5.2)')

   endif

   ;----------
   ; Write info to the log file, and echo only some to the terminal

   if not keyword_set(simple) then begin
      splog, file=logfile
      splog, 'FILE1 = ' + filename1
      splog, 'FILE2 = ' + filename2
      splog, 'MJD = ', mjdstr
      splog, camname , ' mean offset = ' + meanyoffstr + ' pix'
   end

;splog, 'Camera = ', camname
   ;splog, ' '
;splog, 'IMAGE OF WAVELENGTH-OFFSETS', /no_stdout
;splog, '(in the correct orientation that 0,0 is lower left)', /no_stdout
;splog, ' ', /no_stdout
  ; for iregy=nregy-1, 0, -1 do begin
   ;   format = '(' + string(nregx) + 'f6.2)'
    ;  splog, yoffset[*,iregy], format=format, /no_stdout
   ;endfor
;splog, ' ', /no_stdout
;splog, 'Min offset = ' + minyoffstr + ' pix', /no_stdout
;splog, 'Max offset = ' + maxyoffstr + ' pix', /no_stdout
;splog, 'RMS across CCD = ' + meanydevstr + ' pix', /no_stdout
;splog, ' ', /no_stdout

;splog, ' '

   ;----------
   ; Output the predicted movements in order to collimate

   spectroid = strmid(camname,1,1)

   tolstr = strtrim(string(focustol, format='(f6.2)'),2)
   if keyword_set(simple) then begin
      diags = ['"Out of focus"', '"In focus"']
      print, string(camname, ' MeanOffset ', strtrim(meanyoff,2), ",", $
                    diags[(abs(meanyoff) LT focustol)], $
                    format='(a2,a,f5.2,a,a)')
   endif else begin
      if (abs(meanyoff) LT focustol) then begin
         splog, 'Camera ' + camname $
                + ' appears to be IN-FOCUS (|' + strtrim(meanyoff,2) + '| < ' + tolstr + ' pix)'
      endif else begin
         splog, 'Camera ' + camname $
                + ' appears to be OUT-OF-FOCUS (|' + strtrim(meanyoff,2) + '| > ' + tolstr + ' pix)'
      endelse
   endelse


;ADD FUNNY FUDGE FACTORS - KD
   if (camcolor EQ 'r') then begin
      pixscale = -15.           ; microns
      if (camname EQ 'r1') then val = long( -9150. * meanyoff*1.07 * pixscale/24. )
      if (camname EQ 'r2') then val = long( -9150. * meanyoff*1.21 * pixscale/24. )

      if keyword_set(simple) then begin
         print, string(camname, " PistonMove ", strtrim(val,2))
      endif else begin
         splog, 'Predict ' + camname + ' piston movement of ' + strtrim(val,2) + ' steps.'
         splog, 'Issue command: boss moveColl ' $
                + 'piston=' + strtrim(val,2) $
                + 'spec=sp' + spectroid

         colla=sxpar(hdr1,"COLLA")
         collb=sxpar(hdr1,"COLLB")
         collc=sxpar(hdr1,"COLLC")
         splog, 'current   motorPosition ', colla, collb, collc
         splog, 'new  expected    motorPosition',  colla+val, collb+val, collc+val
      endelse
   endif else if (camcolor EQ 'b') then begin
      pixscale = 15.                        ; microns
      val = -35.69 * meanyoff * pixscale/24.

      if keyword_set(simple) then begin
         print, string(camname, " RingMove ", strtrim(val,2), $
                       format='(a,a,f5.1)')
      endif else begin
         splog, 'Predict ' + camname + ' camera ring movement of ' $
                + string(val, format='(f6.1)') + ' degrees. (if red is already in focus)'
      endelse
   endif

   if not keyword_set(simple) then begin
      splog, /close
   end
end
;------------------------------------------------------------------------------
