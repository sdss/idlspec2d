;+
; NAME:
;   combsmallcollimate
;
; PURPOSE:
;   Compute the spectrograph collimation focus from Hartmann mask exposures.
;
; CALLING SEQUENCE:
;   smallcollimate, expnum1, [ expnum2, docam1s=,docams2= indir=, nregx=, nregy=, $
;    maxshift=, /nocheck, /debug, /simple,spec= ,/test ]
;
; INPUTS:
;   expnum1    - First exposure number of raw sdR file.
;
; OPTIONAL KEYWORDS:
;   expnum2    - Second exposure number; default to EXPNUM1+1.
;   docams1     - Cameras to analyze; default to ['b1','r1'].
;   docams2     - Cameras to analyze; default to ['b2','r2'].
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
;   simple     - generate simpler output, for the benefit of external
;                programs. [SEE NOTES BELOW]
;
;   spec       - to do them by spectrograph, use 'sp1' or 'sp2' or
;                ['sp1','sp2'](default)
;
;   badres     - bad residual for blue ring.  Greater than this,
;                observers will need to move the ring.  Default to 6
;                degrees
;   test       - flag that prints out 1 line for comparison with smallcollimate
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
; /SIMPLE COMMENTS:
;   The /simple output is sometimes read by an external program which
;   expects or passes on the following output lines, given ${specname}
;   in ['sp1', 'sp1'] and ${camname} in ['b1','b2',r1',r2'].
;
;    We will likely turn this output into a .par file, but in the meanwhile:
;
;    ${specname} AverageMove <steps>
;      - the move to apply to the ${specname} collimator.
;        If this line does not appear, the collimation process fails.
;
;    ${specname} Residuals <steps>,<angle>,<string>
;      - the residual error after applying the AverageMove. 
;          . steps for red,
;          . angle for blue,
;          . if the string is "OK" the blue ring position after the
;            move is taken to be within spec. Otherwise collimation 
;            is treated as failed, and the string is passed on to the humans.
;         
;    ${camname} MeanOffset <pixelError>
;      - For all cameras, the measured pixel error. Not used, but recorded.
;
;    ${camname} PistonMove <steps>
;      - For the red cameras, the collimator move required to correct
;        MeanOffset. Not used, but recorded.
;
;    ${camname} RingMove <angle>
;      - For the blue cameras, the ring rotation required to correct
;        MeanOffset. Not used, but recorded.
;
;    Those outputs all generate obvious keywords for the APO archiver,
;    etc. Any other outputs on stdout and stderr are passed on as
;    warnings. 
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
;   fits_wait()
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
       splog, 'WARNING: Only ' + strtrim(hgcd_sum,2) + ' of 4 HgCd lamps turned on'
      if (ne_sum + hgcd_sum EQ 0) then qbad = 1
   endif else begin
      splog, 'WARNING: NE or HGCD not in file header'
      qbad = 1
   endelse

   return, qbad
end
;------------------------------------------------------------------------------
pro combsmallcollimate, expnum1, expnum2=expnum2, docams1=docams1,docams2=docams2, indir=indir, $
  maxshift=maxshift, nocheck=nocheck,debug=debug, $
 simple=simple,spec=spec,badres=badres,test=test
 

   if (n_params() LT 0) then begin
      print, 'Syntax - collimate, expnum1, [ expnum2, docams1=,docams2=,spec=, indir=, $'
      print, ' nregx=, nregy=, maxshift=, /nocheck, /debug, /simple ]'
      return
   endif

   ;----------
   ; Set defaults
   if (NOT keyword_set(spec)) then spec=['sp1','sp2']
   if (NOT keyword_set(expnum2)) then expnum2 = expnum1 + 1
   if (NOT keyword_set(maxshift)) then maxshift = 2
   if (NOT keyword_set(badres)) then badres=6       ;bad residual on blue ring
   ngrow = 2 ; Grow bad pixels by this number of pixels in every dimension
   focustol = 0.20 ;if yoffset less than this, in focus

   ;----------
   ; If SPEC is an array, then call this routine recursively for each spectrograph

   nspec=n_elements(spec)
   if (nspec gt 1) then begin
      for ispec=0, nspec-1  do begin
          combsmallcollimate,expnum1,expnum2=expnum2,docams1=docams1, $
            indir=indir,maxshift=maxshift, nocheck=nocheck,debug=debug, $ 
            simple=simple,spec=spec[ispec],docams2=docams2, $ 
            badres=badres,test=test
      endfor
      return
   endif

if (NOT keyword_set(docams1)) then docams1=['b1','r1']
if (NOT keyword_set(docams2)) then docams2=['b2','r2']

    ;----------
   ; Locate the input files (either compressed or un-compressed)   
   if (NOT keyword_set(indir)) then begin
      indir = getenv('BOSS_SPECTRO_DATA')
      if (NOT keyword_set(indir)) then $
         indir = '/data/spectro'
      indir = indir + '/*'
   endif
   
if (spec eq 'sp1') then docams=docams1 ;picking appropriate camera
if (spec eq 'sp2') then docams=docams2 ; for each spectrograph

ncam=n_elements(docams)
bad_camera_exposure=0     ;if an exposure does not work (no light, no arcs, not enough light,etc) will =1

for icam=0,ncam-1 do begin     ;do focus for each camera
   filename1pat = djs_filepath('sdR-' + docams[icam] + '-' + string(expnum1, format='(i8.8)') $
                               + '.fit*', root_dir=indir)
   filename2pat = djs_filepath('sdR-' + docams[icam] + '-' + string(expnum2, format='(i8.8)') $
                               + '.fit*', root_dir=indir)
   
   filename1 = (findfile(filename1pat, count=ct1))[0]
   filename2 = (findfile(filename2pat, count=ct2))[0]
   
   if (ct1 EQ 0 OR ct2 EQ 0) then begin
      print, 'All files not found: ' + filename1pat + ", " + filename2pat
      return
   endif
   
                                ;----------
                                ; Read the two images
   
  if (fits_wait(filename1, deltat=5, tmax=30) EQ 0) then begin
     print, 'Failure reading file ' + filename1
     return
  endif
  if (fits_wait(filename2, deltat=5, tmax=30) EQ 0) then begin
     print, 'Failure reading file ' + filename2
     return
  endif
  ; sdssproc, filename1, bigimg1, ivar1, hdr=hdr1, /silent
                                ; sdssproc, filename2, bigimg2, ivar2,
                                ; hdr=hdr2, /silent 
  ;Not using  sdssproc because of subregion arrays affect bias
   bigimg1=mrdfits(filename1,0,hdr1)
   bigimg2=mrdfits(filename2,0,hdr2)

      bias1 = [ median(bigimg1[10:100,950:1338]), median(bigimg1[4250:4340,950:1338]), 1, 1 ]
      bias2 = [ median(bigimg2[10:100,950:1338]), median(bigimg2[4250:4340,950:1338]), 1, 1 ]

      if docams[icam] eq 'b1' then gain = [1.048, 1.048, 1.018, 1.006]  
      if docams[icam] eq 'b2' then gain = [1.040, 0.994, 1.002, 1.010] 
      if docams[icam] eq 'r1' then gain = [1.966, 1.566, 1.542, 1.546]
      if docams[icam] eq 'r2' then gain = [1.598, 1.656, 1.582, 1.594]

      if (docams[icam] eq 'b1' or docams[icam] eq 'b2') then begin
         bigimg1p=fltarr(4096,4112)
         bigimg2p=fltarr(4096,4112)
         bigimg1p[0:2047,0:2055]=gain[0]*(bigimg1[128:2175,56:2111]-bias1[0])
         bigimg1p[2048:4095,0:2055]=gain[1]*(bigimg1[2176:4223,56:2111]-bias1[1])
         ;bigimg1p[0:2047,2056:4111]=gain[2]*(bigimg1[128:2175,2112:4167]-bias1[2])                       ;not needed since only use bottom portion of ccd
         ;bigimg1p[2048:4095,2056:4111]=gain[3]*(bigimg1[2176:4223,2112:4167]-bias1[3])
         bigimg2p[0:2047,0:2055]=gain[0]*(bigimg2[128:2175,56:2111]-bias2[0])
         bigimg2p[2048:4095,0:2055]=gain[1]*(bigimg2[2176:4223,56:2111]-bias2[1])
         ;bigimg2p[0:2047,2056:4111]=gain[2]*(bigimg2[128:2175,2112:4167]-bias2[2])
         ;bigimg2p[2048:4095,2056:4111]=gain[3]*(bigimg2[2176:4223,2112:4167]-bias2[3])
         bigimg1=fltarr(4096,4112)
         bigimg2=fltarr(4096,4112)
         bigimg1[*,*]=bigimg1p[0:4095,0:4111]
         bigimg2[*,*]=bigimg2p[0:4095,0:4111]
      endif else begin
         bigimg1p=fltarr(4114,4128)
         bigimg2p=fltarr(4114,4128)
         bigimg1p[0:2056,0:2063]=gain[0]*(bigimg1[119:2175,48:2111]-bias1[0])
         bigimg1p[2057:4113,0:2063]=gain[1]*(bigimg1[2176:4232,48:2111]-bias1[1])
         bigimg2p[0:2056,0:2063]=gain[2]*(bigimg2[119:2175,48:2111]-bias2[0])
         bigimg2p[2057:4113,0:2063]=gain[3]*(bigimg2[2176:4232,48:2111]-bias2[1])
         bigimg1=fltarr(4114,4128)
         bigimg2=fltarr(4114,4128)
         bigimg1[*,*]=bigimg1p[0:4113,0:4127]
         bigimg2[*,*]=bigimg2p[0:4113,0:4127]
      endelse
      bigimg1p=0
      bigimg2p=0
      bias1=0
      bias2=0
   if (collimate_bad_hdr(hdr1) OR collimate_bad_hdr(hdr2)) then begin
      return
   endif
   
                                ;region to look at
   yl=850
   yh=1300
   xl=1500
   xh=2500
 
   im=bigimg1[xl:xh,yl:yh]

   for i=0, 10 do begin         ;Apodization step, to minimize edge effects
       im[i,*]=i/10.*im[i,*]
       im[450-i,*]=i/10.*im[450-i,*]
   endfor

   if (docams[icam] eq 'b1' or docams[icam] eq 'b2') then var=variance(im[*,300:450])  else var=variance(im[*,0:150]) ; find variance near bright lines
  
         ; Remember original image size for contour plots
   dims_orig = size(bigimg1, /dimens)
   nx_orig = dims_orig[0]
   ny_orig = dims_orig[1]
   camname = strtrim(sxpar(hdr1, 'CAMERAS'),2)
   camcolor = strmid(camname,0,1)
   
;first check that the camera is capturing light by requiring
;that the variance is greater than 100
   
   cam_flag=0
   
   if (var lt 100 ) then begin 
      if keyword_set(simple) then begin
         print,"THERE DOES NOT APPEAR TO BE ANY LIGHT FROM THE ARCS IN ", camname, "!!!"
         cam_flag=1
         return
      endif else begin
         print,""
         print,""
         print,""
         print,"ERROR FOUND!!"
         print,""
         print,""
         wait,3
         print,"THERE DOES NOT APPEAR TO BE ANY LIGHT FROM THE ARCS IN THIS CAMERA!!!"
         print,""
         print,""
         print,""
         wait,3
         cam_flag=1
         return
      endelse
   endif

   dims = size(bigimg1, /dimens)
   nx = dims[0]
   ny = dims[1]


   if (keyword_set(nocheck)) then begin  ;if left right order of exposures switched
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
                 ; any bad pixels by NGROW in each dimension
                                ; Also mask any pixels within MAXSHIFT from the edge of the CCD.
;djs_iterstat,bigimg1,ivar=ivar1 ;can't use these since subarray
;djs_iterstat,bigimg2,ivar=ivar2

;;;;;;;;;;;;;;;;;;;
   width = ngrow*2 + 1
;   gmask = (ivar1 NE 0) AND (ivar2 NE 0)
   gmask = bigimg1*0.+1.
   gmask = smooth(gmask*width^2, width, /edge_truncate) EQ width^2
   gmask[0:maxshift,*] = 0
   gmask[nx-1-maxshift:nx-1,*] = 0
   gmask[*,0:maxshift] = 0
   gmask[*,ny-1-maxshift:ny-1] = 0
   
   bigimg1 = bigimg1 * gmask
   bigimg2 = bigimg2 * gmask
   
              ;use one subregion                 
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

       if docams[icam] eq 'b1' then begin  ;using collimate and not collimate,/debug with data from 56434; focus ring change
         m=1.00
         b=0.028
         yoffsetb1=yoffset*m + b
      endif

      if docams[icam] eq 'b2' then begin
         m=1.00
         b=-0.05
         yoffsetb2=yoffset*m + b
      endif

      if docams[icam] eq 'r1' then begin
         m=1.00
         b=-0.25
         yoffsetr1=yoffset*m + b
      endif

      if docams[icam] eq 'r2' then begin
         m=1.00
         b=-0.013
         yoffsetr2=yoffset*m + b
     endif

     if keyword_set(test) then b=0. ;for determining above offsets


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
      splog, 'Camera = ', camname
   splog, ' '
splog, 'IMAGE OF WAVELENGTH-OFFSETS', /no_stdout
splog, '(in the correct orientation that 0,0 is lower left)', /no_stdout
splog, ' ', /no_stdout
splog, ' ', /no_stdout
;splog, 'Min offset = ' + minyoffstr + ' pix', /no_stdout
;splog, 'Max offset = ' + maxyoffstr + ' pix', /no_stdout
;splog, 'RMS across CCD = ' + meanydevstr + ' pix', /no_stdout
splog, ' ', /no_stdout
splog, camname , ' mean offset = ' + meanyoffstr + ' pix'
if keyword_set(test) then print,'test '+string(expnum1,format='(i8.8)')+ ' ' +camname+ $
  ' mean offset = '+ meanyoffstr +' pix'
splog, ' '
   endif
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
      if (camname EQ 'r1') then val = long( -9150. * meanyoff*1.12* pixscale/24. )
      if (camname EQ 'r2') then val = long( -9150. * meanyoff*1.12* pixscale/24. ); 1.21*1.09 is new motor factor from Wednesday 11-10, looks like 1.12 from 55513
      valr=val
      if keyword_set(simple) then begin
         print, string(camname, " PistonMove ", strtrim(val,2))
      endif else begin
         splog, 'Predict ' + camname + ' piston movement of ' + strtrim(val,2) + ' steps.'
         splog, 'Issue command: boss moveColl' $
                + ' piston=' + strtrim(val,2) $
                + ' spec=sp' + spectroid

         colla=sxpar(hdr1,"COLLA")
         collb=sxpar(hdr1,"COLLB")
         collc=sxpar(hdr1,"COLLC")
         splog, 'current   motorPosition ', colla, collb, collc
         splog, 'new  expected    motorPosition',  colla+val, collb+val, collc+val
      endelse
   endif else if (camcolor EQ 'b') then begin
      pixscale = 15.                        ; microns
      if (camname EQ 'b1') then val = -31.87 * meanyoff * pixscale/24.;-2.35; 2.35 is offset from collimate to combsmallcollimate
      if (camname EQ 'b2') then val = -28.95 * meanyoff * pixscale/24.;+2.26
      if camname eq 'b1' then valb=-val*292;274;steps per degree
      if camname eq 'b2' then valb=-val*292;310
      if keyword_set(simple) then begin
         print, string(camname, " RingMove ", strtrim(val,2), $
                       format='(a,a,f5.1)')
      endif else begin
         splog, 'Predict ' + camname + ' camera ring movement of ' $
                + string(val, format='(f6.1)') + ' degrees. (if red is already in focus)'
      endelse
      val=valb
   endif
   if cam_flag eq 1 then begin
      bad_camera_exposure=1
      return
   endif
endfor                          ;for camera loop

bres=0
rres=0
;average move
if bad_camera_exposure eq 1 then ncam=1 ;don't count bad camera anymore
if ncam eq 2 then  begin 
   val=(valr+valb)/2.
   bres=-(valb-val)/292.             
   rres=valr-val
endif

if abs(rres) gt 1700 then begin  ;previously 1000; 1700 corresponds to ~.22 pix
   if NOT keyword_set(simple) then begin 
      splog,'    '
      splog,'      !!!!!!!!!!!!!!!!!!!! '
      splog, 'RESIDUALS OUT-OF-FOCUS (|' + strtrim(rres,2) + '| > 1700 pix)'
      splog,''
      splog,bres, ' residual for blue in degrees'
      splog,rres, ' residual for red in steps'
   endif
endif

if abs(bres) lt badres then begin
   resq='"OK"' 
endif else begin
   resq=string('"Bad angle !!!!!!! move blue ring ', bres*2, $
               ' degrees and then run gotoField noSlew"', $
               format='(a,f0.1,a)')
endelse

if keyword_set(simple) then begin
   print, string(spec, " Residuals ", strtrim(rres,2), ',', $
                 strtrim(bres,2),',', resq, $
                 format='(a,a,i0,a,f0.1,a,a)')
   print, string(spec, " AverageMove ", strtrim(val,2), format='(a,a,i)')
endif else begin 
   val=ceil(val)
   splog,' '
   splog,' To split the difference: boss moveColl piston='+ strtrim(val,2) $
                + ' spec=sp' + spectroid
   splog,' '
   splog,' The blue camera has a residual move of ',bres,' degrees' 
   splog,' The red camera has a residual  move of ',rres,' steps'
if abs(bres) gt badres then begin
   splog,'Bad angle !!!!!!!!!!.  Move blue ring ', bres*2, ' degrees and then run gotoField noSlew'
endif

   splog, /close
endelse
end

;------------------------------------------------------------------------------


