;+
; NAME:
;   boss_gain
;
; PURPOSE:
;   Compute the gain in all amplifiers from a pair of flat-field exposures.
;
; CALLING SEQUENCE:
;   boss_gain, expnum1, [ expnum2, docams=, indir=, $
;    gain_range=, thresh=, /doplot ]
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
;   gain_range - Range of gain values to test; default to [0.5,2.0,0.002]
;                to test the range [0.5,2.0] spaced every 0.002
;   thresh     - Threshold of pixel values to use; if a scalar, then this is
;                a lower limit only; if this is a 2-element vector, then it
;                is a lower and upper bound; default to [500,10000]
;   doplot     - If set, then generate a plot with the chi distributions
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   A number of gain values are explicitly tested, with the correct
;   gain assumed to be that which gives a Gaussian normal distribution.
;   That test is simply that if we difference the two flat-field images,
;   68.3% of the points fall within one standard deviation.
;
;   Only use pixels where the value in the first image is above THRESH counts.
;
;   A multiplicative scaling between each flat-field image is computed
;   in each row for each amplifier.  This is to account for the lamp
;   temperature or brightness changing.
;
;   If these are flats:
;   An X shift of the image (flexure) is also measured and reported
;   for each amplifier region.  The reported gains shouldn't be trusted
;   if the flexure is any more than a tiny fraction of a pixel.
;   Flexure is only tested in the range [-0.1,0.1] pix.
;
; EXAMPLES:
;   Solve for the gain of all 4 CCD's from exposures 101702+101703
;   on MJD 55093 (assuming the files exist in /data/spectro/55093):
;     IDL> boss_gain, 101702
;
; BUGS:
;   This routine will also attempt to compute the gain from warm darks,
;   but that appears to give very different answers.  From MJD 55106:
;     IDL> boss_gain, 102134, docams=['b2','r2']
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_filepath()
;   djs_iterstat
;   djs_oplot
;   djs_plot
;   hogg_iter_linfit
;   splog
;   mrdfits
;   sshift2d()
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   boss_gain_flexure()
;
; REVISION HISTORY:
;   10-Oct-2009  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
; Check that the flexure between these two images is small
function boss_gain_flexure, image1, image2

   best_corr = -1
   best_shift = 0
   for xshift=-0.100, 0.100, 0.005 do begin
      corr = total(image1 * sshift2d(image2, [xshift,0]))
      if (corr GT best_corr) then begin
         best_corr = corr
         best_xshift = xshift
      endif
   endfor

   return, best_xshift
end
;------------------------------------------------------------------------------
pro boss_gain, expnum1, expnum2, docams=docams, indir=indir, $
 gain_range=gain_range, thresh=thresh, doplot=doplot

   ;----------
   ; Set defaults

   if (NOT keyword_set(expnum2)) then expnum2 = expnum1 + 1
   if (NOT keyword_set(docams)) then docams = ['b1','b2','r1','r2']
   if (NOT keyword_set(gain_range)) then gain_range = [0.5,3.0,0.002]
   if (n_elements(thresh) EQ 0) then thresh = [500.,10000.]

   quiet = !quiet
   !quiet = 1

   ;----------
   ; If DOCAMS is an array, then call this routine recursively

   ncam = n_elements(docams)
   if (ncam GT 1) then begin
      for icam=0, ncam-1 do begin
         boss_gain, expnum1, expnum2, docams=docams[icam], indir=indir, $
          gain_range=gain_range, thresh=thresh, doplot=doplot
      endfor
      !quiet = quiet
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

   splog, 'Image #1 = ', filename1
   splog, 'Image #2 = ', filename2
   splog, 'Threshold = ', thresh[0], thresh[1], ' ADU'

   im1 = mrdfits(filename1, 0, hdr1, /fscale, /silent)
   im2 = mrdfits(filename2, 0, hdr2, /fscale, /silent)
   if (sxpar(hdr1,'MJD') LT 55113) then im1 = rotate(im1,2)
   if (sxpar(hdr2,'MJD') LT 55113) then im2 = rotate(im2,2)
   mjdstr = strtrim(sxpar(hdr1,'MJD'),2)

   flavor = strtrim(sxpar(hdr1, 'FLAVOR'),2)
   camcolor = strmid(docams[0],0,1)

   case camcolor of
   'b': begin
      xbias0 = [  0, 4284,   0,  4284]
      xbias1 = [ 67, 4340,  67,  4340]
      x0 = [ 129, 2176,  129,  2176]
      x1 = [2175, 4222, 2175,  4222]
      y0 = [  56,   56, 2112,  2112] 
      y1 = [2111, 2111, 4167,  4167]
      end
   'r': begin
      xbias0 = [  0, 4241,  0,  4241]
      xbias1 = [110, 4351,110,  4351]
      x0 = [ 120, 2176,  120,  2176]
      x1 = [2175, 4231, 2175,  4231]
      y0 = [ 100,  100, 2112,  2112] 
      y1 = [2111, 2111, 4123,  4123]
      end
   endcase

; Replace the first region with a simulated image...
;nn = 20. ; read noise in e-
;gg = 1.3 ; gain
;im1[0:110,100:2111] = nn * randomn(1111,111L*2012L) / gg
;im2[0:110,100:2111] = nn * randomn(2222,111L*2012L) / gg
;im1[120:2175,100:2111] = (300. + sqrt(300.+nn^2)*randomn(1234,2056L*2012L))/gg
;im2[120:2175,100:2111] = (800. + sqrt(800.+nn^2)*randomn(5678,2056L*2012L))/gg

   ndegree = 4 ; polynomial fit to 'dark' exposures
   minpts = 1000 ; At least minpts pixels

   if (keyword_set(doplot)) then begin
      !p.multi = [0,2,2]
      plotfile = 'Gain-'+(str_sep(fileandpath(filename1),'.'))[0]+'.ps'
      dfpsplot, plotfile, /color, /square
   endif

   namp = n_elements(x0)
   if (flavor EQ 'flat') then xshift = fltarr(namp)
   for iamp=0, namp-1 do begin
      nx = x1[iamp] - x0[iamp] + 1
      ny = y1[iamp] - y0[iamp] + 1
      xvec = lindgen(nx)
      numer = fltarr(nx,ny)
      denom1 = fltarr(nx,ny)
      denom2 = fltarr(nx,ny)

      subim1 = im1[x0[iamp]:x1[iamp],y0[iamp]:y1[iamp]]
      subim2 = im2[x0[iamp]:x1[iamp],y0[iamp]:y1[iamp]]

      if (n_elements(thresh) EQ 2) then $
       mask = subim1 GE thresh[0] AND subim1 LE thresh[1] $
      else $
       mask = subim1 GE thresh[0]

      if (flavor EQ 'flat') then $
       xshift[iamp] = boss_gain_flexure(subim1, subim2)

      for iy=0L, ny-1L do begin
         ; Read noise computed on a row-by-row basis in the overscan
         djs_iterstat, im1[xbias0[iamp]:xbias1[iamp],iy+y0[iamp]], $
          mean=bias1, sigma=noise_dn1
         djs_iterstat, im2[xbias0[iamp]:xbias1[iamp],iy+y0[iamp]], $
          mean=bias2, sigma=noise_dn2
         noise_dn1 *= 1.015
         noise_dn2 *= 1.015
         vec1 = im1[x0[iamp]:x1[iamp],iy+y0[iamp]] - bias1
         vec2 = im2[x0[iamp]:x1[iamp],iy+y0[iamp]] - bias2

         if (flavor EQ 'dark') then begin
            splog, 'Fitting polynomial to each row, order=', ndegree
;            res1 = poly_fit(xvec, vec1, ndegree, yfit=yfit1)
;            res2 = poly_fit(xvec, vec2, ndegree, yfit=yfit2)

            ; Fit with sigma-rejection, but assuming a gain of 1
            parr = transpose(poly_array(nx, ndegree))
            hogg_iter_linfit, parr, vec1, 1./((vec1>0) + noise_dn1^2), $
             par1, nsigma=5
            yfit1 = transpose(parr ## par1)

            hogg_iter_linfit, parr, vec2, 1./((vec2>0) + noise_dn2^2), $
             par2, nsigma=5
            yfit2 = transpose(parr ## par2)
         endif else begin
            yfit1 = 1.
            thisthresh = thresh[0] > 1
            ii = where(vec1 GT thisthresh AND vec2 GT thisthresh, ct)
            if (ct GT 10) then begin
               djs_iterstat, vec2[ii]/vec1[ii], mean=yfit2
            endif else begin
               yfit2 = 1.
               mask[*,iy] = 0
            endelse
         endelse

         if (n_elements(yfit1) GT 1) then mask[*,iy] *= (yfit1 GT 0)
         if (n_elements(yfit2) GT 1) then mask[*,iy] *= (yfit2 GT 0)

         numer[*,iy] = vec1/yfit1 - vec2/yfit2
         denom1[*,iy] = (vec1>0)/yfit1^2 + (vec2>0)/yfit2^2
         denom2[*,iy] = noise_dn1^2/yfit1^2 + noise_dn2^2/yfit2^2
      endfor

      imask = where(mask, nmask)

      numer2 = numer^2
      gain_best = 0
      if (nmask GT minpts) then begin
         for gain=gain_range[0], gain_range[1], gain_range[2] do begin
            chi2im = numer2 / (denom1/gain + denom2)
            frac = mean(chi2im[imask] LT 1)
            if (frac GT 0.683) then gain_best = gain
         endfor
      endif

      splog, 'Camera = '+docams[0], ' Amp = ', strtrim(iamp,2), $
       '  Gain = ', string(gain_best,format='(f6.3)'), $
       ' from ', string(100*float(nmask)/n_elements(mask),format='(f6.2)'), $
       ' % of pix'

      if (keyword_set(doplot)) then begin
         chiim = numer / sqrt(denom1/gain_best + denom2)

         xbin = 0.01
         xrange = [-3,3]
         nbin = (xrange[1]-xrange[0])/xbin
         xx = xrange[0] + xbin * findgen(nbin)
         dist1 = fltarr(nbin)
         if (nmask GT minpts) then $
          populate_image, dist1, (chiim[imask]-xrange[0])/xbin
         if (iamp EQ 0) then $
          title=fileandpath(filename1)+' MJD='+mjdstr+' '+flavor $
          else title=''
         thiswin = ([2,1,0,3])[iamp]
         !p.multi = [thiswin,2,2]
         djs_plot, xx, dist1, psym=10, title=title, $
          xrange=[-3,3], xstyle=1, xtitle='Chi', ytitle='Npix', /noerase
         dist2 = nmask * xbin * exp(-xx^2/2) / sqrt(2*!pi)
         djs_oplot, xx, dist2, color='red'
         xyouts, -2.5, total([0.05,0.35]*!y.crange), $
          'Amp #'+strtrim(iamp,2), charsize=2
         xyouts, -2.5, total([0.05,0.25]*!y.crange), $
          'Gain ='+string(gain_best,format='(f6.3)'), charsize=2
      endif
   endfor

   if (flavor EQ 'flat') then $
    splog, 'Measured flexures = ', string(xshift,format='(4f7.3)'), ' pix'

   if (keyword_set(doplot)) then begin
      !p.multi = [0,1,1]
      dfpsclose
   endif

   print, ''
   !quiet = quiet

   return
end
;------------------------------------------------------------------------------
