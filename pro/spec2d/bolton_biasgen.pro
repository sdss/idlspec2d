;+
;
; NAME:
;  bolton_biasgen
;
; PURPOSE:
;  Generate BOSS 4-amp master bias frames given an input file list.
;
; USAGE:
;  bolton_biasgen, bias_files [, outfile=, niter=, sigthresh=, outbias=]
;
; ARGUMENTS:
;   bias_files: list of full path qualified sdR files to use
;   outfile: output filename, default is, e.g., 'boss_pixbias-r1-55300.fits'
;     where the camera ID and MJD are taken from the header of the
;     first file in the list of inputs.
;   niter: number of modeling iterations, default is 10
;   sigthresh: sigma theshold for mask rejection, default is 5 (sigmas)
;
; OUTPUTS:
;   outbias (optional): the master bias image (also written to outfile)
;
; NOTES:
;   This routine builds a bias model given by
;     c_ijk = b_ij + d_k,
;   where c_ijk are the counts in pixel i,j of bias image k,
;   b_ij is the master bias, and d_k is the offset of image k.
;   The model is built independently in each of the four amp
;   quadrants, and the b_ij image is written out including overscan.
;
;   The model is iteratively optimized in a minimum chi^2 sense,
;   weighting by the estimated readnoise and doing model-based
;   masking of significantly discrepant pixels.
;
;   The output is intended to work with BOLTON_BIASSUB, which
;   determines d_k for arbitrary frames by analogous minimum chi^2
;   modeling with rejection in the overscan region for each amp.
;
; WRITTEN:
;  A. Bolton, U. of Utah, 2011aug
;
;-

pro bolton_biasgen, bias_files, outfile=outfile, niter=niter, sigthresh=sigthresh, outbias=outbias

; Set some defaults:
if (n_elements(sigthresh) eq 0) then sigthresh = 5.
if (n_elements(niter) eq 0) then niter = 10L

; Make sure we have something to work with:
nf = n_elements(bias_files)
if (nf lt 1) then begin
  print, 'Need at least one bias filename'
  return
end

; Get basic info about the images:
hdr0 = headfits(bias_files[0])
nxfull = sxpar(hdr0, 'NAXIS1')
nyfull = sxpar(hdr0, 'NAXIS2')
xhw = nxfull / 2
yhw = nyfull / 2
mjd = strtrim(string(sxpar(hdr0, 'MJD')), 2)
cam = strtrim(sxpar(hdr0, 'CAMERAS'), 2)

print, 'Constructing ' + cam + ' bias for MJD ' + mjd

; Construct filename if necessary:
if (n_elements(outfile) ne 1) then outfile = 'boss_pixbias-' + mjd + '-' + cam + '.fits'
print, ' Output file is ' + outfile

; Initialize arrays to hold output:
outbias = fltarr(nxfull, nyfull)
outoffset = fltarr(nf, 2, 2)
outrchi2 = fltarr(niter+1, 2, 2)

; Begin the various loops:

; We'll work with quadrant flags in x and y:
; 0 for the low numbers, 1 for the high numbers.

for xflag = 0, 1 do begin
   for yflag = 0, 1 do begin
      print, ' Initializing for quadrant ', xflag, yflag
; Indices to pick out the quadrant:
      xlo = xflag * xhw
      ylo = yflag * yhw
      xhi = (xflag + 1) * xhw - 1
      yhi = (yflag + 1) * yhw - 1
; Indices to pick out the inner quadrant of the quadrant:
      dxlo = (1-xflag) * xhw/2
      dylo = (1-yflag) * yhw/2
      dxhi = ((1-xflag) + 1) * xhw/2 - 1
      dyhi = ((1-yflag) + 1) * yhw/2 - 1
; Array to hold quadrant bias:
      bstack = replicate(0., xhw, yhw, nf)
; Loop to read in:
      for i = 0L, nf-1 do bstack[*,*,i] = (mrdfits(bias_files[i], 0, /fscale))[xlo:xhi,ylo:yhi]
; Initialize frame-wise offset and pixel mask:
      offset = fltarr(nf)
      bmask = replicate(1B, xhw, yhw, nf)
; Initialize the offset with a global median over the image quadrant:
      for i = 0L, nf-1 do offset[i] = median(bstack[*,*,i])
; Initialize the master bias with a median over frames,
; and also initialize a working array to hold the offset-subtracted
; individual bias frames:
      workdata = fltarr(xhw, yhw, nf)
      for i = 0L, nf-1 do workdata[*,*,i] = bstack[*,*,i] - offset[i]
      mbias = median(workdata, dimen=3)
; Initialize error estimate:
      rmsfluc = fltarr(nf)
      for i = 0L, nf-1 do begin
         djs_iterstat, bstack[*,*,i] - mbias - offset[i], sigma=sigma
         rmsfluc[i] = sigma
      endfor
      sigma = sqrt(mean(rmsfluc^2))
; Initialize the mask:
      for i = 0L, nf-1 do bmask[*,*,i] = (abs(bstack[*,*,i] - mbias - offset[i]) / sigma) le sigthresh
; Initialize rchi2:
      chi2 = 0.d0
      for i = 0L, nf-1 do chi2 = chi2 + total((bstack[*,*,i] - mbias - offset[i])^2 * bmask[*,*,i] / sigma^2, /double)
      rchi2 = chi2 / (total(bmask, /double) - double(xhw)*double(yhw) - double(nf))
      print, '   Initial rchi2 = ', rchi2
      outrchi2[0,xflag,yflag] = rchi2

; Now we iterate:
      for iiter = 1, niter do begin
         print, '   Iteration ', iiter, ' of ', niter
; First is a step to find the best offset, given the master bias:
         for i = 0L, nf-1 do offset[i] = total((bstack[*,*,i] - mbias) * bmask[*,*,i], /double) / total(bmask[*,*,i], /double)
; Next is a step to update the mask:
         for i = 0L, nf-1 do bmask[*,*,i] = (abs(bstack[*,*,i] - mbias - offset[i]) / sigma) le sigthresh
; Next is a step to update the offset-subtracted working array:
         for i = 0L, nf-1 do workdata[*,*,i] = bstack[*,*,i] - offset[i]
; Next is a step to update the master bias given the offsets:
         mbias = total(workdata * bmask, 3, /double) / (total(bmask, 3, /double) > 1.)
; Next is a step to update the mask:
         for i = 0L, nf-1 do bmask[*,*,i] = (abs(bstack[*,*,i] - mbias - offset[i]) / sigma) le sigthresh
; Next we compute something like a reduced chi-squared to track convergence:
         chi2 = 0.d0
         for i = 0L, nf-1 do chi2 = chi2 + total((bstack[*,*,i] - mbias - offset[i])^2 * bmask[*,*,i] / sigma^2, /double)
         rchi2 = chi2 / (total(bmask, /double) - double(xhw)*double(yhw) - double(nf))
         print, '           rchi2 = ', rchi2
         outrchi2[iiter,xflag,yflag] = rchi2
; Finally we enforce that the master bias has zero median
; in the inner quadrant of the quadrant,
; just to take out that constant d.o.f. degeneracy
; that might lead us to wander:
         mbias = mbias - median(mbias[dxlo:dxhi,dylo:dyhi])
      endfor

; Fill in for this quadrant:
      outbias[xlo:xhi,ylo:yhi] = mbias

   endfor
endfor

; Recycle the zeroth header, and add on cards to log the exposure IDs:
for i = 0L, nf-1 do begin
   spos = strpos(bias_files[i], 'sdR-', /reverse_search)
   this_id = strmid(bias_files[i], spos+4, 11)
   this_num = string(i, format='(i2.2)')
   sxaddpar, hdr0, 'EXPID'+this_num, this_id, 'ID string for exposure ' + this_num
endfor

; Write out final pixel bias (probably want more header cards at some point):
mwrfits, outbias, outfile, hdr0, /create
;writefits,'wf-'+outfile,outbias,hdr0
;stop
return
end
