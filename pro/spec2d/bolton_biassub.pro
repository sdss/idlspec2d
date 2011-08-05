;+
;
; NAME:
;  bolton_biassub
;
; PURPOSE:
;  Perform BOSS 4-amp bias subtraction using master pixel bias.
;  also trim off overscan.
;
; USAGE:
;  img = bolton_biassub(fname, biasname [, $
;         sigthresh=sigthresh, rnoise=rnoise, notrim=notrim])
;
; ARGUMENTS:
;   fname: full path qualified srR filename to read
;   biasname: name of bias file to subtract:
;   sigthresh (optional) # of sigmas threshold for pixel masking in
;     the overscan region used to solve for offset relative to pixbias.
;   notrim: set this keyword to bypass trimming of overscan regions.
;
; OUTPUTS:
;   img: bias-subtracted and (normally) overscan-trimmed image
;   rnoise (optional): 2 x 2 array of estimated RMS readnoise in ADU,
;     one value for each quadrant, in the arrangement you'd expect.
;
; NOTES:
;   Applied bias file should be output from BOLTON_BIASGEN,
;   with name like 'boss_pixbias-r?-MJDXX.fits.gz'
;
;   See algorithmic notes in BOLTON_BIASGEN.
;
; BUGS:
;   Does not complain if you mix up fname and biasname, whereas
;   it probably ought to.
;
; WRITTEN:
;  A. Bolton, U. of Utah, 2011 aug.
;
;-

function bolton_biassub, fname, biasname, sigthresh=sigthresh, rnoise=rnoise, notrim=notrim

if (n_elements(sigthresh) eq 0) then sigthresh = 4.

hdr = headfits(fname)
cam = strtrim(sxpar(hdr, 'CAMERAS'), 2)

if ((cam eq 'r1') or (cam eq 'r2')) then begin
   bx0 = 10
   bx1 = 100
   by0 = 52 ; 48 ; ASB: I don't like row 48
   by1 = 2110 ; 2111 ; ASB: middle rows suspect.
   dx0 = 119
;   dx1 = 2175 ; ASB: don't need this
   dy0 = 48
;   dy1 = 2111 ; ASB: don't need this
endif

if ((cam eq 'b1') or (cam eq 'b2')) then begin
   bx0 = 10
   bx1 = 67
   by0 = 60 ; 56 ; ASB: boosting this while boosting same for r above
   by1 = 2111
   dx0 = 128
;   dx1 = 2175 ; ASB: don't need this
   dy0 = 56
;   dy1 = 2111 ; ASB: don't need this
endif

data_img = mrdfits(fname, 0, /fscale)
bias_img = mrdfits(biasname)
rnoise = fltarr(2, 2)
nxfull = (size(data_img))[1]
nyfull = (size(data_img))[2]
xhw = nxfull / 2
yhw = nyfull / 2
bsize = size(bias_img)
if ((bsize[1] ne nxfull) or (bsize[2] ne nyfull)) then begin
   print, 'Data and bias image sizes mismatched!'
   return, 0
endif

; Loop over quadrants:

for xflag = 0, 1 do begin
   for yflag = 0, 1 do begin
; Indices to pick out the quadrant:
      xlo = xflag * xhw
      ylo = yflag * yhw
      xhi = (xflag + 1) * xhw - 1
      yhi = (yflag + 1) * yhw - 1
; Indices to pick out the bias region:
      bxlo = xflag ? nxfull - bx1 - 2 : bx0
      bxhi = xflag ? nxfull - bx0 - 2 : bx1
      bylo = yflag ? nyfull - by1 - 2 : by0
      byhi = yflag ? nyfull - by0 - 2 : by1
; Get data and bias subimages:
      data_sub = data_img[bxlo:bxhi,bylo:byhi]
      bias_sub = bias_img[bxlo:bxhi,bylo:byhi]
; Initialize offset with median:
      offset = median(data_sub - bias_sub)
; Estimate error scale:
      djs_iterstat, data_sub - bias_sub - offset, sigma=sigma, sigrej=sigthresh
; Mask discrepant pixels:
      pmask = (abs(data_sub - bias_sub - offset) / sigma) le sigthresh
; Find best offset with mask:
      offset = float(total((data_sub - bias_sub) * pmask, /double) / total(pmask, /double))
; Iterate the mask once and re-find offset:
      pmask = (abs(data_sub - bias_sub - offset) / sigma) le sigthresh
      offset = float(total((data_sub - bias_sub) * pmask, /double) / total(pmask, /double))
; Compute the RMS fluctuation:
      rnoise[xflag,yflag] = sqrt(total((data_sub - bias_sub - offset)^2 * pmask, /double) / total(pmask, /double))
; Do the bias subtraction:
      data_img[xlo:xhi,ylo:yhi] = data_img[xlo:xhi,ylo:yhi] - bias_img[xlo:xhi,ylo:yhi] - offset
   endfor
endfor

; Trim the overscan unless otherwise requested:
if (not keyword_set(notrim)) then data_img = data_img[dx0:nxfull-dx0-1,dy0:nyfull-dy0-1]

return, data_img
end
