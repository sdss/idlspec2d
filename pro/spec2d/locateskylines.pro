;------------------------------------------------------------------------------
;+
; NAME:
;   locateskylines
;
; PURPOSE:
;   Tweak wavelength calibration using skylines
;
; CALLING SEQUENCE:
;   locateskylines, skylinefile, fimage, ivar, $
;	tset_arc, invset, tset_tweak, invset_tweak, $
;       xsky, ysky, skywaves, lambda=lambda, errcode=errcode
;
; INPUTS:
;   skylinefile - filename of skyline file
;   fimage      - flattened image containing skylines (2k x 2k)
;   ivar        - inverse variance of fimage
;   tset_arc    - traceset from arclines (lambda -> pix)
;   invset      - inverse arc traceset (lambda -> pix)
;
; OPTIONAL KEYWORDS:
;   lambda      - override skyline file (log10 Angstroms)
;
; OUTPUTS:
;   xsky        - pixel position of lines [nfiber, nline]
;   ysky        - fiber number [nfiber, nline]
;   skywaves    - wavelengths used (Angstroms)
;   tset_tweak  - traceset (pix -> lambda) tweaked to sky
;   invset_tweak- inverse traceset (lambda -> pix) tweaked to sky
;
; OPTIONAL OUTPUTS:
;   errcode     - returns errcode (see below)
;
; COMMENTS:
;   Error codes:
;  -1 - <unknown error>
;   0 - Everything fine
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   15-Oct-1999  Written by S. Burles, D. Finkbeiner, & D. Schlegel, APO
;-
;------------------------------------------------------------------------------

PRO locateskylines, skylinefile, fimage, ivar, $
	tset_arc, invset, tset_tweak, invset_tweak, $
        xsky, ysky, skywaves, lambda=lambda, errcode=errcode

	if keyword_set(lambda) then begin 
	  skywaves=10.d^lambda
        endif else begin 
          print,'Reading ',skylinefile
          readcol,skylinefile,skywaves,/silent,form='d'
 	endelse

	nskyline = n_elements(skywaves)
	nfiber = (size(fimage))[2]
	ncoeff = (size(invset.coeff))[1]
        npix   = (size(fimage))[1]

;---------------------------------------------------------------------------
;	Make xsky, ysky, pairs from skywaves
;---------------------------------------------------------------------------
;       xarc contains skyline positions based on arc fit

	ysky = dindgen(nfiber)#(dblarr(nskyline)+1)
	xarc = traceset2pix(invset, alog10(skywaves))
;---------------------------------------------------------------------------
; Discard lines that are out of bounds
;---------------------------------------------------------------------------
	good = total(((xarc LE 0) OR (xarc GE npix)),1) EQ 0 
        gind = where(good)
        if gind[0] eq -1 then message, 'No sky lines!!!'

; trim list
        xarc=xarc[*,gind]
	ysky=ysky[*,gind]
        skywaves=skywaves[gind]


	fmed=fimage-fimage
	for i=0,nfiber-1 do fmed[*,i]=median(fimage[*,i],17)

	xskytmp = trace_fweight(fimage, xarc, ysky, invvar=ivar) 
; iterate trace_fweight ?!?
	xsky =  trace_fweight(fimage, xskytmp, ysky, invvar=ivar) 

;	xsky =  trace_gauss(fimage, xarc, ysky, invvar=ivar) 

;---------------------------------------------------------------------------
;  Bonehead Statistics
;---------------------------------------------------------------------------

; xsky now contains positions of skylines.  We have no idea if they 
;   are any good unless we do some checks.


	nskyline=(size(xsky))[2]
	mean  = fltarr(nskyline)
	sigma = fltarr(nskyline)
	xres  = xarc-xsky

; Amp 1
	for i=0,nskyline-1 do begin 
	  djs_iterstat,xres[3:150,i],mean=mn,sigma=sig
	  mean[i]=mn
	  sigma[i]=sig
        endfor 
        mean1=mean & sigma1=sigma

; Amp 2
	for i=0,nskyline-1 do begin 
	  djs_iterstat,xres[170:316,i],mean=mn,sigma=sig
	  mean[i]=mn
	  sigma[i]=sig
        endfor 
        mean2=mean & sigma2=sigma

; Both amps
	for i=0,nskyline-1 do begin 
	  djs_iterstat,xres[*,i],mean=mn,sigma=sig
	  mean[i]=mn
	  sigma[i]=sig
        endfor 

;---------------------------------------------------------------------------
;  Quality assurance
;---------------------------------------------------------------------------

        pmulti = !p.multi
        !p.multi = [0,1,2]
	plot,skywaves,mean1,xr=[5500,9500],yr=[-.3,.3]+median(mean),/yst, $
	    xtit='Wavelength (A)', ytit='Delta x (Pix)', /xst, $
            title='Sky line residual'
	errplot,skywaves,mean1-sigma1,mean1+sigma1
	oplot,skywaves,mean2,line=1

; Reject outliers with "drop dead" requirements
	good = (abs(mean) LT 3) AND (abs(sigma) LT 0.2)
	gind = where(good, nline)	
        if gind[0] eq -1 then message, 'No sky lines!!!'
	print,nline, ' good lines - ',fix(total(good EQ 0)), ' rejected'


; Trim list again 
        mean  = mean[gind]
        sigma = sigma[gind]
        xarc  = xarc[*,gind]
	ysky  = ysky[*,gind]
        skywaves=skywaves[gind]
      
        if nline GE 1 then nord=0
        if nline GE 3 then nord=1
	if nline GE 7 then nord=2

        print,'Fit order:',nord

	flexcoeff=polyfitw(skywaves,mean,1./sigma^2,nord,yfit)
        oplot,skywaves,yfit,line=2

	infostr=string('Dispersion:', stdev(mean-yfit)*1E3,' mpix', $
	   format='(A,F7.1,A)')

	plot,skywaves,mean-yfit,xr=[5500,9500],yr=[-.2,.2], /xst, $
	    xtit='Wavelength (A)', ytit='Delta x (Pix)', $
            title='After flexure correction'
	errplot,skywaves,mean-yfit-sigma,mean-yfit+sigma
        xyouts, 0.95, 0., systime(), /normal, align=1, chars=0.5
        xyouts, 0.05, 0., infostr, /norm
	print, infostr
	

stop


	traceset2xy, tset_arc, xpos, ypos


        fit_tset, xnew, ycen, lambda, goodlines, tset_tweak, invset_tweak


	!p.multi=pmulti


  return

END
;---------------------------------------------------------------------------
;  locateskylines.pro EOF
;---------------------------------------------------------------------------
