;+
; NAME:
;   fitarcimage
;
; PURPOSE:
;   Determine wavelength calibration from arclines
;
; CALLING SEQUENCE:
;   fitarcimage, arc, arcinvvar, side, linelist, xnew, ycen, tset, invset, $
;    func=func, ncoeff=ncoeff, ans=ans, lambda=lambda, $
;    thresh=thresh, row=row, $
;    xdif_lfit=xdif_lfit, xdif_tset=xdif_tset
;
; INPUTS:
;   arc        - Extracted arc spectra with dimensions [NY,NFIBER]
;   arcinvvar  - Inverse variance of ARC
;   side       - 'red' or 'blue'; not required if ANS is set
;   linelist   - Array with wavelengths of possible arc lines
;
; OPTIONAL KEYWORDS:
;   func       - Name of fitting function; default to 'legendre'
;   ncoeff     - Number of fitting coefficients; default to 5
;   ans        - First guess for wavelength solution; what is it???
;   thresh     - Threshhold counts for significant lines;
;                default to 200 if side='blue' or 500 if side='red'
;   row        - Row to use in initial guess of wavelength solution;
;                default to (NFIBER-30)/2
;
; OUTPUTS:
;   xnew       - pixel position of lines [nfiber, nline]
;   ycen       - fiber number [nfiber, nline]
;   tset       - traceset (pix -> lambda)
;   invset     - inverse traceset (lambda -> pix)
;
; OPTIONAL OUTPUTS:
;   lambda     - returns alog10(wavelength) of good lines
;   xdif_lfit  - fit residual for individual arclines
;   xdif_tset  - fit residual of traceset
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Check for case where initial ROW is an un-illuminated fiber.
;
; INTERNAL PROCEDURES:
;   gaussprofile()
;   tset_struc()
;   lampfit()
;   corrlamps()
;   fullfit()
;
; PROCEDURES CALLED:
;   djs_median
;   djsig()
;   fit_tset
;   fitwithmx()
;   trace_crude()
;   trace_fix() ???
;   trace_fweight()
;   traceset2pix()
;
; REVISION HISTORY:
;   15-Oct-1999  Written by S. Burles, D. Finkbeiner, & D. Schlegel, APO
;-
;------------------------------------------------------------------------------

; This routine just makes a 2D kernel with a simple
; gaussian shape.   I can easily believe this can be done in easier.

function gaussprofile, width, halfsize

   if (NOT keyword_set(width)) then return, 0
   if (NOT keyword_set(halfsize)) then halfsize = long(4 * width)

   fullsize = 2*halfsize + 1
   profile = fltarr(fullsize)
   ii = lindgen(halfsize+1) 
   submm = fltarr(halfsize+1,5)
   for i=0,4 do submm[*,i] = ii + (i-2)*0.2
   denom = 1.06448 * width
   ee = (submm*submm)/(width*width/2.77254) 
   
   subprofile = exp(-ee)/denom
   for i=0,halfsize do begin
           profile[-i + halfsize] = total(subprofile[i,*])/5.0
           profile[i + halfsize] = total(subprofile[i,*])/5.0
   endfor

   return, profile
end

;------------------------------------------------------------------------------

; Define traceset structure
function tset_struc, func, ncoeff, nTrace

   tset = $      
    { func    :    func              , $
      xmin    :    0.0d               , $
      xmax    :    0.0d               , $
      coeff   :    dblarr(ncoeff, nTrace) $
    }

   return,tset
end

;------------------------------------------------------------------------------

;   Lampfit sets up data for amoeba with gaussian profile
;   At exit it resets guess coefficients by correcting
;   for best lag

function lampfit, spec, linelist, guess0, guesshi, width=width, $
 lagwidth=lagwidth, lag=lag, ftol=ftol

   common lampstuff, start, loglamlist, intensity, speccorr
   common lagstuff, llag, bestlag, bestcorr, zeroterm
   common pixstuff, pixarray, profile, middle

   zeroterm = guess0

   if (NOT keyword_set(width)) then width = 2.5
   if (NOT keyword_set(lagwidth)) then lagwidth = 300
   if (NOT keyword_set(ftol)) then ftol = 1.0e-4

   profilesize = 2 * long(width*3) + 1

   profile = gaussprofile(width)
   middle = n_elements(profile)/2

   
   ;
   ;   Guess is an array of polynomial coefficients
   ;

   npix = n_elements(spec)
   pixarray = 2.0*findgen(npix)/(npix-1) - 1.0


   nlines = (size(linelist))[1]
   loglamlist = (alog10(linelist[*,0]))[*]
   intensity = (linelist[*,1])[*]
   pix = fltarr(nlines)

   ;
   ;   Cap strong lines at 3000 counts
   ;
   nsmooth = npix+2*lagwidth
   specsmooth = convol(spec < 3000,profile)
   speccorr = fltarr(nsmooth)
   start = lagwidth

   speccorr[start:start+npix-1] = specsmooth 
   lag = lindgen(2*lagwidth) - lagwidth
   llag = lag

   p0 = (guesshi[0,*])[*]
   scale = (guesshi[1,*])[*]
   res = amoeba(ftol, p0 = p0, scale=scale, function_name='corrlamps', nmax=200)
   
   ans = [guess0,res]
   nparams = n_elements(ans)
   goffset = fltarr(nparams)

   pixoff = 2.0*bestlag/(npix-1) 
   goffset[0] = ans[1]*pixoff
   if(nparams LT 3) then return, ans-goffset

    goffset[0] = goffset[0] - ans[2]*pixoff*pixoff   
    goffset[1] = 2.0*ans[2]*pixoff   

   return, ans-goffset
end   

;------------------------------------------------------------------------------

;   Corrlamps is the raw function called by amoeba which calculates
;   the best cross_correlation and returns the corresponding maximum
;   value of the cross correlation.  Construction of a new model
;   spectrum is done using the coefficients in (a).   

function corrlamps, a

   common lampstuff, start, loglamlist, intensity, speccorr
   common lagstuff, lag, bestlag, bestcorr, zeroterm
   common pixstuff, pixarray, profile, middle

   guess = [zeroterm, a]
   loglam = poly(pixarray,guess)
   model = speccorr*0.0
   nlines = n_elements(loglamlist)
   for i=0,nlines-1 do begin
           diff = loglam - loglamlist[i]
      val = min(abs(diff), place)

      if(place GT 0 AND place LT 2047) then begin   
        x1 = place+start-middle
        x2 = place+start+middle
        model[x1:x2] = model[x1:x2] + profile*intensity[i]
        endif
     endfor

   bestcorr = 0.0
   bestlag = 0
   model = sqrt(model)


   if (total(model) GT 0) then begin
      res = c_correlate(model,speccorr, lag)
      bestcorr = max(res,val)
      bestlag = lag[val]
   endif
   print, guess, bestcorr, bestlag

   return, 1.0 - bestcorr
end

;------------------------------------------------------------------------------

; Fullfit performs multiple (at least 2) iterations of lampfit
; One might think of this as a first coarse iterations,
; ans subsequent fine iterations.

function fullfit, spec, linelist, guess

   common lagstuff, lag, bestlag, bestcorr, zeroterm

   guess0 = guess[0]
   p0 = guess[1:*]
   scale = abs(p0) * 0.3d
   scale[0] = scale[0] * 0.05d

   first = lampfit(spec, linelist, guess0, transpose([[p0],[scale]]), $
      width = 25.0d, lagwidth=250d, ftol=1.0d-2)
   
   final = first

   while (abs(bestlag) GT 5)  do begin
      guess0 = final[0]
      p0 = final[1:*]
      scale = abs(p0)*0.1
      scale[0] = scale[0]*0.2
      final = lampfit(spec, linelist, guess0, transpose([[p0],[scale]]), $
       width = 10.0d, lagwidth=100d, ftol=1.0d-4)
   endwhile

   return, final
end

;------------------------------------------------------------------------------
function fitarcimage, arc, arcinvvar, side,linelist, xnew, ycen, tset, invset, $
 func=func, ncoeff=ncoeff, ans=ans, lambda=lambda, $
 thresh=thresh, row=row, $
 xdif_lfit=xdif_lfit, xdif_tset=xdif_tset

   common lagstuff, lag, bestlag, bestcorr, zeroterm

   ;---------------------------------------------------------------------------
   ; Preliminary stuff
   ;---------------------------------------------------------------------------

   if (NOT keyword_set(ans)) then begin
      if (side NE 'blue' AND side NE 'red') then $
       message, "SIDE must be set to 'blue' or 'red' if ANS not specified"
   endif
   if (NOT keyword_set(func)) then func = 'legendre'
   if (NOT keyword_set(ncoeff)) then ncoeff = 5
   if (NOT keyword_set(ans)) then ans = 0
   if (func EQ 'legendre') then function_name = 'flegendre'
   if (func EQ 'chebyshev') then function_name = 'fchebyshev'

   t_begin = systime(1)

   icoeff = ncoeff + 2

   ndim = size(arc, /n_dim)
   if (ndim NE 2) then $
    message, 'Expecting 2-D arc image'
   dims = size(arc, /dim)
   npix = dims[0]
   nfiber = dims[1]
   if (total(dims NE size(arcinvvar, /dim))) then $
    message, 'ARC and ARCINVVAR must have same dimensions'

   if (NOT keyword_set(row)) then row = (nfiber-30)/2

   if (NOT keyword_set(thresh)) then begin
      if (side EQ 'blue') then thresh = 200
      if (side EQ 'red') then thresh = 500
   endif

   ;---------------------------------------------------------------------------
   ; First trace
   ;---------------------------------------------------------------------------
   ; One might want to change nave and nmed for first pass???
   ; Now find the initial wavelength solution if ANS is not passed

   if (ans[0] EQ 0) then begin

      ; Extract one spectrum from the 5 spectra around fiber number ROW
      ; by taking the median value at each wavelength

      spec = djs_median(arc[*,row-2:row+2], 2)

      ; guess is needed to give fullfit (and lampfit) and initial starting
      ; point for wavelength solutions.

      guess = 0
      if (side EQ 'blue') then guess = [3.68, -0.106, -0.005, 0.005]
      if (side EQ 'red') then guess = [3.87, 0.10, -0.003, -0.001]   

      guess = double(guess)
      ans   = fullfit(double(spec), linelist, guess)

      if (bestcorr LT 0.5) then begin
            print, 'Initial wavelength solution looks suspicious'
	    return, 1
      endif

    endif

   ; Now store best log-lambda solutions in tset

   tset   = tset_struc(func, ncoeff, nfiber)
   invset = tset_struc(func, icoeff, nfiber) 
   tset.xmin = 0.0
   tset.xmax = 1.0*(npix - 1)

   if (side EQ 'blue') then begin
       invset.xmin = 3.57
       invset.xmax = 3.80
   endif
   if (side EQ 'red') then begin
       invset.xmin = 3.75
       invset.xmax = 3.97
   endif

   ymid = 0.5*(invset.xmax + invset.xmin)
   yrange = invset.xmax - invset.xmin
   xx = dindgen(2048)
   pixarray = 2.0d0 * dindgen(npix) / (npix-1) - 1.0d0

   ; Evaluate polynomial returned by fullfit     

   yy = poly(pixarray,ans)

   yynorm = 2.0d0 * (yy-ymid) / yrange
   invans = func_fit(yynorm, xx, icoeff, function_name=function_name)

   uselines = where(linelist[*,2] GT 0.0)
   if (uselines[0] EQ -1) then $
    message, 'No unblended lines in linelist'

   uselambda = alog10(linelist[uselines,0])
   uselambdanorm = 2.0d * (uselambda - ymid) / yrange

   if (func EQ 'legendre') then $
    xstart = flegendre(uselambdanorm,icoeff) # invans
   if (func EQ 'chebyshev') then $
    xstart = fchebyshev(uselambdanorm,icoeff) # invans

   inimage = where(xstart GE 0.0 AND xstart LE 2047.0)
   if (inimage[0] NE -1) then xstart = xstart[inimage]

   ; Allow for a shift of up to 2 pixels in the initial centers,
   ; but only 0.3 pixels while tracing

   xcen = trace_crude(arc, yset=ycen, nave=1, nmed=1, xstart=xstart, $
    ystart=row, maxshifte=0.3d, maxshift0=2.0d)
 
   ;
   ; Iterate with trace_crude
   ;
   xcen = trace_crude(arc, yset=ycen, nave=1, nmed=1, xstart=xcen[row+30,*], $
    ystart=row+30, maxshifte=0.3d, maxshift0=2.0d)

   ; Iterate the flux-weighted centers
   xnew = trace_fweight(arc, xcen, ycen)

   ; Make use of the errors??? - Seems to just mess things up???
   ; Well... the reason for that is satured lines, which return infinite errors
;   xnew = trace_fweight(arc, xcen, ycen, invvar=arcinvvar, xerr=xerr)

   bad = where(abs(xnew-xcen) GT 3.0)
   if (bad[0] NE -1) then xnew[bad] = xcen[bad]

   ntempTrace = (size(xnew))[2]
   y2 = findgen(nfiber)

   ;---------------------------------------------------------------------------
   ; Reject saturated lines
   ;---------------------------------------------------------------------------
   ; Reject lines with more than 2% of pixels saturated within 1 pixel 
   ; of all line centers (in dispersion direction). 

   isbad = lonarr(ntempTrace)
   for i=0, ntempTrace-1 do begin
      inew = round(xnew[*,i])
      arcvarlist=[arcinvvar[inew,y2],arcinvvar[inew+1,y2], $
      arcinvvar[inew-1,y2]]
      wheresat=where(arcvarlist EQ 0, nsat)
      fracbad = nsat/float(nfiber*3)
      isbad[i] = (fracbad GT 0.02)
      if isbad[i] then print, 'Discarding trace',i, $
       '    fraction bad', fracbad
   endfor

   goodind = where(1-isbad, numgood)
   print, 'I found some lines: ', numgood
   if goodind[0] EQ -1 then $
    message,'No good arc lines!!!'

   ; Trim linelist

   xnew = xnew[*,goodind]
   ycen = ycen[*,goodind]
   lambda = uselambda[inimage[goodind]]

   ;---------------------------------------------------------------------------
   ; Do the first traceset fit
   ;---------------------------------------------------------------------------

   fit_tset, xnew, ycen, lambda, goodlines, tset, invset ;, /linesearch
   print,'Pass 1 complete'

   ; Keep only "good" lines
   ; Currently this means only 3 lines can be rejected per bundle
   possible = (size(goodlines))[2]
   nlines = (size(goodlines))[1]
   if (possible NE 320) then $
    message, "can't figure out bundle test"
   testg = reform(goodlines, nlines, 20, 16)   
   gind = where(total(total(testg EQ 0,2) GT 3,2) EQ 0, gindct)
   if (gindct LT 5) then begin
    print, 'FIT_TSET failed'
    return, 2
   endif
   goodlines = 0

   nline = n_elements(gind)    ; number of good lines
   xnew = xnew[*,gind]
   lambda = lambda[gind]
   ycen = ycen[*,gind]

   fit_tset, xnew, ycen, lambda, goodlines, tset, invset
   print, 'Pass 2 complete'

   ; Fit arc lines subtracting out scatter term

   xmeasured = xnew
   xnew = fitwithmx(invset, lambda, xmeasured)

   fit_tset, xnew, ycen, lambda, goodlines, tset, invset
   print, 'Pass 3 complete'
   ; now tset, invset contain "all" the information in the data.

   ;---------------------------------------------------------------------------
   ; Quality Assurance
   ;---------------------------------------------------------------------------

   ; pixel positions derived from the traceset

   tset_pix = traceset2pix(invset, lambda)

   xdif_tset = (xmeasured-tset_pix)  ; difference between measured line 
                                     ;  positions and fit positions
   xdif_lfit = (xmeasured-xnew)      ; dif between measured line positions
                                     ;  and best fit for each line

   print
   print, 'Arcline fit summary'
   print, 'All sigma values are in millipixels'
   print, format='(71("-"))'
   for k=0, n_elements(gind)-1 do $
      print,'Arcline',k,':  lambda =',10.^lambda[k], $
            '    sig_lfit =', djsig(1e3*xdif_lfit[*,k]), $
            '    sig_tset =', djsig(1e3*xdif_tset[*,k]), $
            format='(A,I3,A,F8.2,A,F7.2,A,F7.2)'

   highones = where(lambda GT alog10(8000.), highct)
   if (nline LT 6) then $
     message, 'WARNING: only '+string(nline)+ ' good arclines found',/cont

   print, 'Found ', nline, ' good arc lines'
   if (highct GT 0) then print, '----', highct, ' are above 8000 A'
   print,'> FITARCIMAGE: ',systime(1)-t_begin, ' seconds elapsed', $
    format='(A,F8.2,A)'

   return, 0
end
;------------------------------------------------------------------------------
