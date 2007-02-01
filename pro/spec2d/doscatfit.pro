;+
; NAME:
;   doscatfit
;
; PURPOSE:
;   Fit for the coefficients in the red or blue scattering models, or
;   generate a scattering surface based on such a model for a given image.
;
;
; CALLING SEQUENCE:
;   coeffs = doscatfit(image, arcfile, trace=, keeprows=,
;                      addrows=, notchcols=, outfile=
;
; INPUTS:
;   image     - Raw SDSS file name
;   arcfile   - SDSS
;
; OPTIONAL KEYWORDS:
;    trace    - Only fit around the Nth (1-based!) trace.
;    keeprows - Only consider the given rows in the image.
;    addrows  - account for N pseudo-rows of illumination at the top of the image.
;    notchcols- zero out the weights for the 2N+1 columns from the
;               middle of each trace.
;    outfile  - name of file to write the interpolated per-row
;               coefficients to.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; BUGS:
;
; PROCEDURES CALLED:
;                        FIXME - CPL
;   djs_filepath()
;   djs_iterstat
;   fileandpath()
;   findopfile()
;   fits_purge_nans
;   fits_wait()
;   headfits()
;   idlspec2d_version()
;   idlutils_version()
;   lookforgzip()
;   rdss_fits()
;   sphdrfix
;   splog
;   sxaddpar
;   sxpar()
;   writefits
;   yanny_free
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; DATA FILES:
;   $IDLSPEC2D_DIR/examples/opConfig*par
;   $IDLSPEC2D_DIR/examples/opECalib*par
;   $IDLSPEC2D_DIR/examples/opBC*par
;   $SPECFLAT_DIR/biases/pixbias*.fits
;   $SPECFLAT_DIR/flats/pixflat*.fits

function doscatfit, imgfile, arcfile, flatfile=flatfile,trace=trace, $
                    addrows=addrows, notchcols=notchcols,$
                    ymodel=ymodel, evalwith=evalwith, img=img, ximg=ximg, wset=wset, spec=spec, $
                    extract_only=extract_only, maxkernr=maxkernr, $
                    sigma=sigma, bandrows=bandrows, dosmall=dosmall,doplot=doplot, $
                    outfile=outfile

    if not keyword_set(notchcols) then notchcols=4
    if not keyword_set(bandrows) then bandrows=1
    if not keyword_set(sigma) then sigma=0.8

    sdssproc,imgfile,img,imgivar,indir='/u/dss/rawdata',hdr=hdr
    sdssproc,arcfile,arc,arcivar,indir='/u/dss/rawdata'
    cam = strmid(sxpar(hdr, 'CAMERAS'), 0, 2)
    color = (strmid(cam,0,1) EQ 'r') ? 'red' : 'blue'

    if keyword_set(flatfile) then begin
       sdssproc,flatfile,flat,flativar,/nobleed,indir='/u/dss/rawdata'
    end else begin
       flat = img
       flativar = imgivar
    end
                                ; Place the spectrum and measure the
                                ; flux in the same we would measure it
                                ; in a real extraction. Or close.
    ximg = trace_crude(flat, flativar, thresh=1000., yset=yimg, xerr=xerr)
    extract_image, flat, flativar, ximg, 1.0, spec, specivar, $
                   proftype=3, npoly=3, wfixed=[1,1,1], $
                   ansimage=ansimage, highrej=30., ymodel=ymodel0

                                ; Get a reasonably close wavelength
                                ; solution. Note that fitarcimage only
                                ; really work with 320 fibers, so we
                                ; have to hack it a bit.
    extract_image, arc, arcivar, ximg, 1.0, arcspec, arcspecivar, $
                  proftype=3, npoly=3, wfixed=[1,1], relative=1, $
                  highrej=30.

    if keyword_set(trace) then begin 
       fitarcimage, arcspec, arcspecivar, xpeak, ypeak, wset, ncoeff=5, $
                    color=color,nmed=1,row=trace-1,/dofirst
    end else begin
       nspec = (size(ximg, /dim))[1]
       dofirst = nspec NE 320 
       fitarcimage, arcspec, arcspecivar, xpeak, ypeak, wset, ncoeff=5, $
                    color=color,nmed=1,dofirst=dofirst;,maxdev=1e-4
    end

    traceset2xy, wset, wx, wy
    splog, "wavelength sanity mins: ", 10^wy[0,*]
    splog, "wavelength sanity maxs: ", 10^wy[2047,*]

    if keyword_set(extract_only) then return, 0

                                ; This is an important but imperfect
                                ; step: the image flux drives the
                                ; scattering surface, so saturated
                                ; pixels and bad columns need to have
                                ; some reasonable values.
                                ; [ The various reject_cr routines
                                ; don't appear to work well on
                                ; spectro frames. ]
    ;fitimg = psf_reject_cr ????
    fitimg = djs_maskinterp(img, imgivar EQ 0, iaxis=0)

                                ; Seed with some plausible defaults.
    if color EQ 'red' then begin
       evalfunc = 'redhalo'
       if not keyword_set(maxkernr) then maxkernr=256
       fitrows = 2048
       dobands = 1
       dosmall = 1

       tp = replicate({value:0.0, fixed:0, limited:[0,0], limits:[0.0,0.0]}, 4)
       if cam EQ 'r1' then begin
          tp[2].value=78.4
          tp[3].value=8.22
       end else if cam eq 'r2' then begin
          tp[2].value=80.0
          tp[3].value=9.37
       end
    end else begin
       evalfunc = 'bluescat'
       if not keyword_set(maxkernr) then maxkernr=32
       fitrows = 32
       dobands = 0
       if not keyword_set(evalwith) then bandrows = fitrows + 2*maxkernr

       tp = replicate({value:0.0, fixed:0, limited:[0,0], limits:[0.0,0.0]}, 7)
                                ; tp[0] -> DC
                                ; tp[1] -> if non-zero, core sigma
       tp[2].value = sigma
       tp[2].fixed = 1
       if cam EQ 'b1' then begin
          tp[3].value = 0.005   ; b - wing scale
          tp[4].value = 1.5     ; p - wing exponent
          tp[4].limited[1] = 1
          tp[4].limits[1] = 2.0

          tp[5].value = 10.0    ; sigma2 - second core sigma
          tp[6].value = 0.00005 ; a2 - second core scale
       end else begin
          tp[3].value = 0.001
          tp[4].value = 1.5
          tp[4].limited[1] = 1
          tp[4].limits[1] = 2.0

          tp[5].value = 10.0
          tp[6].value = 0.0001
       end
    end
                                ; If a single trace is specified,
                                ; extract the smallest region we can. 
    ntrace = (size(ximg, /dim))[1]
    if keyword_set(trace) then begin
       trace -= 1
       if trace eq 0 then begin
          x0 = 0
          x1 = min((ximg[*,1] + ximg[*,0])/2)
       end else if trace eq (size(spec, /dim))[1] - 1 then begin
          x0 = long(max((ximg[*,ntrace-1] + ximg[*,ntrace-2])/2))
          x1 = 2047
       end else begin
          x0 = long(max((ximg[*,trace] + ximg[*,trace-1])/2))
          x1 = long(min((ximg[*,trace+1] + ximg[*,trace])/2))
       end
       ntrace = 1
    end else begin
       trace = indgen(ntrace)
       x0 = 0
       x1 = 2047
    end

    if x0 mod 2 eq 1 then x0 += 1
    if x1 mod 2 eq 0 then x1 -= 1
    splog, "trimming to trace, x0, x1: ", trace, x0, x1
    
    fitreg = fitimg[x0:x1, *]
    fitivar = imgivar[x0:x1, *]
    fitx = ximg[*,trace] - x0
    fitlam = wy[*,trace]
    fitflux = spec[*,trace]

                                ; Flatten to a single wavelength per row.
    if size(fitlam, /n_dim) GT 1 then $
       fitlam = djs_median(fitlam,2)

                                ; Scattered light comes in from beyond
                                ; the detector region of the chip. We
                                ; only _really_ care about the red end
                                ; of red, and slightly care about the
                                ; blue end of blue. Fabricate a lie
                                ; for the evaluator by extending the
                                ; wavelength solution and copying the
                                ; last row of the detector. Actually,
                                ; we want to extend on the sides as
                                ; well, and probably attenuate the
                                ; extended flux.
    if keyword_set(addrows) then begin
       splog, "adding rows: ", addrows

       imgsize = size(fitreg,/dim)
       nr = imgsize[1]
       nc = imgsize[0]

       fitreg = [[fitreg], [rebin(fitreg[*,nr-1], [nc, addrows])]]
       fitivar = [[fitivar], [rebin(fitivar[*,nr-1], [nc, addrows])]]
       tfitx = transpose(fitx)
       fitx = transpose([[tfitx], [rebin(tfitx[*,nr-1], [ntrace, addrows])]])

       if keyword_set(evalwith) and size(evalwith, /n_dim) GT 1 then begin
          nparams = (size(evalwith, /dim))[0]
          evalwith = [[evalwith], [rebin(evalwith[*,nr-1], [nparams, addrows])]]
       end

       dlam = fitlam[nr-1] - fitlam[nr-2]
       elam = findgen(addrows)/(addrows-1) * (dlam * addrows)
       elam += fitlam[nr-1]+dlam
       fitlam = [fitlam, elam]
    end

    if keyword_set(keeprows) then begin
       splog, "trimming to keeprows: ", keeprows

       r0 = keeprows[0] - maxkernr
       rr0 = (r0 LT 0) ? (maxkernr+r0) : maxkernr
       r0 = r0 > 0
       r1 = (keeprows[1] + maxkernr) < 2047
       
       realkeeprows = [rr0, rr0 + keeprows[1]-keeprows[0]]
       fitreg = fitreg[*, r0:r1]
       fitivar = fitivar[*, r0:r1]
       fitx = fitx[r0:r1, *]
       fitlam = fitlam[r0:r1]
       fitflux = fitflux[r0:r1, *]
    end

                                ; We have a very poor idea of what the
                                ; cores look like, and have no chance
                                ; of modelling them. So deweight them
                                ; completely. 
    tivar = fitivar
    if keyword_set(notchcols) then begin
       nr = (size(tivar, /dim))[1]
       for t=0,ntrace-1 do begin
          for l=0,nr-1 do begin
             c = round(fitx[l, t])
             tivar[(c-notchcols)>0:(c+notchcols)<2047, l] = 0
          endfor
       endfor
    end

                                ; Pass the evaluator a single
                                ; wavelength per row, in microns. 
    fitmicrons = 10^fitlam / 10000
    
                                ; Fit or evaluate the surface....
    if keyword_set(evalwith) then begin
       keeprows = (size(img, /dim))[1]
       ymodel = gen_scatter(fitreg, fitmicrons, evalwith, $
                            evalfunc=evalfunc, bandrows=bandrows, $
                            maxkernr=maxkernr, keeprows=[0,keeprows-1])
       allres = evalwith
    end else begin
       nfitr = (size(fitreg,/dim))[1]
       nfitbands = nfitr / fitrows
       allres = fltarr(size(tp,/dim)+2, nfitbands)

                                ; Seed the 0-th order background. Ugh.
       djs_iterstat,fitreg, invvar=tivar, median=dc, sigrej=2
       tp[0].value = dc * 0.5   ; DC 
       tp[0].limited[0] = 1
       tp[0].limits[0] = 0.0
       tp[1].value = 0          ; Do not notch core.
       tp[1].fixed = 1          ; 

                                ; Fit a set of coefficients for each
                                ; band of rows. We might actually want
                                ; to do this for, say, each bundle
                                ; region, or some other 2d gridding.
       for f_i=0,nfitbands-1 do begin
          r0 = f_i*fitrows - maxkernr
          rr0 = (r0 LT 0) ? (maxkernr+r0) : maxkernr
          r0 = r0 > 0
          r1 = ((f_i+1)*fitrows + maxkernr-1) < (nfitr-1)

          allres[0, f_i] = r0
          allres[1, f_i] = r1
          bandkeeprows = [rr0, rr0 + fitrows-1]
          rfitreg = fitreg[*, r0:r1]
          rtivar = tivar[*, r0:r1]
          rfitmicrons = fitmicrons[r0:r1, *]

          if dobands EQ 0 then begin
             rbandrows = r1-r0+1
             rmaxkernr = rbandrows/2
          end else begin
             rbandrows = bandrows
             rmaxkernr = maxkernr
          end

          fa = {dosmall:keyword_set(dosmall),bandrows:rbandrows, $
                evalfunc:evalfunc, maxkernr:rmaxkernr, keeprows:bandkeeprows}
          res = mpfit2dfun('gen_scatter', rfitreg, rfitmicrons, rfitreg, 1/sqrt(rtivar), $
                           functargs=fa, parinfo=tp, yfit=ymodel, $
                           status=status, dof=dof, bestnorm=bestnorm, errmsg=errmsg, $
                           /quiet, niter=niter)
          allres[2, f_i] = res
          splog, "fit1: ", r0, r1, niter, bestnorm, bestnorm/dof
          splog, "fit2: ", res

          if keyword_set(doplot) then begin
             tpe = res
             tpe[1] = tp[2].value
             smod = gen_scatter(rfitreg, rfitmicrons, tpe, bandrows=4, $
                                evalfunc=evalfunc, maxkernr=rbandrows/2, keeprows=bandkeeprows)
             pres,rfitreg[*,rr0:rr0+fitrows-1],smod,ivar=rtivar[*,rr0:rr0+fitrows-1],$
                  rows=(bandkeeprows[1]-bandkeeprows[0])/2,xrange=[800,1200]
          end

                                ; Pass our fit on to our neighbor.
          tp.value = res
          if color eq 'blue' then begin
             tp[0].value *= 0.75
             tp[3].value *= 0.75
             tp[4].value = 1.0
             tp[6].value *= 2
          end

       end
    end

    if keyword_set(outfile) then begin
                                ; Save per-row coefficients. Do not
                                ; bother with the DC term or the
                                ; "bookkeeping" terms.
       if size(allres, /n_dim) GT 1 then begin
          s = size(allres[4:*,*], /dim)
          ncoeff = s[0]
          nrows = s[1]
          
          outres = fltarr(ncoeff, 2048)
          for i=0,ncoeff-1 do begin
             outres[i,*] = interpol(allres[4+i,*], (allres[0,*]+allres[1,*])/2, findgen(2048))
          end
       end else begin
          outres = allres[4:*]
       end

       sxaddpar,hh,'imgname',imgfile,'Filename of image.'
       sxaddpar,hh,'arcname',arcfile,'Filename of arc.'
       sxaddpar,hh,'notchrad',notchcols,'Columns/trace deweighted.'
       sxaddpar,hh,'fitrows',fitrows,'Rows used for fit.'
       if keyword_set(addrows) then $
          sxaddpar,hh,'addrows',addrows,'Rows added to image end'
       
       mwrfits, outres, outfile, hh, /create
    end

    return, allres
 end
