;+
; NAME:
;   rm_local_qsofit
;
; PURPOSE:
;   This routine is much like rm_qsofit, but instead of a global 
;   continuum fit, it only fits a simple local continuum underneath 
;   the line. 
;
; CALLING SEQUENCE:
;   rm_local_qsofit, obs_wave, flux, err, z, [ra=,dec=,/psplot,/fits]
; ------------------------------------------------------------------

pro rm_local_qsofit, lam0, flux0, err0, z, ivar0=ivar0 $
      , ra=ra,dec=dec, deredden=deredden $
      , emparfile=emparfile $
      , input_fitmask=input_fitmask $ ; input fitmasks, e.g., custom absorption masks
      , poly_ord=poly_ord, rej_iter=rej_iter $
      , add_noise=add_noise, auto_abs_rej = auto_abs_rej $
      , conti_auto_abs_rej = conti_auto_abs_rej $
      , psplot=psplot,fits=fits,output_name=output_name,outdir=outdir,objtag=objtag $
      , plot_vel=plot_vel, nsmooth=nsmooth, append=append, silent=silent $
      , sdss_name=sdss_name, para=para, diet=diet, noplot=noplot, _extra=extra


   ; Define line fitting tolerance in mpfitfun, 
   ; which helps for some objects in line fits (i.e., Hbeta and Halpha)
   if ~keyword_set(xtol) then xtol=1d-20
   if ~keyword_set(ftol) then ftol=1d-15
  
   if n_elements(plot_vel) eq 0 then plot_vel=1
   if ~keyword_Set(outdir) then outdir = '/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/local_fit/'

   ; Define constants
   cs = 2.9979246d5  ; speed of light, km/s
   if n_elements(z) eq 0 then z = 0. ; default is in restframe

   ; default is to auto reject absorption from the first linefit
   if n_elements(auto_abs_rej) eq 0 then auto_abs_rej = 1L

   ; Number of linefit absorption rejection iteration
   if n_elements(rej_iter) eq 0 then rej_iter = 2L

   lam = lam0 & flux = flux0
   if keyword_set(ivar0) then begin
      err = 0.*ivar0
      ind = where(ivar0 NE 0.)
      err[ind] = 1./sqrt(ivar0[ind])
   endif else err = err0

   npix = n_elements(flux)

   if keyword_set(add_noise) then begin ; add gaussian noise
      flux = flux + randomn(seed, npix)*err
   endif

      ; deredden Galactic extinction if asked
   if keyword_set(deredden) then begin
      if not keyword_Set(silent) then splog, $
         'Dereddening spectrum using the SFD map and CCM extinction curve'
      dereddening_spec, lam, flux, err=err, ra=ra, dec=dec $
         , dered_flux = dered_flux, dered_err = dered_err
      flux = dered_flux & err = dered_err
   endif

   ivar = dblarr(n_elements(err))
   indd = where(err gt 1d-6)
   ivar[indd] = 1./(err[indd])^2

   ; mask out 3-sigma lower outlier of the 20-pix smoothed spectrum, 
   ; to reduce the effects of absorption
   fitmask = ivar NE 0.0

   ; Mask outliers and absorption in the emission lines by finding 5-sigma
   ; lower outliers from 20 pix smoothed spectrum
   if n_elements(outlier_mask) eq 0 then outlier_mask = 1L
   if keyword_set(outlier_mask) then begin
      disp_vec = (alog10(lam)-shift(alog10(lam), 1))[1:*]
      disp = djs_median(disp_vec)
      NRES = 20L   ; 20L
      bkspace = NRES*disp
      spec_set = bspline_iterfit(alog10(lam), flux $
                             , invvar = ivar*fitmask $
                             , bkspace = bkspace $
                             , yfit = flux_spline, maxiter = 10 $
                             , upper = 0, lower = 5, nord = 3 $
                             , maxrej = 50, outmask = outmask, /silent)
      ; exclude redward of restframe 3000A, since little absorption expected there.
      ind_exclude = where(lam/(1.+z) ge 3000.)
      if ind_exclude[0] ne -1 then outmask[ind_exclude]=1
      fitmask = fitmask*outmask
   endif

   ; Add additional mask bits from input
   if keyword_set(input_fitmask) then begin
      if n_elements(input_fitmask) eq n_elements(lam) then $
         fitmask = fitmask*input_fitmask else $
         splog, 'input_fitmask dimension does not match lam grid.'
   endif

   ; Shift to restframe
   wave = lam/(1.0+z)

   ; Trim the spectrum in the windows specified by keywords wave_range
   if keyword_Set(wave_range) then begin
      ind = where(wave ge wave_range[0] and wave le wave_range[1])
      wave = wave[ind] & flux = flux[ind] & err = err[ind]
      fitmask = fitmask[ind]
   endif

   ; generate the output structure
   output = create_struct('z',z,'wave', wave, 'flux', flux, 'err', err, $
      'fitmask', fitmask)

   ; ----------------
   ; Now start to fit the emission/absorb lines defined in emparfile
   line_flux = flux ; since we don't subtract the continuum in this fit

   ; Setup line fitting parameters
   if keyword_set(emparfile) then emline_file=emparfile else $
    emline_file = getenv('IDLRM_DIR')+'/etc/qsoline_few.par'
   yanny_read, emline_file, pdata
   linelist = *pdata[0]
   yanny_free, pdata
   nline_all = n_elements(linelist)
   ncomp_all = n_elements(uniq(linelist.compname))
   ngauss_tot = total(linelist.ngauss)
   parinfo = replicate({value:0.D, fixed:0, limited:[0,0],limits:[0.D,0],tied:''}, $
          ngauss_tot*3 + ncomp_all*2)
   npar = n_elements(parinfo)
   line_fit_all = dblarr(npar) & perror_all = dblarr(npar)
   linecomp_all = (linelist.compname)[uniq(linelist.compname)]
   linename = linelist.linename
   parindx = lonarr(2,nline_all)  ; this only index to lines
   parfitflag = lonarr(npar)
   parlinename = replicate('',npar)
   line_redchi2 = dblarr(npar) & line_status = lonarr(npar)
   line_npix_fit = lonarr(npar)
   inext = 0
   ; assign initial values and limitations
   nline_p=0L ; # of lines in the previous complex
   for icomp=0L, ncomp_all - 1 do begin
      ; find all lines in each line complex
      ind1 = where(linelist.compname eq linecomp_all[icomp], nline)
      for i=0L, nline - 1L do begin
         ngauss = linelist[ind1[i]].ngauss
         parindx[*,i+nline_p] = [inext, inext+3*ngauss-1] 

         for j=0L, ngauss - 1L do begin

            ; add a small, random offset of each broad gaussian component
            if ngauss gt 1 then $
            lnlamini = alog(linelist[ind1[i]].lambda) + (2*randomu(seed,1)-1.)*1d-3 else $
            lnlamini = alog(linelist[ind1[i]].lambda)
            parlinename[inext:inext+2] = linelist[ind1[i]].linename
            parinfo[inext:inext+2].value = [linelist[ind1[i]].fvalue, lnlamini $
                  , linelist[ind1[i]].inisig]
            ; restrict the Gaussian to be positive or negative depending on the initial flux guess
            if parinfo[inext].value gt 0 then begin
              parinfo[inext].limited = [1,0] & parinfo[inext].limits=[0,1d10]
            endif else begin
              parinfo[inext].limited = [0,1] & parinfo[inext].limits=[-1d10,0]
            endelse
            parinfo[inext+1].limited = [1,1]
            parinfo[inext+1].limits=[lnlamini-linelist[ind1[i]].voff, lnlamini+linelist[ind1[i]].voff]
            parinfo[inext+2].limited = [1,1]
            parinfo[inext+2].limits = [linelist[ind1[i]].minsig, linelist[ind1[i]].maxsig]

            inext = inext + 3
        endfor
      endfor
      ; assign the PL parameters
      parinfo[inext:inext+1].value=[1., -0.5]

      ; proceed to the next line complex by skipping the 2 elements used for the local power-law continuum
      inext = inext + 2
      nline_p = nline_p + nline
   endfor

   minwave = min(wave,max=maxwave)
   ; Find which lines are covered
   ind1 = where(linelist.lambda ge minwave and linelist.lambda le maxwave, $
         complement = indd1)
   if indd1[0] ne -1 then linelist[indd1].fitflag = 0
   ind_fit = where(linelist.fitflag eq 1)
   if ind_fit[0] eq -1 then $
     splog, 'No line to fit, return' $
   else begin

      ; loop over each line complex
      uniq_linecomp = uniq(linelist[ind_fit].compname)
      ncomp = n_elements(uniq_linecomp)
      compname = (linelist[ind_fit].compname)[uniq_linecomp]

      nline_p=0L ; # of lines in the previous complex
      for ii=0L, ncomp - 1L do begin
          ind_line = where(linelist.compname eq compname[ii] and $
             linelist.fitflag eq 1, nline_fit)
          linelist_fit = linelist[ind_line]
          ngauss_fit = linelist_fit.ngauss
          if not keyword_set(silent) then splog, 'Fit line complex: ', compname[ii]
          iline_start = ind_line[0] & iline_end = ind_line[nline_fit-1]
          parindx_sub = parindx[0,iline_start]+indgen(3*linelist[iline_start].ngauss)
          for iitmp=1, nline_fit - 1L do begin
              parindx_sub = [parindx_sub, $
               parindx[0,ind_line[iitmp]]+indgen(3*linelist[ind_line[iitmp]].ngauss)]
          endfor
          ; add the PL parameters to the line complex
          parindx_sub = [parindx_sub, max(parindx_sub)+[1,2]]
          parinfo1 = parinfo[ parindx_sub ]

          ; Now setup constraints of parameters
          ; Tie velocity
          arr = linelist_fit.vindex
          arr = arr[sort(arr)]
          arr = arr[uniq(arr)]
          for iuniq=0L, n_elements(arr) - 1L do begin
             ind_fix = where(linelist_fit.vindex eq arr[iuniq],ntmp)
             if ntmp gt 0 and arr[iuniq] ne 0 then begin

                ; For a broad line, enforce symmetric profile if required
                if ntmp eq 1 and linelist_fit[ind_fix[0]].ngauss gt 1 then begin
                   ind_fix = ind_fix[0]
                   pstr0='P(' + string(1+3*( total(ngauss_fit[0:ind_fix]) $
                    - ngauss_fit[ind_fix]),format='(i0)')+')'
                   for igauss=1, ngauss_fit[ind_fix]-1L do begin
                     parinfo1[1+3*(total(ngauss_fit[0:ind_fix]) $
                         - ngauss_fit[ind_fix])+igauss*3].tied = pstr0
                   endfor
                endif

                pstr0='(P(' + string(1+3*( total(ngauss_fit[0:ind_fix[0]])-1) $
                          , format='(i0)')+')'
                for iitmp=1,ntmp-1 do begin
                   parinfo1[1+3*(total(ngauss_fit[0:ind_fix[iitmp]])-1)].tied=pstr0 $
                    + ' - alog('+string(linelist_fit[ind_fix[0]].lambda,format='(f7.2)') $
                    + ')) + alog('+string(linelist_fit[ind_fix[iitmp]].lambda,format='(f7.2)') $
                    + ')'
                endfor
             endif
          endfor
          ; Tie line width
          arr = linelist_fit.windex
          arr = arr[sort(arr)]
          arr = arr[uniq(arr)]
          for iuniq=0L, n_elements(arr) - 1L do begin
             ind_fix = where(linelist_fit.windex eq arr[iuniq],ntmp)
             if ntmp gt 1 and arr[iuniq] ne 0 then begin
                pstr0='P('+string(2+3*(total(ngauss_fit[0:ind_fix[0]])-1) $
                 ,format='(i0)')+')'
                for iitmp=1,ntmp-1 do begin
                   parinfo1[2+3*( total(ngauss_fit[0:ind_fix[iitmp]])-1)].tied = pstr0
                endfor
             endif
          endfor
          ; Tie line flux
          arr = linelist_fit.findex
          arr = arr[sort(arr)]
          arr = arr[uniq(arr)]
          for iuniq=0L, n_elements(arr) - 1L do begin
             ind_fix = where(linelist_fit.findex eq arr[iuniq],ntmp)
             if ntmp gt 1 and arr[iuniq] ne 0 then begin

                pstr0='P('+string(3*( total(ngauss_fit[0:ind_fix[0]])-1) $
                 ,format='(i0)')+')'
                for iitmp=1,ntmp-1 do begin
                   parinfo1[3*(total(ngauss_fit[0:ind_fix[iitmp]])-1)].tied=pstr0 $
                    + '*' + string(linelist_fit[ind_fix[0]].lambda,format='(f7.2)') $
                    + '/' + string(linelist_fit[ind_fix[iitmp]].lambda,format='(f7.2)') $
                    + '*' + string(linelist_fit[ind_fix[iitmp]].fvalue/linelist_fit[ind_fix[0]].fvalue $
                    , format='(f3.1)')
                endfor
             endif
          endfor

          ; message, 'stop'
          ; Do the fit
          ind = where(wave ge linelist_fit[0].minwav $
                  and wave le linelist_fit[0].maxwav $
                  and fitmask ne 0., nnn)
          if nnn gt 10 then begin
             line_fit = mpfitfun('local_conti_line',alog(wave[ind]),line_flux[ind],err[ind] $
                  , MAXITER = 500L, parinfo = parinfo1, xtol=xtol, ftol=ftol $
                  , perror = perror, yfit = yfit2, /quiet $
                  , nfev = nfev, niter = niter, status = status, bestnorm = chi2)
             ttt = where(parinfo1.fixed eq 1 or strlen(parinfo1.tied) gt 0, nfix)
             dof = n_elements(ind) - (n_elements(parinfo1) - nfix)
             redchi2 = chi2/dof
             if not keyword_set(silent) then splog,'LINE-FITTING: MPFIT nfev=',nfev $
                 , ' niter=', niter, ' status=', status, ' redchi2=', redchi2

             ; Do a further step to remove 3sigma absorption below the first fit
             ; Do not do this for the Halpha or Hbeta line complex
             ; first setup fitmask2 and the new line fitting windows
             if keyword_set(auto_abs_rej) and redchi2 lt 100. $
                and compname[ii] ne 'Halpha' and compname[ii] ne 'Hbeta' and compname[ii] ne 'CaII3934' $
             then begin

                n_iter = 1L

       iter1:
                fitmask2 =  replicate(1L, n_elements(wave))
                ind_abs = where(line_flux[ind] lt yfit2 - 3.*err[ind])
                if ind_abs[0] ne -1 then fitmask2[ind[ind_abs]] = 0L
                fitmask2 = fitmask*fitmask2

                ind2 = where(wave ge linelist_fit[0].minwav $
                         and wave le linelist_fit[0].maxwav $
                         and fitmask2 ne 0.)

                if (n_elements(ind2) gt 10L) and (n_elements(ind2) lt n_elements(ind)) then begin

                   line_fit2 = mpfitfun('local_conti_line', alog(wave[ind2]) $
                        , line_flux[ind2], err[ind2] $
                        , MAXITER = 500L, parinfo = parinfo1, xtol=xtol,ftol=ftol $
                        , perror = perror2, yfit = yfit2_2, /quiet $
                        , nfev = nfev, niter = niter, status = status2, bestnorm = chi2_2)
                   ttt = where(parinfo1.fixed eq 1 or strlen(parinfo1.tied) gt 0, nfix)
                   dof2 = n_elements(ind2) - (n_elements(parinfo1) - nfix)
                   redchi2_2 = chi2_2/dof2
                   if not keyword_set(silent) then splog $
                     , 'LINE-FITTING: MPFIT nfev=', nfev $
                     , ' niter=', niter, ' status=', status2, ' redchi2=', redchi2_2 $
                     , ' iter=', n_iter
                   ; replace the original fit if justified
                   if status2 gt 0 and redchi2_2 gt 0 and redchi2_2 lt redchi2 then begin
                      if not keyword_set(silent) then $
                       splog, 'Improved fits by rejecting absorptions. '
                      ind = ind2 & line_fit = line_fit2 & perror = perror2
                      redchi2 = redchi2_2 & status = status2 & yfit2 = yfit2_2
                      fitmask = fitmask2
                      if n_iter lt rej_iter then begin
                         n_iter = n_iter + 1L
                         goto, iter1
                      endif
                   endif

                endif
             endif

             if n_elements(ind_line_fit) eq 0 then ind_line_fit=ind else $
                ind_line_fit = [ind_line_fit, ind]
             line_fit_all[ parindx_sub ] = line_fit
             perror_all[ parindx_sub ] = perror
             parfitflag[parindx_sub] = 1L
             line_redchi2[parindx_sub] = redchi2
             line_status[parindx_sub] = status
             line_npix_fit[parindx_sub] = n_elements(ind)

         endif
         ind_tmp = where(linelist.compname eq compname[ii], nline_n)
         nline_p = nline_p + nline_n

      endfor
      output.fitmask = fitmask
      temp = create_struct('linename',parlinename, 'fitflag', parfitflag $
        , 'ind_line_fit', n_elements(ind_line_fit) gt 0 ? ind_line_fit : [-1L] $
        , 'line_fit', line_fit_all $
        , 'line_fit_err', perror_all, 'line_redchi2', line_redchi2 $
        , 'line_status', line_status, 'line_npix_fit', line_npix_fit)
      output = struct_addtags(output, temp)
   endelse

   ; --------------------------- End of fitting emission lines ---------------

   if keyword_set(diet) then begin
      tag_rej=['WAVE','FLUX','ERR','FITMASK','IND_LINE_FIT']
      remove_tags, output, tag_rej, output1
      output = output1
   endif
   para = output

   ; Write the fit results
   if n_elements(outdir) eq 0 then cd, current=outdir
   if n_elements(output_name) eq 0 then output_name = 'qsofit'
   if keyword_Set(fits) then begin ; output the fits result
      fitsfile = outdir + '/fits/' + output_name + '.fits'
      mwrfits, output, fitsfile, /create
   endif

   if keyword_set(noplot) then return

   ; ------------------------  Plotting block ----------------------------
      ; now make a plot
   if keyword_set(psplot) then begin
      if not keyword_set(append) then begin
         figfile = outdir + '/QA/' + output_name + '.ps'
         ;begplot, name=figfile, xsize = 40, ysize = 25,/color
         begplot, name=figfile, /landscape, /color
      endif
      charsize = 1. & thick = 4. & xticks = 2L & xminor = 5L
      linethick = 0.1 & symsize = 3.
      ang = string(197B)
   endif else begin
      linethick = 1. & symsize = 3.
      xticks = 2L & xminor = 5L
      ang = textoidl('\AA')
   endelse

   uniq_linecomp = uniq(linelist.compname)
   ncomp = n_elements(uniq_linecomp)
   compname = (linelist.compname)[uniq_linecomp]
   ; determine how many line complex are fitted and should be plotted
   plotkey=lonarr(ncomp)
   for ii=0L, ncomp - 1L do begin
      ind_line = where(linelist.compname eq compname[ii] and linelist.fitflag eq 1, $
                nline_fit)
      if nline_fit gt 0 then begin
         linelist_fit = linelist[ind_line]
         pxrange=[linelist_fit[0].minwav,linelist_fit[0].maxwav]
         ind = where(wave ge pxrange[0] and wave le pxrange[1], nnn)
         if nnn gt 10 then plotkey[ii] = 1L ; plot this line complex
      endif
   endfor
   temp=where(plotkey eq 1, nplot)
   ; determine layout
   if nplot gt 0 then plot_layout, nplot, xypos=xypos, omargin=[0.08, 0.05, 0.98, 0.94], pmargin=[0.1,0.1]

   len = 0.04
   if not keyword_set(nsmooth) then nsmooth = 1L
   ytitle = textoidl('Flux Density f_\lambda [10^{-17} erg s^{-1} cm^{-2} ') $
      + ang + textoidl('^{-1}]')
   if ~keyword_set(plot_vel) then xtitle = textoidl('Rest Wavelength [')+ang+']' $
     else xtitle = textoidl('Velocity [km s^{-1}]')


   ; now plot each fitted line complex
   ; first refresh the plot
   plot, [0],[0], /nodata, xsty=5, ysty=5  

   iuse=0 ; index of panel to plot
   for ii=0L, ncomp - 1L do begin

       ind_line = where(linelist.compname eq compname[ii] and linelist.fitflag eq 1, $
                nline_fit)

       if nline_fit gt 0 then begin

          pos=xypos[*, iuse]
          iuse=iuse + 1
          linelist_fit = linelist[ind_line]
          ; splog, 'Plotting line complex: ', compname[ii]
          iline_start = ind_line[0] & iline_end = ind_line[nline_fit-1]
          parindx_sub = parindx[0,iline_start] + indgen(3*linelist[iline_start].ngauss)
          for iitmp=1, nline_fit - 1L do begin
             parindx_sub = [parindx_sub, parindx[0,ind_line[iitmp]] $
              + indgen(3*linelist[ind_line[iitmp]].ngauss)]
          endfor
          parindx_sub_lineonly = parindx_sub
          ; add the PL parameters to the line complex
          parindx_sub = [parindx_sub, max(parindx_sub)+[1,2]]

          pxrange=[linelist_fit[0].minwav,linelist_fit[0].maxwav]
          ind = where(wave ge pxrange[0] and wave le pxrange[1], nnn)
          ind_bad=where(wave ge pxrange[0] and wave le pxrange[1] and fitmask eq 0, nbad)
          ind_good=where(wave ge pxrange[0] and wave le pxrange[1] and fitmask ne 0, ngood)
          if keyword_set(plot_vel) then pxrange=[-3000.,3000.]

          if nnn gt 10 then begin
             ;yrange=[-0.5,1.05]*max((smooth(line_flux[where(flux*sqrt(ivar) gt 5.)],20))[ind])
             yrange=[min(median(line_flux[ind_good],5)), max(median(line_flux[ind_good],5))*1.05]
             if ~keyword_set(plot_vel) then xarr=wave else $
                xarr=(wave/linelist_fit[0].lambda - 1.)*cs
             plot, xarr[ind], smooth(line_flux[ind],nsmooth), /noerase $
              , xrange=pxrange, yrange=yrange $
              , pos = pos,/xsty,thick=linethick, charsize=charsize,xticklen=len $
              , yticklen=len, xthick=thick, ythick=thick, charthick=thick $
              , xticks = xticks, xminor = xminor, psym=10
             oplot, xarr[ind],err[ind],color=cgcolor('gray'),thick=linethick
             if nbad gt 0 then oplot, [xarr[ind_bad]], [line_flux[ind_bad]] $
               , psym=4, thick =thick, color=cgcolor('cyan'),symsize=0.5
              
              ; plot the model  
              pp = line_fit_all[ parindx_sub ]
                ;pp_broad = pp[0:3*linelist_fit[0].ngauss-1]
                ;if n_elements(pp) gt n_elements(pp_broad) then begin
                ;   pp_narrow = pp[3*linelist_fit[0].ngauss:*]
                ;   for jj=1L, n_elements(pp_narrow)/3L do begin
                ;      pp_narrow_1gauss = pp_narrow[(jj-1)*3L:jj*3L-1]
                ;      oplot, wave[ind], onegauss(alog(wave[ind]), pp_narrow_1gauss) $
                ;        , color=cgcolor('cyan')
                ;   endfor
                ;endif
                ;oplot,wave[ind], manygauss(alog(wave[ind]),pp_broad) $
                ; , color=cgcolor('green')
                ;; now plot the total model
               oplot,xarr[ind],local_conti_line(alog(wave[ind]),pp) $
                  , color=cgcolor('red')

               xyouts, pos[2]-0.05,pos[3]-0.03, compname[ii], /norm,charsize=charsize
               ; plot floating linenames
               for jline=0,n_elements(linelist_fit)-1 do begin
                 if ~keyword_set(plot_vel) then $
                 xyouts, linelist_fit[jline].lambda, yrange[1]*(0.9-(jline mod 2)*0.1), $
                    linelist_fit[jline].linename,color=cgcolor('magenta'),charsize=0.7,align=0.5 $
                 else xyouts, (linelist_fit[jline].lambda/linelist_fit[0].lambda - 1.)*cs, $
                   yrange[1]*(0.9-(jline mod 2)*0.1), $
                    linelist_fit[jline].linename,color=cgcolor('magenta'),charsize=0.7,align=0.5
               endfor

          endif
      endif
   endfor

   xyouts,0.5,0.01,xtitle,/norm,charsize=charsize,charthick=thick,align=0.5
   xyouts,0.02,0.5,ytitle,/norm,charsize=charsize,charthick=thick,align=0.5,orien=90

   if keyword_set(objtag) then $
     xyouts, 0.15, 0.96, objtag, /norm, charsize = charsize, charthick=thick
   if keyword_set(SDSS_name) then xyouts, 0.5, 0.96, /norm, SDSS_name $
     , charsize = charsize, charthick = thick
   xyouts, 0.85, 0.96, 'z='+string(z,format='(f5.3)'), /norm $
        , charsize=charsize, charthick=thick, color=cgcolor('red')

   if keyword_set(psplot) and not keyword_set(append) then endplot

end
