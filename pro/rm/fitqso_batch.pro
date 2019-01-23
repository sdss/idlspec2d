; Perform rm_qsofit to a list of plate-fiber-mjd

pro fitqso_batch, plate,fiber,mjd,zfit=zfit,ntrial=ntrial, $
     silent=silent, outdir=outdir,tag=tag,ra=ra,dec=dec,$
     emparfile=emparfile, linename=linename,contiwave=contiwave, $
     noqa=noqa,everynqa=everynqa,istart=istart,resumeplot=resumeplot

   if not keyword_set(emparfile) then $
     emparfile='/home/yshen/products/Linux/idlrm/etc/qsoline_qsovar.par'

   if n_elements(noqa) eq 0 then noqa=0L
   if keyword_set(noqa) then psplot=0 else psplot=1

   if ~keyword_set(everynqa) then everynqa=1L ; make QA for every object

   if not keyword_set(tag) then tag = 'qsofit'
   if not keyword_set(outdir) then outdir='/data1/quasar/yshen/DR12_QSO_fits/' 
   ; '/data3/quasar/yshen/work/md_qso/'

   if n_elements(silent) eq 0 then silent = 1L

   nqso = n_elements(plate)

   if file_test(outdir) eq 0 then spawn, 'mkdir '+ outdir
   fitsdir = outdir + 'fits/'
   if file_test(fitsdir) eq 0 then spawn, 'mkdir ' + fitsdir
   qadir = outdir + 'QA/'
   if file_test(qadir) eq 0 then spawn, 'mkdir ' + qadir
   errdir = outdir + 'err/'
   if file_test(errdir) eq 0 then spawn, 'mkdir ' + errdir

   ; Default number of Monte Carlo trials of error estimation
   if n_elements(ntrial) eq 0 then ntrial=50L

   ; setup the plot
   if psplot eq 1 and ~keyword_set(resumeplot) then begin
      figfile=outdir+'QA-'+ tag + '.ps'
      begplot, name=figfile, /landscape, /color
   endif

   ; Setup error handeling
   errfile = outdir + tag + '.err.out.mjd' + string(current_mjd(), format='(f9.3)')
   openw, lun, errfile, /get_lun

   ;message, 'stop'
   if ~keyword_set(istart) then istart=0L
   for i=istart, nqso - 1L do begin

      ;; Estabilish error handlers
      Catch, error_status
      if error_status ne 0 then begin
         ;print, 'Error index: ', error_status
         catch, /cancel
         printf, lun, 'obj plate mjd fiber'
         printf, lun, plate[i], mjd[i], fiber[i]
         printf, lun, 'Error message: ', !error_state.msg
         printf, lun, ' '
         continue
      endif

      ; read in the spectrum
      readspec_sdss, plate[i],fiber[i],mjd[i],wave=lam0,flux=flux0,invvar=ivar0, $
        plug_ra=plug_ra, plug_dec=plug_dec

      err0=ivar0*0
      ind=where(ivar0 ne 0)
      if ind[0] ne -1 then err0[ind]=1./sqrt(ivar0[ind])

      ; Deredden the spectrum
      if ~keyword_set(ra) then begin
        ra_use=plug_ra & dec_use=plug_dec
      endif else begin
        ra_use=ra[i] & dec_use=dec[i]
      endelse
      dereddening_spec, lam0, flux0, err=err0 $
         , ra=ra_use, dec=dec_use $
         , dered_flux=dered_flux, dered_err=dered_err
      flux0 = dered_flux & err0 = dered_err

      objtag1=string(plate[i],format='(i4.4)')+'-'+string(mjd[i],format='(i5.5)') $
         +'-'+string(fiber[i],format='(i4.4)')
      objtag=tag +' '+string(i,format='(i0)')+': ' + objtag1
      radec2string,ra_use,dec_use,dummy,rahr=rahr,ramin=ramin,rasec=rasec $
        , decdeg=decdeg,decmin=decmin,decsec=decsec
      sdss_name=rahr+ramin+rasec+decdeg+decmin+decsec

      z0 = zfit[i]
      ; check whether or not to plot this object
      if (i mod everynqa) eq 0 and psplot then noplot=0 else noplot=1
      rm_qsofit,lam0,flux0,err0,z0,outdir=outdir,/append,objtag=objtag,silent=silent, $
        /fits, psplot=psplot, noplot=noplot, output_name=objtag1, $
        sdss_name=sdss_name,nsmooth=3,emparfile=emparfile
      splog, 'Fitting: ', objtag
      ; proceed to compute error array
      if n_elements(err_array) ne 0 then tmp=temporary(err_array)
      for jj=0,ntrial-1 do begin
         rm_qsofit,lam0,flux0,err0,z0,/noplot,/silent,/diet,/add_noise,para=para,emparfile=emparfile
         tag_rej=['CONTI_FIT_ERR','LINE_FIT_ERR']
         remove_tags, para, tag_rej, para_keep

         splog, '  Error trial: ',jj+1
         ; Only add this trial if the continuum fit is OK
         ; The large redchi2 is to keep fits on very high SN spectrum
         if para_keep.conti_redchi2 lt 100. and para_keep.conti_redchi2 gt 0 then begin
            if n_elements(err_array) eq 0 then err_array=para_keep $
            else err_array=[err_array,para_keep]
         endif
      endfor
      fits_errfile=errdir+objtag1+'.fits'
      if n_elements(err_array) ne 0 then mwrfits, err_array, fits_errfile, /create
   endfor

   if psplot eq 1 then begin
      endplot
      ; gzip the QA file
      spawn, 'gzip ' + figfile
   endif

   ; output the complication of line measurements
   ; set rm_plate=0 to avoid the default setting for the RM plates
   outfile='qso_prop-' + tag + '.fits'
   if not keyword_set(linename) then $
      linename=['Halpha_br', 'Hbeta_br', 'MgII_br', 'CIII_br', 'CIV_br', 'OIII5007']
   if not keyword_set(contiwave) then contiwave=[1350.,1700.,3000.,5100.]
   rm_compile_qsofit,plate,mjd,input_fiber=fiber,outdir_in=outdir,outfile=outfile,rm_plate=0, $
     linename=linename,contiwave=contiwave,/add_bhmass

   ;; close error log
   close, lun
   free_lun, lun

end



