;+
; NAME:
;   rm_fitplate
;
; PURPOSE:
;   Fit the RM quasars on one plate
;
; CALLING SEQUENCE:
;   rm_fitplate,7338,56660,zlist=zlist,calibdir='recalib/'
;   rm_fitplate,[7338,7339],[56660,56683]
;   rm_fitplate,0,56837,calibdir='wh_skysub/',outdir=outdir, $
;    outfile=outfile,emparfile=emparfile
;
; INPUTS:
;   plate	- Plate number
;   mjd		- MJD
;
; OPTIONAL INPUTS:
;   calibdir - 'recalib/' [path of the spPlate files]
;   zlist    - Enforce the fits using the input redshift; same order as in fibermap
;   range    - range of fitted QSO index, e.g., [0,10]
;   ntrial   - number of trials in MC error estimation; default 50
;
;-----------------------------------------------------------------------------
pro rm_fit1fiber,plate,fiberid,mjd,calibdir=calibdir,ntrial=ntrial,silent=silent $
      , ra=ra,dec=dec, zfit=zfit, deredden=deredden, pause=pause $
      , rm_indd=rm_indd, _extra=extra, para=para $
      , outdir=outdir,errdir=errdir,fits=fits,psplot=psplot,append=append $
      , err_scal=err_scal

   if n_elements(calibdir) eq 0 then calibdir = 'wh_skysub/'  ;'recalib/'

   if n_elements(fiberid) gt 1 then begin
      nobj = n_elements(fiberid)
      if n_elements(plate) eq 1 then plate=replicate(plate,nobj)
      if n_elements(plate) ne nobj then message, 'Plate dim does not match fiberid'
      if n_elements(mjd) eq 1 then mjd=replicate(mjd,nobj)
      if n_elements(mjd) ne nobj then message, 'mjd dim does not match fiberid'

      if n_elements(ra) ne nobj then ra=replicate(0., nobj)
      if n_elements(dec) ne nobj then dec=replicate(0., nobj)
      if n_elements(rm_indd) ne nobj then rm_indd=replicate(0,nobj)     
 
      for i=0L, nobj-1L do begin
         rm_fit1fiber,plate[i],fiberid[i],mjd[i],calibdir=calibdir,ntrial=ntrial $
           ,silent=silent, ra=ra[i],dec=dec[i], zfit=zfit[i], deredden=deredden $
           ,rm_indd=rm_indd[i], _extra=extra $
           ,outdir=outdir,errdir=errdir,fits=fits,psplot=psplot,append=append
         if keyword_set(pause) then pause
         if i eq nobj-1 then return
      endfor
   endif

   rm_readspec,plate,fiberid,mjd=mjd,calibdir=calibdir $
      ,zans=zans,wave=lam0,flux=flux0,invvar=ivar0,/silent,_extra=extra
   npix1=n_elements(flux0)
   
   err0=ivar0*0
   ind=where(ivar0 ne 0)
   if ind[0] ne -1 then err0[ind]=1./sqrt(ivar0[ind])

   ; degrade the spectrum if err_scal is set
   if keyword_set(err_scal) then begin ; 
      ivar0 = ivar0/(err_scal^2)
      err0 = err0*err_scal
      flux0 = flux0 + randomn(seed, npix1)*err0
   endif
      
   ; Deredden the spectrum
   if n_elements(deredden) eq 0 then deredden=1L
   if keyword_set(deredden) then begin
      dereddening_spec, lam0, flux0, err=err0 $
         , ra=ra, dec=dec $
         , dered_flux = dered_flux, dered_err = dered_err
      flux0 = dered_flux & err0 = dered_err
   endif

   objtag1=string(plate,format='(i4.4)')+'-'+string(mjd,format='(i5.5)') $
         +'-'+string(fiberid,format='(i4.4)')
   objtag='RM Target '+string(rm_indd,format='(i3.3)')+': '+objtag1
   if n_elements(ra) eq 0 then ra1=0. else ra1=ra
   if n_elements(dec) eq 0 then dec1=0. else dec1=dec
   radec2string,ra1,dec1,dummy,rahr=rahr,ramin=ramin,rasec=rasec $
        , decdeg=decdeg,decmin=decmin,decsec=decsec
   sdss_name=rahr+ramin+rasec+decdeg+decmin+decsec

   if keyword_set(zfit) then z0=zfit[0] else z0=zans.z 
     
   rm_qsofit,lam0,flux0,err0,z0,outdir=outdir,append=append $
        ,objtag=objtag,silent=silent, _extra=extra $
        ,fits=fits, psplot=psplot, output_name=objtag1 $
        ,sdss_name=sdss_name,nsmooth=3, para=para
   splog, 'Fitting: ', objtag

   ; proceed to compute error array
   if n_elements(err_array) ne 0 then tmp=temporary(err_array)
   if n_elements(ntrial) eq 0 then ntrial = 0L
   for jj=0,ntrial-1 do begin
      rm_qsofit,lam0,flux0,err0,z0,/noplot,/silent,/diet,/add_noise,para=para
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
   if keyword_set(errdir) then begin
      fits_errfile=errdir+objtag1+'.fits'
      if n_elements(err_array) ne 0 then mwrfits, err_array, fits_errfile, /create
   endif

end

;-------------------------------------------------------------------------
pro rm_fitplate, plate,mjd,zlist=zlist,calibdir=calibdir,range=range,ntrial=ntrial, $
      topdir=topdir,silent=silent,coadd=coadd,emparfile=emparfile,outdir=outdir, $
      linename=linename,contiwave=contiwave,add_bhmass=add_bhmass,outfile=outfile, $
      rm_plate=rm_plate,figfile=figfile,trim_edge=trim_edge, $
      err_scal=err_scal

   if not keyword_set(emparfile) then $
     emparfile=getenv('IDLRM_DIR')+'/etc/qsoline.par'

   if (size(plate))[0] ne 0 then begin
      for i=0, n_elements(plate)-1 do $
        rm_fitplate, plate[i],mjd[i],zlist=zlist,calibdir=calibdir,range=range,ntrial=ntrial, $
         silent=silent,coadd=coadd,emparfile=emparfile,outdir=outdir, $
         linename=linename,contiwave=contiwave
      return
   endif

   if ~keyword_set(topdir) then begin
       topdir = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
           + '/' + string(plate[0],format='(i4.4)') + '/'
   endif
   if n_elements(calibdir) eq 0 then begin
      calibdir = 'wh_skysub/'
      if plate[0] eq 0 or keyword_set(coadd) then calibdir='wh_skysub/'
   endif

   if keyword_set(calibdir) then specdir=topdir+calibdir
   
   if n_elements(silent) eq 0 then silent = 1L

   if n_elements(trim_edge) eq 0 then trim_edge=1L ; default is to trim the spec within [3650-10300]

   if n_elements(add_bhmass) eq 0 then add_bhmass=1L ; default is to add BH mass in rm_compile_qsofit

   ; read in the master file
   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)

   ; default to use the zfinal in fibermap as input z in rm_qsofit
   if n_elements(zlist) eq 0 then zlist = fibermap.zfinal

   if not keyword_set(outdir) then outdir = topdir + 'qsofit/'
   if file_test(outdir) eq 0 then spawn, 'mkdir '+ outdir
   fitsdir = outdir + 'fits/'
   if file_test(fitsdir) eq 0 then spawn, 'mkdir ' + fitsdir
   qadir = outdir + 'QA/'
   if file_test(qadir) eq 0 then spawn, 'mkdir ' + qadir
   errdir = outdir + 'err/'
   if file_test(errdir) eq 0 then spawn, 'mkdir ' + errdir

   ; Default number of Monte Carlo trials of error estimation
   if n_elements(ntrial) eq 0 then ntrial=50L

   ; --------
   ; proceed to fit the plate
   if plate ne 0 and not keyword_set(coadd) then begin
      platelist = fibermap[0].plate & mjdlist = fibermap[0].mjd
      ind = (where(platelist eq plate[0] and mjdlist eq mjd[0]))[0]
      indd = where(strmatch(fibermap.sourcetype, 'RM*'), nqso)
      if keyword_set(zlist) then zfit = zlist[indd]
      fiberid = (fibermap.fiberid)[ind,indd]
      fibermap = fibermap[indd]
   endif else begin ; If the plate is the coadded plate
      fiberid = lindgen(1000L) + 1L
      indd = where(strmatch(fibermap.sourcetype, 'RM*'), nqso)
      if keyword_set(zlist) then zfit = zlist[indd]
      fiberid = fiberid[indd]
      fibermap = fibermap[indd]
   endelse
   ; trim the list if required
   if keyword_set(range) then begin
      fiberid = fiberid[range[0]:range[1]]
      fibermap = fibermap[range[0]:range[1]]
      if keyword_set(zlist) then zfit = zfit[range[0]:range[1]]
      indd = indd[range[0]:range[1]]
   endif
   nqso = n_elements(indd)

   ; setup the plot
   if not keyword_set(figfile) then figfile1=outdir+'QA-'+string(plate,format='(i4.4)')+'-' $
     +string(mjd,format='(i5.5)') + '.ps' else figfile1=outdir+figfile
   begplot, name=figfile1, /landscape, /color
   
   ; Setup error handeling
   errfile = outdir + 'err.out.mjd' + string(current_mjd(), format='(f9.3)')
   openw, lun, errfile, /get_lun

   ;message, 'stop'
   
   for i=0L, nqso - 1L do begin
      
       ; Estabilish error handlers
      Catch, error_status
      if error_status ne 0 then begin
         ;print, 'Error index: ', error_status
         catch, /cancel
         printf, lun, 'obj plate mjd fiber'
         if i le nqso-1 then printf, lun, plate, mjd, fiberid[i]
         printf, lun, 'Error message: ', !error_state.msg
         printf, lun, ' '
         continue
      endif

      
      rm_readspec,plate[0],fiberid[i],mjd=mjd[0],path=specdir $
      ,zans=zans,wave=lam0,flux=flux0,invvar=ivar0,/silent
      npix1=n_elements(flux0)

      if keyword_set(trim_edge) then begin ; trim the edge of the spectrum, which is typically bad
        ind_trim=where(lam0 gt 10350. or lam0 lt 3650.)
        if ind_trim[0] ne -1 then ivar0[ind_trim]=0
      endif
      
      err0=ivar0*0
      ind=where(ivar0 ne 0)
      if ind[0] ne -1 then err0[ind]=1./sqrt(ivar0[ind])

      ; degrade the spectrum if err_scal is set
      if keyword_set(err_scal) then begin ; 
         ivar0 = ivar0/(err_scal^2)
         err0 = err0*err_scal
         flux0 = flux0 + randomn(seed, npix1)*err0
      endif
      
      ; Deredden the spectrum
      ra=fibermap[i].ra & dec=fibermap[i].dec
      dereddening_spec, lam0, flux0, err=err0 $
         , ra=ra, dec=dec $
         , dered_flux = dered_flux, dered_err = dered_err
      flux0 = dered_flux & err0 = dered_err
      
      objtag1=string(plate[0],format='(i4.4)')+'-'+string(mjd[0],format='(i5.5)') $
         +'-'+string(fiberid[i],format='(i4.4)')
      objtag='RM Target '+string(indd[i],format='(i3.3)')+': '+objtag1
      ;print, objtag, plate[0],fiberid[i],mjd[0]

      radec2string,ra,dec,dummy,rahr=rahr,ramin=ramin,rasec=rasec $
        , decdeg=decdeg,decmin=decmin,decsec=decsec
      sdss_name=rahr+ramin+rasec+decdeg+decmin+decsec
      
      if keyword_set(zlist) then z0 = zfit[i] else z0 = zans.z
      rm_qsofit, lam0,flux0,err0,z0,outdir=outdir,/append,objtag=objtag,silent=silent $
        , /fits, /psplot, output_name=objtag1, sdss_name=sdss_name,nsmooth=3,emparfile=emparfile
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

   endplot
   ; fix the orientation of the plot to seascape
   print, figfile1
   cgfixps, figfile1

   ; gzip the QA file
   spawn, 'gzip -f ' + figfile1

   ; output the complication of line measurements
   rm_compile_qsofit,plate[0],mjd[0], outdir_in=outdir,outfile=outfile, $
     linename=linename,contiwave=contiwave,add_bhmass=add_bhmass,rm_plate=rm_plate

   ; close error log
   close, lun
   free_lun, lun

end
