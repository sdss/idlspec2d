;+
; NAME:
;   rm_compile_qsofit
;
; PURPOSE:
;   Compile the fitting results from rm_qsofit, including Monte Carlo errors
;
; CALLING SEQUENCE:
;   rm_compile_qsofit,7338,56660,[fiber=], [redshift=], [outdir=]
;   rm_compile_qsofit,[7338,7338],[56660,56664]
;
; INPUTS:
;   plate_in       - Input plate
;   mjd_in         - Input MJD
; OPTIONAL INPUTS:
;   outdir_in      - Path of the fits file and the output file
;   outfile        - Name of the output file
;   linename       - Names of the lines to be compiled; [NLINE]
;                    linename should match those in $IDLRM/etc/qsoline.par,
;                    or whatever line list file is used in the fit
;                    Default: ['Halpha','Hbeta','MgII','CIII','CIV']
;                    or: ['Halpha_br', 'Hbeta_br']
;                    or: ['OIII5007'] ; core+wing [OIII] components
;                    or: ['OIII5007w', 'OIII5007c']
;   contiwave      - Restwave of continuum luminosity 
;                    Default: [1350.,1700.,3000.,5100.]
;
; OUTPUTS:
;   result         - Structure containing the compiled properties
;     LOGLxxxx     - Continuum luminosity, in log (erg/s)
;     LOGLxxxx_ERR - Errors in logL_conti
;     LINENAME     - [peak wave, FWHM (km/s), logL_line (erg/s), rest EW, centroid of top 50% flux] 
;     LINENAME_ERR - Errors in LINENAME
; 
; OPTIONAL OUTPUTS:
;
; INTERNAL ROUTINES CALLED:
;   get_qso_prop
;
; EXTERNAL ROUTINES CALLED:
;   red
;
; COMMENTS:
;
; EXAMPLES:
;
; REVISION HISTORY:
;---------------------------------------------------------------------------------------
pro get_qso_prop, para, linelist=linelist,contiwave=contiwave, $
    z=z, conti_prop=conti_prop, line_prop=line_prop, localfit=localfit

   ; The following feature is obsolete, as the routine will recognize 
   ; whether or not there is a global "conti_fit" tag in para. 
   ; if /localfit, then use different parameter extraction
   ; since the fits in rm_local_qsofit are done differently. 
   ; actually I found this doens't make a difference, as the 
   ; local-conti parameters in the para structure are not populated by
   ; a specific LINENAME. 

   forward_function dluminosity, f_conti_only

   ; Get continuum properties
   if not keyword_set(contiwave) then contiwave=[1350.,3000.,5100.]
   conti_lum_str = 'logL'+string(contiwave,format='(i0)')
   conti_prop = create_struct('null',0.)
   for i=0,n_elements(conti_lum_str)-1 do $
    conti_prop = struct_addtags(conti_prop,create_struct(conti_lum_str[i],0.D))
   ; now add FeII para and conti_fit (default is 14 elements in conti_fit, i.e., 3 poly terms)
   conti_prop = struct_addtags(conti_prop,{conti_fit:dblarr(14), FeII_uv:dblarr(3), FeII_opt:dblarr(3), REW_Fe_4434_4684:0.D, REW_Fe_2250_2650:0.D})
   remove_tags, conti_prop, 'null', conti_prop1
   conti_prop = conti_prop1
   
   ; Get line properties [peak wave, FWHM, flux, EW]
   if not keyword_Set(linelist) then linelist0 = ['Hbeta_br','CIII'] $
   else linelist0 = linelist
   nline = n_elements(linelist0)
   line_prop = create_struct('null',0.)
   for i=0, nline -1 do $
    line_prop = struct_addtags(line_prop,create_struct(linelist0[i],dblarr(5)))
   remove_tags, line_prop, 'null', line_prop1
   line_prop = line_prop1

   cs = 2.9979246d5  ; speed of light, km/s
   dlum = dluminosity(z[0],/cm)
   lumscale = 1d-17*(1.+z[0])*dlum^2*4.*!PI

   ; Assign to the conti_prop struct
   if tag_exist(para, 'CONTI_FIT') then begin
      nconti = n_elements(contiwave)
      conti_fit = para.conti_fit
      tmp_conti = f_conti_only(contiwave, conti_fit[6:*])
      ; Do not extrapolate the continuum fit beyond the BOSS spectral range [3600, 10400]
      ind=where(contiwave*(1.+z[0]) lt 3600. or contiwave*(1.+z[0]) gt 10400.)
      if ind[0] ne -1 then tmp_conti[ind] = 0.D
   
      tmp_conti = alog10(tmp_conti*contiwave*lumscale > 1.) ; ??
      for i=0L, nconti - 1L do conti_prop.(i) = tmp_conti[i]
      ; add FeII fits and conti_fit
      conti_prop.(nconti) = conti_fit
      conti_prop.(nconti+1) = conti_fit[0:2]  ; uv FeII
      conti_prop.(nconti+2) = conti_fit[3:5]  ; opt FeII
      conti_prop.(nconti+3) = get_fe_ew([4434., 4684.],conti_fit,/opt)
      conti_prop.(nconti+4) = get_fe_ew([2250., 2650.],conti_fit,/uv)
   endif 

   ; Assign to the line_prop struct
   if tag_exist(para, 'LINE_FIT') then begin
      line_fit = para.line_fit & line_redchi2 = para.line_redchi2
      fitflag = para.fitflag & linename = strtrim(para.linename)   
      for i=0L, nline - 1L do begin
   
         if strupcase(linelist0[i]) eq 'CIII' then $  ; default is to return CIII+SiIII+AlIII
          ind=where( (strmatch(strupcase(linename),'CIII*') eq 1 or  $
                    strmatch(strupcase(linename),'SIIII*') eq 1 or $
                    strmatch(strupcase(linename),'ALIII*') eq 1) AND fitflag eq 1) $
         else ind=where( strmatch(strupcase(linename), strupcase(linelist0[i]+'*') ) eq 1 $
		 AND fitflag eq 1)
       
         if ind[0] ne -1 then begin
            tmp_line = dblarr(5)
            pp = line_fit[ind]
            tmp = get_multi_gaussian_prop(pp,/diet)
            tmp_line[0:2]= [ exp(tmp[0]), $ ; peak wavelength
                             tmp[1]*cs, $  ; FWHM
                             tmp[2]] ; line flux, unscaled
         
            ; get REW of the line 
            if tag_exist(para, 'CONTI_FIT') then begin 
               tmp_line[3]=tmp_line[2] $
                / f_conti_only(tmp_line[0], conti_fit[6:*]) ; using the global fit
            endif else begin  ; try to see if it is a local power-law conti fit
               conti_fit = line_fit[max(ind) + [1,2] ]
               f_conti = conti_fit[0]*(tmp_line[0]/3000.0)^conti_fit[1]
               if f_conti gt 0 then $
                  tmp_line[3]=tmp_line[2] / f_conti
            endelse
            
            ; add the centroid computed using >0.5*peakflux
            tmp_line[4] = exp(tmp[5])
 
            ; Scale line flux to line luminosity
            tmp_line[2]=alog10(tmp_line[2]*lumscale > 1.)
         
            line_prop.(i) = tmp_line             
         endif
   
      endfor
   endif

end
;---------------------------------------------------------------------------------------
pro rm_compile_qsofit,plate_in,mjd_in,input_fiber=input_fiber, tags=tags, redshift=redshift, $
    outdir_in=outdir_in,outfile=outfile,linename=linename,contiwave=contiwave,result=result, $
    rm_plate=rm_plate,add_bhmass=add_bhmass,rm_ID=rm_ID, plate_subdir=plate_subdir, dr7=dr7

;;; If supply TAGS (list of objects) instead, set PLATE_IN, MJD_IN and input_fiber=lonarr(n)

   !except=0

   if n_elements(rm_plate) eq 0 then rm_plate=1 ; default is to compile all objects on a RM plate

   if keyword_set(rm_plate) then begin
      if (N_elements(plate_in) GT 1) then begin
         for i=0, N_elements(plate_in)-1 do $
          rm_compile_qsofit,plate_in[i],mjd_in[i],linename=linename,contiwave=contiwave
         return
      endif
   endif

   if not keyword_set(outdir_in) then outdir = $
     getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
       + '/' + string(plate_in[0],format='(i4.4)') + '/' $
       + 'qsofit/' else outdir = outdir_in

   if n_elements(input_fiber) eq 0 then begin
      ; read all fibers, using the order specified by fibermap 
      target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
      fibermap = mrdfits(target_file,1,/silent)
      fiber_all = fibermap.fiberid
      ind = where((fibermap.plate)[*,0] eq plate_in[0] and (fibermap.mjd)[*,0] eq mjd_in[0])
      fiber = reform(fiber_all[ind, *])
      if plate_in[0] eq 0 then fiber = 1L+lindgen(1000L)  ; for the coadded plate
   endif else fiber=input_fiber

   ; Keep the QSO targets only [0-848] if this is for the RM plates
   if keyword_set(rm_plate) then fiber=fiber[0:848]

   nnn=n_elements(fiber)
   if n_elements(plate_in) eq nnn then plate=plate_in else plate=replicate(plate_in[0], nnn)
   if n_elements(mjd_in) eq nnn then mjd=mjd_in else mjd=replicate(mjd_in[0], nnn)
   platestr=string(plate,format='(i4.4)')
   if ~keyword_set(dr7) then fiberstr=string(fiber,format='(i4.4)') else fiberstr=string(fiber,format='(i3.3)')
   mjdstr=string(mjd,format='(i5.5)')

   ; get file identification
   if ~keyword_set(tags) then begin
     if ~keyword_set(dr7) then tags = platestr+'-'+mjdstr+'-'+fiberstr else $
        tags = platestr+'-'+fiberstr+'-'+mjdstr
   endif
   
   ; Create the output structure
   result = {tag:'',rm_ID:-1L,plate:0L,mjd:0L,fiberid:0L,z:0.D,conti_redchi2:0.D}
   nheader = n_tags(result)  
 
   if not keyword_Set(contiwave) then contiwave=[1350.,1700.,3000.,5100.]
   nconti = n_elements(contiwave)
   conti_lum_str = 'logL'+string(contiwave,format='(i0)')
   for i=0,nconti-1 do $
    result = struct_addtags(result,create_struct(conti_lum_str[i],0.D))
   ; Add FeII fits and conti_fit
   result=struct_addtags(result,{conti_fit:dblarr(14), FeII_uv:dblarr(3), FeII_opt:dblarr(3),REW_Fe_4434_4684:0.D,REW_Fe_2250_2650:0.D})

   conti_err_str = 'logL'+string(contiwave,format='(i0)') + '_err'
   result_conti_err = create_struct(conti_err_str[0],-1.D)
   for i=1,nconti-1 do $
    result_conti_err=struct_addtags(result_conti_err,create_struct(conti_err_str[i],-1.D))
   ; Add FeII fits err
   result_conti_err=struct_addtags(result_conti_err,{conti_fit_err:dblarr(14) - 1.D, FeII_uv_err:dblarr(3) - 1.D, $
      FeII_opt_err:dblarr(3) - 1.D, REW_Fe_4434_4684_Err:-1.D, REW_Fe_2250_2650_Err:-1.D})

   ; Find what lines were fit
   if not keyword_Set(linename) then begin
      ;file1 = outdir + 'fits/' + tags[0] + '.fits'
      ;tmp=mrdfits(file1,1,/silent)
      ;linename=strtrim(tmp.linename)
      ;linename=linename[uniq(linename)]
      linename = ['SII6718','Halpha', 'Halpha_br', 'Hbeta', 'Hbeta_br', 'HeII4687','HeII4687_BR', 'OIII5007', $
                  'OIII5007c','MgII', 'MgII_br', 'CIII','CIII_br', 'SiIII1892','AlIII1857', 'NIII1750','CIV', 'CIV_br', $
                  'HeII1640','HeII1640_br', 'SiIV_OIV', 'OI1304', 'Lya', 'NV1240']
   endif
   nline = n_elements(linename)
   for i=0, nline -1 do $
    result = struct_addtags(result,create_struct(linename[i],dblarr(5)))
   for i=0L, nline-1 do $
    result = struct_addtags(result,create_struct(linename[i]+'_redchi2',0.D))

   line_err_str = linename + '_err'
   result_line_err=create_struct(line_err_str[0], replicate(-1.D,5))
   for i=1, nline -1 do $
    result_line_err = $
      struct_addtags(result_line_err,create_struct(line_err_str[i],replicate(-1.D,5)))  
   
   result = replicate(result, nnn)
   struct0 = result[0]
   if n_elements(rm_ID) eq nnn then result.rm_ID=rm_ID else result.rm_ID=lindgen(nnn)
   result.plate=plate
   result.fiberid=fiber
   result.mjd=mjd
   if n_elements(redshift) eq nnn then result.z = redshift
   result.tag = tags  
 
   result_conti_err = replicate(result_conti_err,nnn)
   result_line_err = replicate(result_line_err,nnn)
   
   red, omegalambda=0.7,omega0=0.3,h100=0.7
   
   ; These are the original tag names in the fits
   for i=0L, nnn - 1L do begin
   
      ; Get the fitting results on the original spectrum
      if ~keyword_set(plate_subdir) then fitsfile = outdir + 'fits/' + tags[i] + '.fits*' $
         else fitsfile = outdir + 'fits/' + platestr[i] + '/' + tags[i] + '.fits*'
 
      if file_test(fitsfile) eq 1 then begin
      
         para = mrdfits(fitsfile, 1,/silent)
         z0 = para.z
         if not keyword_set(redshift) then result[i].z = z0 

         ;------------------- 
         ; populate the redchi2 for continuum and each line complex fit
         if tag_exist(para, 'conti_redchi2') then result[i].conti_redchi2=para.conti_redchi2
         ; Does the line fit exist?
         if tag_exist(para, 'linename') eq 1 then begin
           for iline=0L, nline-1 do begin
              ind_line=where( strmatch(strupcase(para.linename), strupcase(linename[iline]+'*') ) ) 
              ind_tag = where(strtrim(tag_names(result)) eq strupcase(linename[iline]+'_redchi2') )
              result[i].(ind_tag) = max( para.line_redchi2[ind_line] )
           endfor
         endif
         ;-------------------

         get_qso_prop, para,linelist=linename,contiwave=contiwave, $
          z=z0, conti_prop=conti_prop, line_prop=line_prop
         
         str_tmp = result[i]
         copy_struct, conti_prop, str_tmp
         copy_struct, line_prop, str_tmp
         result[i] = str_tmp
      
      endif
   
      ; Get errors
      if ~keyword_set(plate_subdir) then errfile = outdir + 'err/' + tags[i] + '.fits*' else $
          errfile = outdir + 'err/' + platestr[i] + '/' + tags[i] + '.fits*'
      if file_test(errfile) eq 1 then begin
         para_err = mrdfits(errfile,1,/silent)
         ntrial = n_elements(para_err)
         errarray = replicate(struct0, ntrial)
         for jj=0L, ntrial - 1L do begin
            para = para_err[jj]
            get_qso_prop, para,linelist=linename,contiwave=contiwave, $
             z=z0, conti_prop=conti_prop, line_prop=line_prop
             
            str_tmp = struct0 
            copy_struct, conti_prop, str_tmp
            copy_struct, line_prop, str_tmp 
            errarray[jj] = str_tmp
         endfor
         
         ; Now estimate the error using the 68% percentile
         ; Do continuum
         for jj=0L, nconti-1 do begin
            arr = errarray.(jj+nheader)
            ; Is this continuum covered?
            if (where(arr ne 0))[0] ne -1 then begin
               sigma = 0.5*( quantile_2d(0.84,arr,dim=1) - quantile_2d(0.16,arr,dim=1) )
               result_conti_err[i].(jj) = sigma
            endif
         endfor
         ; Do conti_fit and the uv and opt FeII fits
         for jj=0L, 2 do begin
           arr = errarray.(jj+nheader+nconti)
           sigma = 0.5*( quantile_2d(0.84,arr,dim=2) - quantile_2d(0.16,arr,dim=2) )
           result_conti_err[i].(jj+nconti) = sigma
         endfor
         ; Do the optical FeII REW around the [4434,4684]
         arr = errarray.(nheader+nconti+3)
         sigma = 0.5*( quantile_2d(0.84,arr,dim=1) - quantile_2d(0.16,arr,dim=1) )
         result_conti_err[i].(3+nconti) = sigma
         ; Do the UV FeII REW around the [2250,2650]
         arr = errarray.(nheader+nconti+4)
         sigma = 0.5*( quantile_2d(0.84,arr,dim=1) - quantile_2d(0.16,arr,dim=1) )
         result_conti_err[i].(4+nconti) = sigma

         ; Do the lines
         for jj=0L, nline - 1 do begin
            arr = errarray.(jj+nheader+nconti+5)
            ; Is this line fitted?
            if (where(arr ne 0))[0] ne -1 then begin
               sigma = 0.5*( quantile_2d(0.84,arr,dim=2) - quantile_2d(0.16,arr,dim=2) )
               result_line_err[i].(jj) = sigma
            endif
         endfor
         
      endif
    
      splog, 'Finished: ', i+1, '/', nnn
 
   endfor
   
   result = struct_addtags(result, result_conti_err)
   result = struct_addtags(result, result_line_err)

   ; add SE virial BH mass estimates if required
   if keyword_set(add_bhmass) then begin   
      newstr = replicate( {logBH_CIV_VP06:0.D, logBH_CIV_VP06_err:-1.D, $
                logBH_MgII_S11:0.D, logBH_MgII_S11_err:-1.D, $
                logBH_HB_VP06:0.D, logBH_HB_VP06_err:-1.D }, nnn)
      ;CIV mass
      aa=0.66 & bb=0.53
      ind=where(result.logL1350 gt 0 and (result.CIV_BR)[1,*] gt 0)
      if ind[0] ne -1 then $
        newstr[ind].logBH_CIV_VP06=aa + bb*( result[ind].logL1350 - 44.) $
              + 2.*alog10( (result[ind].CIV_br)[1,*] )
      ind=where(result.logL1350_err gt 0 and $
            (result.CIV_BR_err)[1,*] gt 0 and (result.CIV_BR)[1,*] gt 0 and $
            (result.CIV_BR_err)[3,*] gt 0 and $
            (result.CIV_BR)[3,*]/(result.CIV_BR_err)[3,*] gt 2.  )
      if ind[0] ne -1 then begin 
         newstr[ind].logBH_CIV_VP06_err=sqrt( bb^2*result[ind].logL1350_err^2  $
              + 2.^2*( (result[ind].CIV_BR_err)[1,*]/(result[ind].CIV_BR)[1,*]/alog(10.) )^2 )          
      endif
      ; MgII mass
      aa=0.74 & bb=0.62
      ind=where(result.logL3000 gt 0 and (result.MGII_BR)[1,*] gt 0)
      if ind[0] ne -1 then newstr[ind].logBH_MGII_S11=aa + bb*( result[ind].logL3000 - 44.) $
              + 2.*alog10( (result[ind].MGII_br)[1,*] )
      ind=where(result.logL3000_err gt 0 and $
            (result.MGII_BR_err)[1,*] gt 0 and (result.MGII_BR)[1,*] gt 0 and $
            (result.MGII_BR_err)[3,*] gt 0 and $
            (result.MGII_BR)[3,*]/(result.MGII_BR_err)[3,*] gt 2.  )
      if ind[0] ne -1 then begin 
         newstr[ind].logBH_MGII_S11_err=sqrt( bb^2*result[ind].logL3000_err^2  $
              + 2.^2*( (result[ind].MGII_BR_err)[1,*]/(result[ind].MGII_BR)[1,*]/alog(10.) )^2 )
      endif
      ; Hb mass
      aa=0.91 & bb=0.50
      ind=where(result.logL5100 gt 0 and (result.Hbeta_BR)[1,*] gt 0)
      if ind[0] ne -1 then newstr[ind].logBH_HB_VP06=aa + bb*( result[ind].logL5100 - 44.) $
              + 2.*alog10( (result[ind].Hbeta_br)[1,*] )
      ind=where(result.logL5100_err gt 0 and $
            (result.Hbeta_BR_err)[1,*] gt 0 and (result.Hbeta_BR)[1,*] gt 0 and $
            (result.Hbeta_BR_err)[3,*] gt 0 and $
            (result.Hbeta_BR)[3,*]/(result.Hbeta_BR_err)[3,*] gt 2.  )
      if ind[0] ne -1 then begin
         newstr[ind].logBH_HB_VP06_err=sqrt( bb^2*result[ind].logL5100_err^2  $
              + 2.^2*( (result[ind].Hbeta_BR_err)[1,*]/(result[ind].Hbeta_BR)[1,*]/alog(10.) )^2 )
      endif
      result = struct_addtags(result, newstr)
   endif

   ; message,'stop'
      
   if not keyword_set(outfile) then outfile1=outdir+ 'qso_prop-'+platestr[0]+'-'+mjdstr[0]+'.fits' $
   else outfile1 = outdir + outfile
   
   mwrfits, result, outfile1, /create
      
end
