; routines to make plots for the sample characterization paper


;------ get target properties
pro comp_target_prop

    file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
    target = mrdfits(file,1)

    ; keep only RM QSO targets
    target = target[0:848]
    nobj = n_elements(target)

    ; store the first SDSS confirmation epoch (either DR7 or DR10) in plate_sdss_first
    more_tags = {zsys:0.D, zsys_err:0.D, BOSS_TARGET1:0LL, BOSS_TARGET2:0LL, ANCILLARY_TARGET1:0LL, ANCILLARY_TARGET2:0LL, $
                 PRIMTARGET:0LL, SECTARGET:0LL, OTHER_TARGET:'', $
                 plate_sdss_first:0L, fiberid_sdss_first:0L, mjd_sdss_first:0L, $
                 plate_all:'', fiberid_all:'', mjd_all:'',lambda_eff_all:'', $
                 ALLWISE1234:DBLARR(4)-1, ALLWISE1234_ERR:DBLARR(4)-1, ALLWISE_OFFSET:-1.D,$
                 UNWISE1234:DBLARR(4)-1, UNWISE1234_ERR:DBLARR(4)-1, UNWISE_OFFSET:-1.D, $
                 BAL_FLAG:0L, FIRST_FR_TYPE:-1L, FIRST_Fint_MJY:-1.D, FINT_REST6CM_MJY_OBS:-1.d, $
                 LOGFNU2500A_ERGS_OBS:0.d, R_6CM_2500A:-1.d, $
                 COMMENT:''}
    more_tags = replicate(more_tags, nobj)
    target = struct_addtags(target, more_tags)
    

    ; add BOSS TARGET BITS
    outfile = '/data3/yshen/work/sdssrm_sample_char/dr10_sdssrm.fits'
    if file_test(outfile) ne 1 then begin ; create this file is nonexistant    
      ; file = '/data2/yshen/data/yshen/sdss3/boss/spAll-v5_6_0.fits' ; this was used for the SDSS-RM targeting; the RA/DEC seems incorrect
      file = '/data1/yshen/data/sdss3/boss/v5_7_1_dr12/spAll-v5_7_0.fits'  ; this was DR12. the RA/DES seems incorrect either
      cols0 = ['ra', 'dec', 'plug_ra', 'plug_dec', 'plate']
      cols = ['ra', 'dec', 'plug_ra', 'plug_dec', 'plate', 'fiberid', 'mjd', 'BOSS_TARGET1', 'BOSS_TARGET2', 'ANCILLARY_TARGET1', $
              'ANCILLARY_TARGET2','lambda_eff', 'zoffset']
      dr10 = hogg_mrdfits(file,1, columns = cols0)
      ; the spall files have incorrect astrometry up to ~1"; so use a large matching radius
      ; recall the BOSS fiber diameter is 2"
      spherematch, target.ra, target.dec, dr10.ra, dr10.dec, 1.0/3600.D, maxmatch=0, match0, match1, dist1
      ; remove spectra from the 3 designated SDSS-RM plates
      ind_keep = where( (dr10[match1].plate lt 7338) or (dr10[match1].plate gt 7340) )
      match0 = match0[ind_keep] & match1 = match1[ind_keep]

      dr10 = mrdfits(file,1, rows=match1) ; we can read all columns, columns = cols)
      mwrfits, dr10, outfile, /create
    endif else begin 
      dr10 = mrdfits(outfile, 1)
      spherematch, target.ra, target.dec, dr10.ra, dr10.dec, 1.0/3600.D, maxmatch=0, match0, match1, dist1
      dr10 = dr10[match1]
    endelse

    for i=0L, nobj - 1 do begin
        ind = where(match0 eq i, nmatch)
        if nmatch gt 0 then begin
            mjdarr = dr10[ind].mjd
            minmjd = min(mjdarr, ind_m)
            ind_sort = sort(mjdarr)
            target[i].plate_all = strjoin(string(dr10[ind[ind_sort]].plate, format='(i4.4)'), ' ')
            target[i].fiberid_all = strjoin(string(dr10[ind[ind_sort]].fiberid, format='(i4.4)'), ' ')
            target[i].mjd_all = strjoin(string(dr10[ind[ind_sort]].mjd, format='(i5.5)'), ' ')
            target[i].lambda_eff_all = strjoin(string(dr10[ind[ind_sort]].lambda_eff, format='(i4.4)'), ' ')
            if minmjd eq 0 then print, 'zero mjd detected: '
            target[i].BOSS_TARGET1 = dr10[ind[ind_m]].BOSS_TARGET1
            target[i].BOSS_TARGET2 = dr10[ind[ind_m]].BOSS_TARGET2
            target[i].ANCILLARY_TARGET1 = dr10[ind[ind_m]].ANCILLARY_TARGET1
            target[i].ANCILLARY_TARGET2 = dr10[ind[ind_m]].ANCILLARY_TARGET2
            target[i].plate_sdss_first = dr10[ind[ind_m]].plate
            target[i].fiberid_sdss_first = dr10[ind[ind_m]].fiberid
            target[i].mjd_sdss_first = dr10[ind[ind_m]].mjd
            ; no need for these two below, since lambda_eff_all provides all fiber offset info
            ;target[i].lambda_eff = dr10[ind[ind_m]].lambda_eff
            ;target[i].zoffset = dr10[ind[ind_m]].zoffset
        endif
    endfor    

    ; add DR7 TARGET BITS
    outfile = '/data3/yshen/work/sdssrm_sample_char/dr7_sdssrm.fits'
    if file_test(outfile) ne 1 then begin ; create the file if nonexistant
      file = '/data1/yshen/DR7_spectra/spAll-allquality.fits'
      cols0 = ['ra', 'dec', 'plug_ra', 'plug_dec']
      cols = ['ra', 'dec', 'plug_ra', 'plug_dec', 'plate', 'fiberid', 'mjd', 'PRIMTARGET', 'SECTARGET']
      dr7 = hogg_mrdfits(file,1, columns = cols0)
      spherematch, target.ra, target.dec, dr7.ra, dr7.dec, 1./3600.D, maxmatch=0, match0, match1, dist1
      dr7 = mrdfits(file,1, rows=match1) ; read all columns
      mwrfits, dr7, outfile, /create
    endif else begin
      dr7 = mrdfits(outfile, 1)
      spherematch, target.ra, target.dec, dr7.ra, dr7.dec, 1./3600.D, maxmatch=0, match0, match1, dist1
      dr7 = dr7[match1]
    endelse
    for i=0L, nobj - 1 do begin
        ind = where(match0 eq i, nmatch)
        if nmatch gt 0 then begin
            mjdarr = dr7[ind].mjd
            minmjd = min(mjdarr, ind_m)
            ind_sort = sort(mjdarr)
            target[i].plate_all = strjoin(string(dr7[ind[ind_sort]].plate,format='(i4.4)'), ' ')+ ' ' + target[i].plate_all
            target[i].fiberid_all = strjoin(string(dr7[ind[ind_sort]].fiberid,format='(i4.4)'), ' ')+ ' ' + target[i].fiberid_all
            target[i].mjd_all = strjoin(string(dr7[ind[ind_sort]].mjd,format='(i5.5)'), ' ')+ ' ' + target[i].mjd_all
            target[i].lambda_eff_all = strjoin(string( replicate(5400., nmatch), format='(i4.4)'), ' ') + ' ' + target[i].lambda_eff_all

            target[i].PRIMTARGET = dr7[ind[ind_m]].PRIMTARGET
            target[i].SECTARGET = dr7[ind[ind_m]].SECTARGET
            if (target[i].mjd_sdss_first eq 0) or (minmjd lt target[i].mjd_sdss_first) then begin
                target[i].plate_sdss_first = dr7[ind[ind_m]].plate
                target[i].fiberid_sdss_first = dr7[ind[ind_m]].fiberid
                target[i].mjd_sdss_first = dr7[ind[ind_m]].mjd
            endif
        endif
    endfor

    ; add the 2 missed dr7_qcat objects (which should have been included in dr7 spall-allquality but for some reason didn't
    ; but last time I checked, they are included in dr7 spall-allquality (probably because these 2 objects had bad astrometry in DR7spall. 
    ;file = '/data1/yshen/Project/DR7_BH_catalog/catalog/dr7_bh_more.fits'
    ;bh = mrdfits(file,1)
    ;ind = where(strtrim(target.release) eq 'dr7_qcat')
    ;spherematch, target[ind].ra, target[ind].dec, bh.ra, bh.dec, 0.2/3600.D, maxmatch=0, match0, match1, dist1
    ;target[ind[match0]].PRIMTARGET = bh[match1].TARGET_FLAG_TARGET
    ;target[ind[match0]].plate_sdss_first = bh[match1].plate
    ;target[ind[match0]].fiberid_sdss_first = bh[match1].fiber
    ;target[ind[match0]].mjd_sdss_first = bh[match1].mjd
    ;for i=0L, 1 do begin
    ;   target[ind[match0[i]]].plate_all = string(bh[match1[i]].plate,format='(i4.4)')+ ' ' + target[ind[match0[i]]].plate_all
    ;   target[ind[match0[i]]].fiberid_all = string(bh[match1[i]].fiber,format='(i4.4)')+ ' ' + target[ind[match0[i]]].fiberid_all
    ;   target[ind[match0[i]]].mjd_all = string(bh[match1[i]].mjd,format='(i5.5)')+ ' ' + target[ind[match0[i]]].mjd_all
    ;   target[ind[match0[i]]].lambda_eff_all = string(5400.,format='(i4.4)')+ ' ' + target[ind[match0[i]]].lambda_eff_all
    ;endfor 

    ; add OTHER TARGETS (AGEIS, MMT)
    ind = where( strtrim(target.release) eq 'egs' )
    target[ind].OTHER_TARGET = target[ind].release

    ind = where( strtrim(target.release) eq 'mmt' )
    target[ind].OTHER_TARGET = 'mmt-var'

    ; add the spectral fits    
    qsofit_file = '/data3/yshen/work/lineshifts/lineshift_civ3gauss.fits'
    qso=mrdfits(qsofit_file,1)
    tags = tag_names(qso)
    tags_keep = [tags[0], tags[5:*]]
    qso=mrdfits(qsofit_file,1, column=tags_keep)
    target = struct_addtags(target, qso)


    outfile = '/data3/yshen/work/sdssrm_sample_char/sample_char.fits'
    mwrfits, target, outfile, /create

end

; add matching properties from PS1, wise, FIRST, GALEX GR7 etc
pro add_target_prop

    outfile = '/data3/yshen/work/sdssrm_sample_char/sample_char.fits'
    target = mrdfits(outfile, 1)

    ; match with FIRST
    file = '/data1/yshen/Semester_Projects/Yen-Ting/FIRST_coverage/first_14dec17.fits'
    first = mrdfits(file, 1)
    ; first determine if all sdss-rm quasars are in first footprint using the same 14dec17 first coverage maps
    temp_flag = is_in_FIRST(ra = target.ra, dec = target.dec, /quiet)
    ind = where(temp_flag ne 0)
    target[ind].FIRST_FR_TYPE = 0 ; these are in FIRST footprint
    ; first match within 30" of the FIRST sources
    spherematch, target.ra, target.dec, first.ra, first.dec, 30./3600.D, match1, match2, distance12, maxmatch=0
    ; for each matched SDSS source, determine the 20cm flux and FR type
    uniq_ind = uniq(match1)
    n_match = n_elements(uniq_ind)
    print, 'SDSS object matched: ', n_match
    for i=0L, n_elements(uniq_ind) - 1L do begin
       ind_SDSS = match1[uniq_ind[i]] ; the index of the matched SDSS quasar
       indd = where(match1 eq match1[uniq_ind[i]], n_first)
       ind_FIRST = match2[indd]  ; the indices of the matched FIRST sources, could be multiple
       if n_first ge 2 then begin ; if multiple FIRST sources detected, add all their flux together
         target[ind_SDSS].FIRST_FR_type = 2L
         ;target[ind_SDSS].FIRST_Fpeak_mJy = max(first[ind_FIRST].fpeak)
         target[ind_SDSS].FIRST_Fint_mJy = total(first[ind_FIRST].fint)
       endif else begin ; if only one source detected, see if it is within 5", otherwise, reject
         spherematch, target[ind_SDSS].ra, target[ind_SDSS].dec, first[ind_FIRST].ra, first[ind_FIRST].dec, $
             5./3600.D, match1_5arcsec, match2_5arcsec, dist_5arcsec, maxmatch = 0
         if match1_5arcsec[0] ne -1 then begin ; the SDSS obj does have a match within 5"
           target[ind_SDSS].FIRST_FR_type = 1L
           ;target[ind_SDSS].FIRST_Fpeak_mJy = first[ind_FIRST].fpeak
           target[ind_SDSS].FIRST_Fint_mJy = first[ind_FIRST].fint
         endif
       endelse
       print, 'finished:', i+1, '/', n_match
    endfor

    ; now compute the *observed* restframe-6cm flux density, assuming a power-law slope -0.5, f_nu ~ nu^alpha_nu
    alpha_nu = -0.5 ; this slope comes from Jiang+2007 and Ivezic+2004b
    ind = where(target.FIRST_Fint_mJy gt -1d-5)
    target[ind].Fint_rest6cm_mJy_obs = (20./6.)^alpha_nu*(1. + target[ind].zfinal)^(-alpha_nu)*target[ind].FIRST_Fint_mJy

   
    ;-------------- WISE matches --------------------
    file = '/data3/yshen/work/sdssrm_sample_char/rmtarget_allwise_10arcsec.tbl' 
    readcol, file, format='x,x,x,x,x,x,d,d,x,x,x,d,d,x,x,d,d,x,x,d,d,x,x,d,d', $
      skipline = 23L, ra, dec, w1,w1sig,w2,w2sig, w3,w3sig,w4,w4sig
    spherematch, target.ra, target.dec, ra, dec, 6./3600.D, match1, match2, distance, maxmatch=1  ; get the closest match only
    ; first do a clearning of all ALLWISE data in target
    target.allwise1234 =  -1.D & target.allwise1234_err = -1.D & target.allwise_offset = -1.D
    ; now assign the allwise data
    nnn = n_elements(match1)
    for i=0L, nnn - 1L do begin
       ind = match2[i]
       target[match1[i]].allwise1234 = [w1[ind], w2[ind], w3[ind], w4[ind]]
       target[match1[i]].allwise1234_err = [w1sig[ind], w2sig[ind], w3sig[ind], w4sig[ind]]
       target[match1[i]].allwise_offset = distance[i]*3600.D ; convert to arcsec
    endfor
    
    ; now do unWISE
    file = '/data3/yshen/work/sdsswise/md07_unwise_all.fits'
    if file_test(file) eq 0 then begin ; create this file if non-existent
       topdir = '/data3/yshen/work/sdsswise/'
       ntile=6
       for i=1,ntile do begin
          file1 = topdir + 'sdsswise_tile' + string(i,format='(i0)')+'.fits'
          tmp = mrdfits(file1,1)
          if n_elements(unwise_all) eq 0 then unwise_all = tmp else unwise_all = [unwise_all, tmp]
       endfor
       ; remove duplicate objects
       nnn=n_elements(unwise_all)
       rejflag = lonarr(nnn)
       spherematch, unwise_all.ra, unwise_all.dec, unwise_all.ra, unwise_all.dec, 0.1/3600.D, match1, match2, dist1
       ind = where(match1 gt match2)
       if ind[0] ne -1 then begin 
          rejflag[match1[ind]]=1
          indd = where(rejflag eq 0)
          unwise_all = unwise_all[indd]
       endif
       mwrfits, unwise_all, file, /create
    endif 
    unwise = mrdfits(file, 1)
    spherematch, target.ra, target.dec, unwise.ra, unwise.dec, 6./3600.D, match1, match2, distance, maxmatch=1  ; closest match only
    ; first do a clearning of all ALLWISE data in target
    target.unwise1234 =  -1.D & target.unwise1234_err = -1.D & target.unwise_offset = -1.D
    ; now assign the unwise data
    nnn = n_elements(match1)
    for i=0L, nnn - 1L do begin
       ind = match2[i]
       target[match1[i]].unwise1234 = [unwise[ind].w1_mag, unwise[ind].w2_mag, unwise[ind].w3_mag, unwise[ind].w4_mag]
       target[match1[i]].unwise1234_err = [unwise[ind].w1_mag_err, unwise[ind].w2_mag_err, unwise[ind].w3_mag_err, unwise[ind].w4_mag_err]
       target[match1[i]].unwise_offset = distance[i]*3600.D ; convert to arcsec
    endfor
    ; ------------------------------------------------

    ; ------ GALEX GR6/7; Feb 27, 2013 -----------------
    ; the GALEX cross catalog is retrived from CALEX casjob at https://galex.stsci.edu/casjobs/ with yshen:xxxxxx



    ; overwrite the output
    mwrfits, target, outfile, /create

end


pro add_sample_zsys_ps1

    outfile = '/data3/yshen/work/sdssrm_sample_char/sample_char.fits'
    target = mrdfits(outfile, 1)

    nnn = n_elements(target)
    alltag = tag_names(target)

    ; assign zsys and zsys_err using the recipes in Shen++(2016)
    ; do a clean first
    target.zsys = -1. & target.zsys_err = -1.
    line=['Hbeta_br', 'OIII5007', 'CaII3934', 'OII3728', 'NeV3426', 'MgII', 'CIII', 'HeII1640', 'CIV_br', 'SIIV_OIV']
    wave=[4862.68, 5008.24, 3934.78, 3728.48, 3426.84, 2798.75, 1908.73, 1640.42, 1549.06, (1396.76 + 1402.06)*0.5]
    line=strupcase(line)
    nline = n_elements(line)
    ; range of wavelengths where I still trust the line fit
    maxwave = 10350. & minwave = 3650.
    for i=0l, nnn - 1 do begin
       zsys = dblarr(nline) - 1. & zsys_err = dblarr(nline) - 1.
       rejflag = lonarr(nline) ; rejection falg when the line center falls out [minwave, maxwave]
       for j=0, nline - 1 do begin
           ind1 = where(alltag eq line[j])
           ind2 = where(alltag eq line[j] + '_ERR')
           obs_lam = (target[i].(ind1))[0] & obs_lam_err = (target[i].(ind2))[0]
           zz = target[i].zfinal
           obs_lam = obs_lam*(1. + zz) & obs_lam_err = obs_lam_err*(1. + zz)

           ; the bestfit line peak falls outside the reliable spectral range, so reject
           if obs_lam gt maxwave or obs_lam lt minwave then rejflag[j] = 1L
           logL1700 = target[i].logl1700 & logL1700_err = target[i].logL1700_err
           
           if obs_lam gt 0 then begin
              
              tmp = get_zsys(obs_lam, obs_lam_err, line_use = line[j], logL=logL1700)
              zsys[j] = tmp[0] & zsys_err[j] = tmp[1]
           endif
       endfor
       indd = where(zsys_err gt 0 and rejflag eq 0)
       if indd[0] ne -1 then begin
          zsys = zsys[indd] & zsys_err = zsys_err[indd]
          ;print, line[indd]
          ;print, zsys
          ; find the best zsys with the smallest zsys_err
          tmp2 = min(zsys_err, ind_min)
          target[i].zsys = zsys[ind_min] & target[i].zsys_err = zsys_err[ind_min]
       endif
       splog, 'update zsys for RMID: ', i, target[i].zfinal
       ; pause
    endfor 

    ;; ---------- add PS1 properties ------------
    file = '/data3/yshen/work/PS1/yue_all.fits'
    ps1 = mrdfits(file,1)
    if (where(alltag eq 'PS1_NMAG_OK'))[0] eq -1 then begin ; add the PS1 tags
       newtag = {PS1_Nmag_OK:lonarr(5), PS1_RMS_Mag:dblarr(5)-1.}
       newtag = replicate(newtag, nnn)
       target = struct_addtags(target, newtag)
    endif else begin ; clear all PS1 values
       target.ps1_nmag_ok = 0 & target.ps1_rms_mag = -1.
    endelse
   
    spherematch, target.ra, target.dec, ps1.ra, ps1.dec, 1./3600.D, match0, match1, dist1
    nmatch = n_elements(match0)
    for i=0, nmatch - 1 do begin
        ind1 = match0[i] & ind2 = match1[i]
        ; what magnitude (default is PSF, or Aperture, Kron?)
        lc_mag = ps1[ind2].LC_mag & lc_err = ps1[ind2].LC_err

        nmag_OK = lonarr(5) & rms_mag = dblarr(5)-1.
        for j=0, 4 do begin ; for the 5 PS1 bands
           mag = lc_mag[j, *] & err = lc_err[j, *]
           ind_OK = where(err gt 0 and mag gt 0 and mag lt 30., nok)
           nmag_ok[j] = nok
           if nok gt 5 then begin ; only if there are more than 5 epochs
               mag = mag[ind_ok] & err = err[ind_ok]
               ; using Eqn 2 in Sun++15 (from Sesar+07) to compute intrinsic rms
               term1 = total((mag - mean(mag))^2, /double)/double(nok - 1)
               term2 = mean( err^2 )
               if term1 ge term2 then $
                  rms_mag[j] = sqrt( term1 - term2 ) else rms_mag[j] = 0.
           endif           
        endfor
        target[ind1].ps1_nmag_ok = nmag_ok
        target[ind1].ps1_rms_mag = rms_mag
    endfor

    ;; ------------- add absorber info from Guangtun Zhu ------
    file = '/data3/yshen/work/absorb/QSObased_Expanded_SDSSRM_107.fits' ; this includes both intervene and associated narrow absorbers, 
                                        ; although Guangtun noted the AAL identification may be less secure
    qsoabs = mrdfits(file,1)
    if (where(alltag eq 'NABS'))[0] eq -1 then begin ; add the narrow absorber tags
       newtag = {nabs:0L, zabs:dblarr(10)-1.}
       newtag = replicate(newtag, nnn)
       target = struct_addtags(target, newtag)
    endif else begin ; clear all narrow abs tags
       target.nabs = 0 & target.zabs = dblarr(10) - 1.
    endelse
    ; using index_qso included in qsoabs to locate the RM target
    ind = where(qsoabs.index_qso lt 849) ; remove the object with incorrect index=999
    qsoabs = qsoabs[ind]
    target[qsoabs.index_qso].nabs = qsoabs.nabs
    target[qsoabs.index_qso].zabs = qsoabs.zabs

    ; overwrite the output file
    mwrfits, target, outfile, /create

end

; clean the characterization sample by including only the useful columns
pro clean_char_sample

    file = '/data3/yshen/work/sdssrm_sample_char/sample_char.fits'

    columns = ['RM_ID', 'RA', 'DEC', 'PSFMAG', 'FIBER2MAG', 'OBJC_TYPE', 'ZFINAL',  $
'ZSYS', 'ZSYS_ERR', 'BOSS_TARGET1',  $
'BOSS_TARGET2', 'ANCILLARY_TARGET1', 'ANCILLARY_TARGET2', 'PRIMTARGET',  $
'SECTARGET', 'OTHER_TARGET', 'PLATE_SDSS_FIRST', 'FIBERID_SDSS_FIRST',  $
'MJD_SDSS_FIRST', 'PLATE_ALL', 'FIBERID_ALL', 'MJD_ALL', 'LAMBDA_EFF_ALL', $
'PS1_NMAG_OK', 'PS1_RMS_MAG', $
'ALLWISE1234', 'ALLWISE1234_ERR', 'ALLWISE_OFFSET', 'UNWISE1234',  $
'UNWISE1234_ERR', 'BAL_FLAG', 'NABS', 'ZABS', 'FIRST_FR_TYPE', $
'FIRST_FINT_MJY', 'FINT_REST6CM_MJY_OBS', 'LOGFNU2500A_ERGS_OBS', 'R_6CM_2500A', $
'CONTI_FIT', 'CONTI_FIT_ERR', 'FEII_UV', 'FEII_UV_ERR', 'FEII_OPT', 'FEII_OPT_ERR', $
'CONTI_REDCHI2', 'LOGL1350', 'LOGL1350_ERR', 'LOGL1700', 'LOGL1700_ERR', 'LOGL3000', 'LOGL3000_ERR', $
'LOGL5100', 'LOGL5100_ERR', 'REW_FE_4434_4684', 'REW_FE_4434_4684_ERR', $
'SII6718', 'HALPHA', 'HALPHA_BR', 'HBETA', 'HBETA_BR', 'HEII4687', 'HEII4687_BR', $
'OIII5007', 'OIII5007C', 'CAII3934', 'OII3728', 'NEV3426', $
'MGII', 'MGII_BR', 'CIII', 'CIII_BR', 'SIIII1892', $
'ALIII1857', 'NIII1750', 'CIV', 'HEII1640', 'HEII1640_BR', 'SIIV_OIV', $
'OI1304', 'LYA', 'NV1240', 'SII6718_REDCHI2', 'HALPHA_REDCHI2',  $
'HALPHA_BR_REDCHI2', 'HBETA_REDCHI2', 'HBETA_BR_REDCHI2', 'HEII4687_REDCHI2', $
'HEII4687_BR_REDCHI2', 'OIII5007_REDCHI2', 'OIII5007C_REDCHI2', 'CAII3934_REDCHI2', 'OII3728_REDCHI2', 'NEV3426_REDCHI2', $
'MGII_REDCHI2', 'MGII_BR_REDCHI2', 'CIII_REDCHI2', 'CIII_BR_REDCHI2', 'SIIII1892_REDCHI2', $
'ALIII1857_REDCHI2', 'NIII1750_REDCHI2', 'CIV_REDCHI2', $
'HEII1640_REDCHI2', 'HEII1640_BR_REDCHI2', 'SIIV_OIV_REDCHI2', 'OI1304_REDCHI2', $
'LYA_REDCHI2', 'NV1240_REDCHI2', 'SII6718_ERR', 'HALPHA_ERR', 'HALPHA_BR_ERR', $
'HBETA_ERR', 'HBETA_BR_ERR', 'HEII4687_ERR', 'HEII4687_BR_ERR', 'OIII5007_ERR', $
'OIII5007C_ERR', 'CAII3934_ERR', 'OII3728_ERR', 'NEV3426_ERR', $
'MGII_ERR', 'MGII_BR_ERR', 'CIII_ERR', 'CIII_BR_ERR', $
'SIIII1892_ERR', 'ALIII1857_ERR', 'NIII1750_ERR', 'CIV_ERR', $
'HEII1640_ERR', 'HEII1640_BR_ERR', 'SIIV_OIV_ERR', 'OI1304_ERR', 'LYA_ERR', $
'NV1240_ERR', 'LOGBH_CIV_VP06', 'LOGBH_CIV_VP06_ERR', 'LOGBH_MGII_S11', $
'LOGBH_MGII_S11_ERR', 'LOGBH_HB_VP06', 'LOGBH_HB_VP06_ERR', 'COMMENT']

    result = mrdfits(file, 1, column=columns)

    ; need to replace some tagnames
    oldtag = ['RM_ID', 'ZFINAL', 'CIII', 'CIII_REDCHI2', 'CIII_ERR']
    newstr = {RMID:0L, ZPIP:0.D, CIII_all:dblarr(5), CIII_all_REDCHI2:0.d, CIII_all_ERR:dblarr(5)-1.}
    newstr = replicate(newstr, n_elements(result))
    result = struct_addtags(result, newstr)
    result.RMID= result.rm_ID & result.zpip = result.zfinal & result.CIII_all = result.CIII 
    result.CIII_all_redchi2 = result.CIII_redchi2 & result.CIII_all_err = result.CIII_err

    outfile = '/data3/yshen/work/sdssrm_sample_char/sample_char_refined.fits'
    mwrfits, result, outfile, /create

    ; now read again with the correct tag order and names:
columns = ['RMID', 'RA', 'DEC', 'PSFMAG', 'FIBER2MAG', 'OBJC_TYPE', 'ZPIP',  $
'ZSYS', 'ZSYS_ERR', 'BOSS_TARGET1',  $
'BOSS_TARGET2', 'ANCILLARY_TARGET1', 'ANCILLARY_TARGET2', 'PRIMTARGET',  $
'SECTARGET', 'OTHER_TARGET', 'PLATE_SDSS_FIRST', 'FIBERID_SDSS_FIRST',  $
'MJD_SDSS_FIRST', 'PLATE_ALL', 'FIBERID_ALL', 'MJD_ALL', 'LAMBDA_EFF_ALL', $
'PS1_NMAG_OK', 'PS1_RMS_MAG', $
'ALLWISE1234', 'ALLWISE1234_ERR', 'ALLWISE_OFFSET', 'UNWISE1234',  $
'UNWISE1234_ERR', 'BAL_FLAG', 'NABS', 'ZABS', 'FIRST_FR_TYPE', $
'FIRST_FINT_MJY', 'FINT_REST6CM_MJY_OBS', 'LOGFNU2500A_ERGS_OBS', 'R_6CM_2500A', $
'CONTI_FIT', 'CONTI_FIT_ERR', 'FEII_UV', 'FEII_UV_ERR', 'FEII_OPT', 'FEII_OPT_ERR', $
'CONTI_REDCHI2', 'LOGL1350', 'LOGL1350_ERR', 'LOGL1700', 'LOGL1700_ERR', 'LOGL3000', 'LOGL3000_ERR', $
'LOGL5100', 'LOGL5100_ERR', 'REW_FE_4434_4684', 'REW_FE_4434_4684_ERR', $
'SII6718', 'HALPHA', 'HALPHA_BR', 'HBETA', 'HBETA_BR', 'HEII4687', 'HEII4687_BR', $
'OIII5007', 'OIII5007C', 'CAII3934', 'OII3728', 'NEV3426', $
'MGII', 'MGII_BR', 'CIII_ALL', 'CIII_BR', 'SIIII1892', $
'ALIII1857', 'NIII1750', 'CIV', 'HEII1640', 'HEII1640_BR', 'SIIV_OIV', $
'OI1304', 'LYA', 'NV1240', 'SII6718_REDCHI2', 'HALPHA_REDCHI2',  $
'HALPHA_BR_REDCHI2', 'HBETA_REDCHI2', 'HBETA_BR_REDCHI2', 'HEII4687_REDCHI2', $
'HEII4687_BR_REDCHI2', 'OIII5007_REDCHI2', 'OIII5007C_REDCHI2', 'CAII3934_REDCHI2', 'OII3728_REDCHI2', 'NEV3426_REDCHI2', $
'MGII_REDCHI2', 'MGII_BR_REDCHI2', 'CIII_ALL_REDCHI2', 'CIII_BR_REDCHI2', 'SIIII1892_REDCHI2', $
'ALIII1857_REDCHI2', 'NIII1750_REDCHI2', 'CIV_REDCHI2', $
'HEII1640_REDCHI2', 'HEII1640_BR_REDCHI2', 'SIIV_OIV_REDCHI2', 'OI1304_REDCHI2', $
'LYA_REDCHI2', 'NV1240_REDCHI2', 'SII6718_ERR', 'HALPHA_ERR', 'HALPHA_BR_ERR', $
'HBETA_ERR', 'HBETA_BR_ERR', 'HEII4687_ERR', 'HEII4687_BR_ERR', 'OIII5007_ERR', $
'OIII5007C_ERR', 'CAII3934_ERR', 'OII3728_ERR', 'NEV3426_ERR', $
'MGII_ERR', 'MGII_BR_ERR', 'CIII_ALL_ERR', 'CIII_BR_ERR', $
'SIIII1892_ERR', 'ALIII1857_ERR', 'NIII1750_ERR', 'CIV_ERR', $
'HEII1640_ERR', 'HEII1640_BR_ERR', 'SIIV_OIV_ERR', 'OI1304_ERR', 'LYA_ERR', $
'NV1240_ERR', 'LOGBH_CIV_VP06', 'LOGBH_CIV_VP06_ERR', 'LOGBH_MGII_S11', $
'LOGBH_MGII_S11_ERR', 'LOGBH_HB_VP06', 'LOGBH_HB_VP06_ERR', 'COMMENT']
    result = mrdfits(outfile, 1, column=columns)
    ; now overwrite
    mwrfits, result, outfile, /create


end

; add radio R parameter and LOGFNU2500A_ERGS_OBS
pro add_radio_R

file = '/data3/yshen/work/sdssrm_sample_char/sample_char_refined.fits'
tt = mrdfits(file,1)

nnn = n_elements(tt)

for i=0, nnn - 1 do begin

  pp = tt[i].conti_fit
  wave_obs = 2500. * (1. + tt[i].zpip)
  if wave_obs gt 3600. and wave_obs lt 10300. then begin
    Fnu2500A_ergs_OBS = f_conti_only(2500., pp[6:*])*1d-9*(2.5d-5*(1. + tt[i].zpip))^2/2.9979246d10
    Fnu2500A_ergs_OBS = Fnu2500A_ergs_OBS > 1d-30 ;; corresponds to ~ 5d-19 erg/s/cm2/A for f_lambda
    if Fnu2500A_ergs_OBS gt 1d-30 then begin
      tt[i].LOGFNU2500A_ERGS_OBS = alog10(Fnu2500A_ergs_OBS)
      if Fnu2500A_ergs_OBS gt 1d-35 and tt[i].FINT_REST6CM_MJY_OBS gt 0 then $
         tt[i].R_6CM_2500A = tt[i].FINT_REST6CM_MJY_OBS*1d-26/Fnu2500A_ergs_OBS
    endif
  endif

endfor

mwrfits, tt, file, /create

end

;; --- compare the improvement in zsys using coadded spectra
pro plot_zsys_comp, vel_plot=vel_plot

    cs = 2.9979246d5

    file = '/data3/yshen/work/sdssrm_sample_char/med_coadd_zpip.fits'
    spec_pip = mrdfits(file,1)
    file = '/data3/yshen/work/sdssrm_sample_char/med_coadd_zsys.fits'
    spec_sys = mrdfits(file,1)
    wave = spec_pip.wave
    spec1 = spec_pip.coadd_flux & err1 = spec_pip.coadd_err
    spec2 = spec_sys.coadd_flux & err2 = spec_sys.coadd_err

    figfile = '/data3/yshen/work/sdssrm_sample_char/figs/zsys_comp1.ps'
    begplot, name=figfile, /color, /landscape

    line = textoidl(['HeII1640', 'CIV', '[NeV]3426', '[OII]3727', 'CaII K', 'H\beta', '[OIII]5007'])
    maxv = [4000., 4000., 2000., 1000., 1000., 2000., 1000.]
    wrange = [ [1610, 1660], $   
               [1450, 1600], $               ; [1700, 1950], CIII] $
               [3310, 3460], $
               [3710, 3770], $
               [3910, 3960], $
               [4800, 4950], $
               [4900, 5050] ]
    linewave = [1640.42, 1549.06, 3426.84, 3728.48, 3934.78, 4862.68, 5008.24]  ; 1908.73

    ; use_line=[0,2,3,4,5,6]
    use_line = [0,1]
    line=line[use_line]
    maxv = maxv[use_line]
    wrange = wrange[*, use_line]
    linewave = linewave[use_line]

    nline = n_elements(linewave)

    nplot = n_elements(wrange[0,*])
    ang = string(197B)
    omargin=[0.1, 0.01, 0.99, 0.95]
    pmargin=[0.1,0.21]
    charsize=1.2
    plot_layout, nplot, xypos=xypos, omargin=omargin, pmargin=pmargin, $
     nrow=2

    for i=0, nplot - 1 do begin
        if i eq 0 then noerase = 0 else noerase = 1L
        ind = where(wave ge wrange[0,i] and wave le wrange[1,i] )
        yrange = [min(spec2[ind]), max(spec2[ind])*1.2]
        if keyword_set(vel_plot) then begin
           xrange = [-maxv[i], maxv[i]]  ;[-1000.,1000.]
           xarr = (wave[ind]/linewave[i] - 1) * cs
        endif else begin
           xrange = wrange[*,i]
           xarr = wave[ind]
        endelse
        
        plot, xarr, spec1[ind], xrange=xrange,/xsty, yrange=yrange, pos=xypos[*,i], noerase=noerase, $
           psym=10, title=line[i], charsize=charsize
        oplot, xarr, spec2[ind], color=cgcolor('red'), psym=10

        if ~keyword_set(vel_plot) then begin
           for ii=0, nline - 1 do begin
               if linewave[ii] gt xrange[0] and linewave[ii] lt xrange[1] then $
                  oplot, [linewave[ii], linewave[ii]], [0,60],line=1
           endfor
        endif else oplot, [0,0],[12,60],line=1
    endfor
    if ~keyword_set(vel_plot) then xtitle = textoidl('Restframe \lambda_{vac}') + ' [' +ang + ']' $
       else xtitle = textoidl('Velocity [km s^{-1}]')
    xyouts, 0.5, 0.02, align=0.5, xtitle,/norm
    xyouts, 0.04, 0.5, align=0.5, orient=90, textoidl('Flux Density f_\lambda [arbitrary units]'),/norm

    endplot
    cgfixps, figfile

end


; plot the PS1 intrinsic RMS magnitudes
pro plot_ps1_rms

file = '/data3/yshen/work/sdssrm_sample_char/sample_char_refined.fits'
tt = mrdfits(file,1)

figfile='/data3/yshen/work/sdssrm_sample_char/figs/ps1_rmsmag.eps'
begplot, name=figfile, /color,/encap, xsize=8, ysize=6

rms = tt.ps1_rms_mag
zz = tt.zsys

grms = rms[0, *] & irms = rms[2,*]
indg = where(grms gt 0) & indi = where(irms gt 0)

print, median(grms), median(irms)

plot, zz[indg], grms[indg], psym=5, xtitle = 'Redshift', ytitle = 'PS1 RMS Magnitude', /ylog, /nodata, $
  yrange=[0.05, 2], pos=[0.11,0.11,0.98,0.98],/ysty
oplot, zz[indg], grms[indg], psym=5, color=cgcolor('green'), thick=3
oplot, zz[indi], irms[indi], psym=1, color=cgcolor('red')

items = textoidl([' g_{PS1}', ' i_{PS1}'])
colors=cgcolor(['green','red'])
legend, items, /norm, box=0, pos=[0.8,0.95], psym=[5,1], color=colors, textcolor=colors, spac=2

endplot

end

; plot Mi(z2)-z distribution
pro plot_Lz

file = '/data3/yshen/work/sdssrm_sample_char/sample_char_refined.fits'
tt = mrdfits(file,1)
imag = (tt.psfmag)[3,*]
Miz2 = get_abs_mag(imag, tt.zpip)

figfile = '/data3/yshen/work/sdssrm_sample_char/figs/L_z.eps'
begplot, name=figfile, /color,/encap, xsize=6, ysize=5.8

plot, [0], [0], xrange = [0, 5], yrange = [-19.9, -30.1], /xsty, /ysty, $
    xtitle = textoidl('Redshift'), ytitle=textoidl('M_i (z=2)'), /nodata, $
    pos = [0.15, 0.12, 0.98, 0.98]
oplot, tt.zpip, Miz2, psym=symcat(9,thick=3), symsize=0.5

zgrid = 0.2 + 0.1*indgen(50)
imag=get_abs_mag(21.7,zgrid)
oplot, zgrid, imag, thick=6., color=cgcolor('red'), linestyle = 2
xyouts, 4.2, -25.25, 'i=21.7',color=cgcolor('red')

endplot
end

; make a collection of example SDSS-RM quasars
pro plot_examples

file = '/data3/yshen/work/sdssrm_sample_char/sample_char_refined.fits'
tt = mrdfits(file,1)

rmid = [3, 8, 34, $
        85, $  ; disk emitter
        90, $
        116, $ ; NALs
        272, 320, $ ; host stellar absorption
        39, 226, 613, $ ; BALs]
        316, $  ; strong optical FeII
        325, $  ; associated LLS from the BLR?
        376, $  ; lots of unusual narrow lines
        23, 44, 484]  ; N-loud

figfile = '/data3/yshen/work/sdssrm_sample_char/figs/examp.eps'
begplot, name=figfile, /color, xsize=8, ysize=5
ang = string(197B)
nnn = n_elements(rmid)
rest_wave = 1000. + findgen(8)*1000.
xtickname = string(rest_wave, format='(i0)')
xrange=[3600,1d4]
for i=0, nnn - 1 do begin
  rm_readspec,0,rmid[i]+1, mjd=56837,wave=wave, flux=flux, invvar=ivar
  ind = where(wave ge 5573 and wave le 5585, complement=indd)
  flux[ind] = interpol(flux[indd], wave[indd], wave[ind])

  if rmid[i] eq 484 then yrange=[0,5] else ttt=temporary(yrange)

  ind=where(ivar gt 0)
  plot, wave[ind], smooth(flux[ind],3), xsty=5, xrange=xrange, $
      ytitle= textoidl('Flux Density [10^{-17 }erg s^{-1}cm^{-2}')+ang+textoidl('^{-1}]'), $
       thick=2, pos=[0.12, 0.13, 0.95, 0.88],yrange=yrange
  xyouts, 0.68, 0.8, /norm, 'RMID'+string(rmid[i],format='(i3.3)') + '   z='+string(tt[rmid[i]].zsys,format='(f0.3)')
  oplot, xrange, [0,0], line=1 

  axis, xaxis=0, xrange=xrange, /xsty, xtitle='Observed Wavelength [' + ang + ']'
  axis, xaxis=1, xrange=xrange, xtickv=rest_wave*(1. + tt[rmid[i]].zsys), xminor=10, xticks=7,/xsty, $
     xtickname=xtickname, xtitle='Rest Wavelength [' + ang + ']'

endfor
endplot

;cgfixps, figfile

end

pro mk_add_linefits_out
; add addtional lines to the qso-prop.fits file
linename = 'NII6585'
outdir = '/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/civ_3gauss/'
outfile = 'qso_prop-0000-56837_NII-only.fits'

rm_compile_qsofit,0,56837, outdir_in=outdir, outfile=outfile,linename=linename, add_bhmass=0

end

pro add_Lbol_fid_BH_Edd

; add bolometric lum, fiducial BH mass and Eddington ratios
file = '/data3/yshen/work/sdssrm_sample_char/sample_char_refined_w_var.fits'
tt = mrdfits(file,1)
nn=n_elements(tt)
str = create_struct('NII6585', dblarr(5), 'NII6585_ERR', dblarr(5)-1.D, 'REW_FE_2250_2650', 0.D, 'REW_FE_2250_2650_err',-1.D, $
       'logLbol', 0.D, 'logLbol_err', -1.D, 'logbh', 0.D, 'logbh_err', -1.D, $
       'logedd_ratio', -99.D, 'logedd_ratio_err', -1.D)
str = replicate(str, nn)
tt = struct_addtags(tt, str)

linefile = '/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/civ_3gauss/qso_prop-0000-56837_NII-only.fits'
tst = mrdfits(linefile,1)
tt.NII6585 = tst.NII6585 & tt.NII6585_err = tst.NII6585_err
tt.REW_FE_2250_2650 = tst.REW_FE_2250_2650 & tt.REW_FE_2250_2650_err = tst.REW_FE_2250_2650_err

; add logLBOL
ind = where(tt.logL3000 gt 0 and tt.logL3000_err gt 0)
tt[ind].logLbol = tt[ind].logL3000 + alog10(5.15D)
tt[ind].logLbol_err = tt[ind].logL3000_err
ind = where( (tt.logL3000 lt 0 or tt.logL3000_err le 0) and tt.logL5100 gt 0 and tt.logL5100_err gt 0)
if ind[0] ne -1 then begin
  tt[ind].logLbol = tt[ind].logL5100 + alog10(9.26D)
  tt[ind].logLbol_err = tt[ind].logL5100_err
endif
ind = where( (tt.logL3000 lt 0 or tt.logL3000_err le 0) and tt.logL1350 gt 0 and tt.logL1350_err gt 0)
if ind[0] ne -1 then begin
  tt[ind].logLbol = tt[ind].logL1350 + alog10(3.81D)
  tt[ind].logLbol_err = tt[ind].logL1350_err
endif

; add default logbh
ind = where(tt.logbh_hb_vp06 gt 0 and tt.logbh_hb_vp06_err gt 0)
tt[ind].logbh = tt[ind].logbh_hb_vp06 & tt[ind].logbh_err = tt[ind].logbh_hb_vp06_err
ind = where( (tt.logbh_hb_vp06 lt 0 or tt.logbh_hb_vp06_err le 0) and tt.logbh_mgii_s11 gt 0 and tt.logbh_mgii_s11_err gt 0)
tt[ind].logbh = tt[ind].logbh_mgii_s11 & tt[ind].logbh_err = tt[ind].logbh_mgii_s11_err
ind = where( (tt.logbh_mgii_s11 lt 0 or tt.logbh_mgii_s11_err le 0) and tt.logbh_civ_vp06 gt 0 and tt.logbh_civ_vp06_err gt 0)
tt[ind].logbh = tt[ind].logbh_civ_vp06 & tt[ind].logbh_err = tt[ind].logbh_civ_vp06_err

; add Eddington ratios
ind = where(tt.logbh gt 0 and tt.loglbol gt 0)
tt[ind].logedd_ratio = tt[ind].logLbol - alog10(1.3d38) - tt[ind].logBH
tt[ind].logedd_ratio_err = sqrt( tt[ind].logLbol_err^2 + tt[ind].logBH_err^2 )

; add BAL flag from Pat Hall
file = '/data3/yshen/work/sdssrm_sample_char/ancil_data/spec_2014_BALrobust.csv'
readcol, file, format='l,l', rmid, bal_flag, delimiter=','
nobj=n_elements(bal_flag)
for i=0, nobj - 1 do begin
  ind = where(tt.rmid eq rmid[i])
  tt[ind].BAL_Flag = bal_flag[i]
endfor

outfile = '/data3/yshen/work/sdssrm_sample_char/sample_char_final.fits'
mwrfits, tt, outfile, /create

end

pro plot_ml_sdssrm
; plot the SDSS-RM M-L plane

file = '/data1/yshen/Project/DR7_BH_catalog/catalog/dr7_bh_more.fits'
bh = mrdfits(file, 1)
ind = where(bh.LOGBH gt 0 and bh.logLbol gt 0)

logbh = bh[ind].logbh & loglbol = bh[ind].loglbol

figfile = '/data3/yshen/work/sdssrm_sample_char/figs/ml_plane.eps'
begplot, name=figfile, xsize=5.5,ysize=5, /color,/encap
plot, logbh, loglbol, xrange=[6.5,10.7], yrange=[44.,48.2], xtitle=textoidl('M_{BH, SE} [M')+sunsymbol()+']', $
   ytitle = textoidl('log L_{bol} [erg s^{-1}]'), /nodata,pos=[0.17,0.14,0.95,0.98],/xsty, /ysty

oplot, [6,11], [6,11] + alog10(1.3d38), line=2
oplot, [6,11], [6,11] + alog10(1.3d38) - 1., line=2
oplot, [6,11], [6,11] + alog10(1.3d38) - 2., line=2
xyouts, 9.,47.7, textoidl('\lambda_{Edd}=1'),/data,charsize=1.2
xyouts, 9.8,47.6, textoidl('\lambda_{Edd}=0.1'),/data, charsize=1.2
xyouts, 9.82,45.8, textoidl('\lambda_{Edd}=0.01'),/data, charsize=1.2
djs_contourpts, logbh, loglbol, bin1=0.2, bin2=0.2, nlevels=8, thick=6, color=cgcolor('blue'), /over,/nopoints

file ='/data3/yshen/work/sdssrm_sample_char/sample_char_final.fits'
tt=mrdfits(file,1)
ind=where(tt.logbh gt 0 and tt.loglbol gt 0)
oplot, tt[ind].LOGBH, tt[ind].logLbol,psym=symcat(9),color=cgcolor('red'),symsize=0.3
endplot


end
