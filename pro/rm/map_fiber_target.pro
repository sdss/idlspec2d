; create the list of objects (targets, std and sky) that received a fiber
; for all the epochs

; output the target_fibermap.fits file

pro map_fiber_target, data_path = data_path

   ; using the plugmap structure to get the fiberid
   ; modify the plate and mjd list when new epochs become available
   plate = [7338, 7338, 7338, 7339, 7339, 7339, 7339, 7339, 7338, 7339, $
            7339, 7340, 7339, 7338, 7339, 7339, 7339, 7339, 7339, 7339, $
            7339, 7339, 7339, 7339, 7339, 7339, 7339, 7339, $
            7340, 7340, 7340, 7340 , $
            7338, 7338, 7338, 7339, 7339, 7340, 7338, 7338, 7338, 7338, 7340, 7340, $  ; 2015 data
            7339, 7339, 7339, 7339, 7339, 7338, 7339, 7339, 7339, 7339, 7339, 7340, 7340, $  ; 2016 data
            7338, 7338, 7338, 7338, 7338, 7338, 7338, 7338, 7339, 7339, 7339, 7339, $ ; 2017 data
            7338, 7338, 7338, 7338, 7340, 7340, 7340, 7340, 7340]   ; 2018 data
   mjd = [56660L, 56664L, 56669L, 56683L, 56686L, 56697L, 56713L, 56715L, 56717L, 56720, $
          56722, 56726, 56739, 56745, 56747, 56749, 56751, 56755, 56768, 56772, $
          56780, 56782, 56783, 56795, 56799, 56804, 56808, 56813, $
          56825, 56829, 56833, 56837, $
          57038, 57050, 57067, 57082, 57097, 57106, 57127, 57135, 57159, 57166, 57185, 57196, $  ; 2015 data
          57428, 57435, 57451, 57463, 57481, 57490, 57492, 57510, 57518, 57544, 57550, 57570, 57576, $  ; 2016 data
          57781, 57789, 57805, 57817, 57832, 57843, 57859, 57874, 57892, 57901, 57918, 57934, $ ; 2017 data
          58127, 58146, 58174, 58201, 58216, 58230, 58258, 58275, 58289] ; 2018 data

   nep_tot=n_elements(plate)

   etcdir = getenv('IDLRM_DIR') + '/etc/'

   if not keyword_set(data_path) then $ ; '/data3/quasar/yshen/spectro/bossredux/v5_6_0/'
      data_path = getenv('BOSS_SPECTRO_REDUX')+'/'+getenv('RUN2D')+'/'
   data_path_eboss = getenv('BOSS_SPECTRO_REDUX')+ '/v5_10_10/' ; eBOSS data since 2015

   file = etcdir + 'RM_targets.fits'

   target = mrdfits(file,1)

   ; keep only those received a fiber
   ind = where(target.tiled eq 1, nnn)
   target = target[ind]
   target = struct_addtags(target, replicate({plate:lonarr(nep_tot),fiberid:lonarr(nep_tot) $
      , mjd:lonarr(nep_tot),med_sn:dblarr(nep_tot),zfinal:0.D,dr_plate:0L,dr_fiberid:0L,dr_mjd:0L},nnn))
   struct = {ra:0.D,dec:0.D,sourcetype:'', priority:999L, psfmag:fltarr(5) $
      , fiber2mag:fltarr(5), z:0.D, objc_type:0L, release:'', indx:-1L $
      , tiled:1L, plate:lonarr(nep_tot), fiberid:lonarr(nep_tot), mjd:lonarr(nep_tot) $
      , med_sn:dblarr(nep_tot),zfinal:0.D,dr_plate:0L, dr_fiberid:0L, dr_mjd:0L}

   ; add previously known SDSS-DR7 and BOSS spectrum plate-fiber-mjd
   file = etcdir + 'targets_final.fits'
   alltarget = mrdfits(file,1)
   target.dr_plate = alltarget[target.indx].plate
   target.dr_fiberid = alltarget[target.indx].fiberid
   target.dr_mjd = alltarget[target.indx].mjd

   ; now add additonal tiled stuff: 1 Lowz target, 70 std and 80 sky
   file = etcdir + 'plateHoles-007338.par'
   yanny_read, file, pdata
   data = *pdata

   ; assign fiberid for the RM targets
   ; do not use the fiberid from plateHoles*.par file [these are designing fiberID]
   ; use the plugmapM information from spPlate HUD5 instead.

   ; fix the 7 egs RM targets of psfmag
   ; NOTE that the DR10 psfmag for these 7 egs quasars are different from 
   ; the i-band magnitudes that Jon Trump sent to me, which are based on DEEP2
   file = etcdir + 'rm_targets_dr10_photo_primary.csv'
   readcol, file,format='x,x,d,d,x,d,d,d,d,d', ra, dec, mu, mg, mr, mi, mz
   spherematch, target.ra, target.dec, ra, dec, 0.1/3600.D, match0, match1, distance
   ind = where(target[match0].psfmag[0] eq 0, nnn)
   for i=0L, nnn - 1L do begin
     print, target[match0[ind[i]]].release, target[match0[ind[i]]].psfmag
     target[match0[ind[i]]].psfmag = [mu[match1[ind[i]]], mg[match1[ind[i]]] $
       , mr[match1[ind[i]]], mi[match1[ind[i]]], mz[match1[ind[i]]] ]
     print, target[match0[ind[i]]].release, target[match0[ind[i]]].psfmag
   endfor

   ; the LRG
   ind = where(data.sourcetype eq 'LRG ', nnn)
   new = replicate(struct, nnn)
   new.ra = data[ind].target_ra & new.dec = data[ind].target_dec
   new.sourcetype = data[ind].sourcetype & new.release = 'LOWZ'
   target = [target, new]

   ; std
   ind=where(data.sourcetype eq 'STD', nnn)
   new = replicate(struct, nnn)
   new.ra = data[ind].target_ra & new.dec = data[ind].target_dec
   new.sourcetype = data[ind].sourcetype & new.release = 'STD'
   ; now populate photometry for standard stars
   file = etcdir + 'cas_dr10_star_match.fits'
   star = mrdfits(file, 1)
   spherematch, new.ra, new.dec, star.ra, star.dec, 0.2/3600.D, match0,match1, distance
   for i=0L, n_elements(match0) - 1L do begin
     new[match0[i]].psfmag = [star[match1[i]].psfmag_u $
        , star[match1[i]].psfmag_g, star[match1[i]].psfmag_r $
        , star[match1[i]].psfmag_i, star[match1[i]].psfmag_z ]
   endfor
   target = [target, new]

   ; sky
   ind = where(data.targettype eq 'SKY', nnn)
   new = replicate(struct, nnn)
   new.ra = data[ind].target_ra & new.dec = data[ind].target_dec
   new.sourcetype = 'SKY' & new.release = 'SKY'
   target = [target, new]

   ; using the plugmap structure to get the fiberid
   spfile = string(plate,format='(i4.4)') + '-' + string(mjd,format='(i5.5)')
   n_epoch = n_elements(spfile)
   for i=0L, n_epoch - 1L do begin
     if mjd[i] lt 56838L then file = $
        data_path + string(plate[i],format='(i4.4)') + '/spPlate-' + spfile[i] + '.fits' $
     else file = data_path_eboss + string(plate[i],format='(i4.4)') + '/spPlate-' + spfile[i] + '.fits'
     plugmap = mrdfits(file,5)
     spherematch, target.ra, target.dec, plugmap.ra, plugmap.dec, 0.1/3600.D, match0,match1,distance
     target[match0].plate[i,*] = plate[i]
     target[match0].mjd[i,*] = mjd[i]
     target[match0].fiberid[i,*] = transpose(plugmap[match1].fiberid)
   endfor

   ; assign the median SN
   med_SN = dblarr(nep_tot, 1000)
   for i=0L, n_epoch-1 do begin
      if i le 31 then begin
        run2d='v5_7_1' & run1d='v5_7_1'
      endif else begin
        run2d='v5_10_10' & run1d='v5_10_10'
      endelse
      rm_readspec,plate[i], (target.fiberid)[i,*], mjd=mjd[i],flux=flux,invvar=ivar,run1d=run1d, run2d=run2d
      med_SN[i,*] = median(flux*sqrt(ivar), dim=1)
   endfor
   target.med_sn = med_sn

   ; NOW assign the final redshift (zfinal) from zans+VI of coadded spectra
   zlistfile=etcdir+'zfinal.list'
   readcol, zlistfile,format='l,d',fiber, zfinal
   target[0:848].zfinal = zfinal

   output = etcdir + 'target_fibermap.fits'
   mwrfits, target, output, /create

   ; add a second extension for other per plate information
   rm_output_spec_obs, result=result
   mwrfits, result, output

end


