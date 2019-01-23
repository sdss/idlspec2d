;+
; NAME:
;   rm_expinfo
;
; PURPOSE:
;  Return individual exposures information for an RM-BOSS epoch
;
; CALLING SEQUENCE:
;  rm_expinfo, Plate, MJD, expinfo=expinfo
;
; INPUTS:
;   plate	- plate for the spPlate coadd file [e.g., 7338]
;   mjd		- mjd for the coadd file [e.g., 56660]
;
; OUTPUTS:
;   expinfo	- structure containing the info of individual exposures
;
; OPTIONAL OUTPUTS:
;   mean_mjd	- mean MJD of all coadded frames
;   exptime	- total exposure time of all coadded frames
;
; REVISION HISTORY:
;   05-Jan-2014  Written by Yue Shen, Carnegie
;-
;------------------------------------------------------------------------
pro rm_expinfo, plate, mjd, expinfo=expinfo, mean_mjd=mean_mjd, exptime=exptime $
       , diet=diet

   topdir = getenv('BOSS_SPECTRO_REDUX')
   if mjd le 56838L then twoddir = getenv('RUN2D') else twoddir = 'v5_10_10' ; eBOSS data uses different RUN1d/2d
   platestr = string(plate,format='(i4.4)')
   mjdstr = string(mjd,format='(i5.5)')

   ; find all the exposures in the spPlancomb file
   filename = 'spPlancomb-' + platestr + '-' + mjdstr + '.par'
   filename = lookforgzip(filepath(filename, root_dir=topdir, $
      subdirectory=[twoddir,platestr]), count=ct)
   if (ct eq 1) then filename = filename[0] else $
         message, 'spPlancomb file not found: ',plate,mjd
   yanny_read, filename, pdata
   result = *pdata
   nexp = n_elements(result)
   struct_tmp = replicate({mjd_beg:dblarr(4), mjd_end:dblarr(4) $
	, airmass:fltarr(500),coadded:lonarr(4), platesn2:0.D, seeing20:0. $
        , seeing50:0., seeing80:0., plate_airmass:0.}, nexp)
   result = struct_addtags(result, struct_tmp)

   ;if ~keyword_set(diet) then begin
   for i=0L, nexp - 1L do begin
      mjd_beg = dblarr(4) & mjd_end = dblarr(4)
      for j=0L, 3L do begin

      spFramefile = topdir+'/'+twoddir+'/'+platestr+'/'+result[i].name[j] + '.gz'
         ;print, spFramefile
         if file_test(spFramefile) ne 0 then begin
  
            if keyword_Set(diet) then hdr1 = headfits(spFramefile,ext=0) else $
              spframe_read,spFramefile, hdr=hdr1, plugmap=plugmap
            mjd_beg[j] = fxpar(hdr1, 'TAI-BEG')/24.D/3600.D
            mjd_end[j] = fxpar(hdr1, 'TAI-END')/24.D/3600.D
             
            if ~keyword_set(diet) then begin
              get_tai, hdr1, tai_beg, tai, tai_end
              for ifib=0, 499L do result[i].airmass[ifib] = $
                tai2airmass(plugmap[ifib].ra, plugmap[ifib].dec, tai=tai)
            endif
         endif     

      endfor
      result[i].mjd_beg = mjd_beg
      result[i].mjd_end = mjd_end
   endfor
  ; endif
   mean_mjd=dblarr(4)

   ; identify all the exposures that were combined in spPlate file
   filename = 'spPlate-' + platestr + '-' + mjdstr + '.fits'
   filename = lookforgzip(filepath(filename, root_dir=topdir, $
      subdirectory=[twoddir,platestr]), count=ct)

   if (ct eq 1) then filename = filename[0] else $
	message, 'spPlate file not found: ',plate,mjd

   head = headfits(filename, ext=0)
   nexp_b1 = fxpar(head, 'nexp_b1')
   nexp_b2 = fxpar(head, 'nexp_b2')
   nexp_r1 = fxpar(head, 'nexp_r1')
   nexp_r2 = fxpar(head, 'nexp_r2')
   sn2_g1=fxpar(head, 'SN2EXT1G')
   sn2_g2=fxpar(head, 'SN2EXT2G')
   sn2=min([sn2_g1,sn2_g2])
   result.platesn2=sn2
   result.seeing50=fxpar(head, 'seeing50')
   result.seeing20=fxpar(head, 'seeing20')
   result.seeing80=fxpar(head, 'seeing80')
   result.plate_airmass=fxpar(head, 'airmass')

   nexp = nexp_b1 + nexp_b2 + nexp_r1 + nexp_r2
   exptags = 'EXPID' + string(indgen(nexp)+1, format='(i2.2)')
   for i=0L, nexp - 1L do begin 
     expstr = 'spFrame-' + strmid(fxpar(head, exptags[i]),0,11) + '.fits'
     ind = where(result.name eq expstr)
     exp_no = floor(ind/4.) & col_no = ind - exp_no*4L
     result[exp_no].coadded[col_no] = 1L
   endfor
  
   ind = where(result.coadded[0,*] eq 1L, ncoadd)

   nexp = n_elements(result)
   icam = lonarr(nexp)
   for i=0, nexp - 1 do icam[i] = (where(strmatch(result[i].name, 'spFrame*') ne 0))[0]
   mean_mjd = 0.D
   for i=0, ncoadd-1 do mean_mjd = mean_mjd + $
     0.5*( (result.mjd_beg)[icam[ind[i]],ind[i]] + (result.mjd_end)[icam[ind[i]],ind[i]])
   mean_mjd = mean_mjd / float(ncoadd)
   ;mean_mjd = mean( 0.5*( (result.mjd_beg)[0,ind] + (result.mjd_end)[0,ind]))

   exptime = total(result[ind].exptime,/double)

   expinfo = result
END
