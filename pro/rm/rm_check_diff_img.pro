; suite to check the difference imaging LCs

; check the error of the std
pro check_std_err, band=band

  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  epoch_info=mrdfits(file,2)
  
  if ~keyword_set(band) then band = 'g'

  result=replicate({rmid:0L, bok_err2stddev:0.D, bok_err2mad:0.D}, 70)
  result.rmid=indgen(70)+850

  ; where the photometric LCs are stored
  topdir='/data3/yshen/ftp/sdssrm/collab/photo_lc/img_diff_lc/by_rmid/'

  figfile='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/bok/std_LC.ps'
  begplot, name=figfile, /color, /landscape

  colors=cgcolor(['black', 'cyan', 'red', 'magenta', 'orange', 'purple'])
  for i=0L, 69L do begin

     rmid=result[i].rmid
     file=topdir+'rm'+string(rmid,format='(i3.3)')+'.fits'
     lc=mrdfits(file,1,/silent)

     ind=where(strmatch(lc.orig,'bok*') and lc.band eq band )
     if ind[0] ne -1 then begin
        lc=lc[ind]
        orig=lc.orig
        uniq_tag=orig[uniq(orig,sort(orig))]
        ntag=n_elements(uniq_tag)
        err2stddev=dblarr(ntag) & err2mad=dblarr(ntag)
        stddev_arr=dblarr(ntag)
        mad_arr=dblarr(ntag) & mederr_arr=dblarr(ntag)
        title='RMID '+string(rmid,format='(i3.3)') + ' g-band'
        for j=0L, ntag - 1 do begin
           indd=where(orig eq uniq_tag[j])
           stddev_arr[j]=stddev(lc[indd].flux)
           mad_arr[j]=median( abs(lc[indd].flux - median(lc[indd].flux) )  )
           mederr_arr[j]=median(lc[indd].flerr)

           if j eq 0 then ploterror, lc[indd].mjd, (lc[indd].flux - median(lc[indd].flux) )/mad_arr[j], lc[indd].flerr/mad_arr[j], $
             xtitle='MJD', ytitle='Flux / MAD (diff img)', psym=4, color=colors[j], errcolor=colors[j],noerase=noerase,title=title, xtickformat='(i5)' $
           else oploterror, lc[indd].mjd, ( lc[indd].flux - median(lc[indd].flux) )/mad_arr[j], lc[indd].flerr/mad_arr[j], $
              psym=4, color=colors[j], errcolor=colors[j]
           oplot, [56600., 56900.], [0.,0.], line=1
        endfor
        legend, uniq_tag, box=0,pos=[0.75,0.92],/norm, color=colors[0:ntag-1],textcolor=colors[0:ntag-1]
       

        ; now normalize median measurement errors to stddev and mad
        err2stddev=mederr_arr/stddev_arr
        err2mad=mederr_arr/mad_arr
        xyouts, 0.15, 0.9, textoidl('(median flerr)/MAD'),/norm
        legend, string(err2mad, format='(f0.3)'), box=0, pos=[0.15, 0.88], textcolor=colors[0:ntag-1],/norm

        result[i].bok_err2stddev=median(err2stddev)
        result[i].bok_err2mad=median(err2mad)

         
     endif


  endfor

  outfile='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/bok/errtest.fits'
  mwrfits, result, outfile, /create

  endplot
  cgfixps, figfile
end

; check the diff img on the fainter SDSS stars generated by Ian
pro compile_faint_star_lc

  ; read the star catalog
  file='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/sdss_star/sdss.fits'
  sdss=mrdfits(file,1)

  band=['g', 'i']

  ; get the diff img
  topdir_std_bok='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/bok/'
  fld=string(indgen(18),format='(i2.2)') & nfld=18

  for iband=0,1 do begin
    for ipos=1,4 do begin
       for ifld=0,nfld-1 do begin
         datdir=topdir_std_bok+strupcase(band[iband])+fld[ifld] $
           +'stnd_ccd'+string(ipos,format='(i0)')+'/'
         orig='bok_'+fld[ifld]+'_ccd'+string(ipos,format='(i0)')

         files=file_search(datdir+'stnd*.dat',count=nfound)
         len=strlen(datdir)
         stdid=long(strmid(strmid(files,len),4, 5))
         nobj=n_elements(stdid)
         for jj=0, nobj-1 do begin
           readcol, files[jj], format='d,d,d',mjd,flux,flerr,/silent
           ntmp=n_elements(mjd)
           if n_elements(mjd_arr) eq 0 then mjd_arr=mjd else mjd_arr=[mjd_arr,mjd]
           if n_elements(flux_arr) eq 0 then flux_arr=flux else flux_arr=[flux_arr,flux]
           if n_elements(flerr_arr) eq 0 then flerr_arr=flerr else flerr_arr=[flerr_arr,flerr]
           if n_elements(stdid_arr) eq 0 then stdid_arr=replicate(stdid[jj],ntmp) $
             else stdid_arr=[stdid_arr, replicate(stdid[jj],ntmp)]
           if n_elements(band_arr) eq 0 then band_arr=replicate(band[iband],ntmp) $
             else band_arr=[band_arr, replicate(band[iband],ntmp)]
           if n_elements(orig_arr) eq 0 then orig_arr=replicate(orig,ntmp) $
             else orig_arr=[orig_arr, replicate(orig,ntmp)]
         endfor
       endfor
    endfor
  endfor
  ; flip the sign of the flux for the original diff img LC
  flux_arr = -flux_arr
  ntot=n_elements(stdid_arr)
  alldata={stdid:0L, mjd:0.D, flux:0.D, flerr:0.D, band:'', orig:''}
  alldata=replicate(alldata, ntot)
  alldata.stdid=stdid_arr & alldata.mjd=mjd_arr & alldata.flux=flux_arr
  alldata.flerr=flerr_arr & alldata.band=band_arr & alldata.orig=orig_arr

  if ~keyword_set(outdir) then outdir='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/'
  outfile=outdir + 'all_lc_data.fits'
  mwrfits, alldata, outfile, /create

end

pro check_faint_star_lc

  ; read the star input catalog
  file='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/sdss_star/sdss.fits'
  sdss=mrdfits(file,1)
  nnn=n_elements(sdss)

  ; read in the LC catalog
  file='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/all_lc_data.fits'
  alldata=mrdfits(file,1)

  ; loop over all stars
  outfile='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/lc_stat.fits'
  if file_test(outfile) eq 0 then begin
  for i=0L, nnn-1 do begin
    ind=where(alldata.stdid eq i)
    if ind[0] gt 0 then begin
      sub=alldata[ind]
      orig=sub.orig
      uniq_tag=orig[uniq(orig,sort(orig))]
      nset=n_elements(uniq_tag)
      for jj=0l, nset - 1 do begin

        ; for gband
        indd=where(sub.orig eq uniq_tag[jj] and sub.band eq 'g')
        if indd[0] gt 0 then begin
           if n_elements(orig_arr) eq 0 then orig_arr=uniq_tag[jj] else orig_arr=[orig_arr,uniq_tag[jj]]
           if n_elements(band_arr) eq 0 then band_arr='g' else band_arr=[band_arr, 'g']
           if n_elements(stdid_arr) eq 0 then stdid_arr=i else stdid_arr=[stdid_arr, i]
           mad=median( abs( sub[indd].flux - median(sub[indd].flux) )  )
           mederr=median(sub[indd].flerr)
           if n_elements(mad_arr) eq 0 then mad_arr=mad else mad_arr=[mad_arr, mad]
           if n_elements(mederr_arr) eq 0 then mederr_arr=mederr else mederr_arr=[mederr_arr, mederr]
           if n_elements(mag_arr) eq 0 then mag_arr=sdss[i].g else mag_arr=[mag_arr,sdss[i].g]
        endif
        ; for iband
        indd=where(sub.orig eq uniq_tag[jj] and sub.band eq 'i')
        if indd[0] gt 0 then begin
           if n_elements(orig_arr) eq 0 then orig_arr=uniq_tag[jj] else orig_arr=[orig_arr,uniq_tag[jj]]
           if n_elements(band_arr) eq 0 then band_arr='i' else band_arr=[band_arr, 'i']
           if n_elements(stdid_arr) eq 0 then stdid_arr=i else stdid_arr=[stdid_arr, i]
           mad=median( abs( sub[indd].flux - median(sub[indd].flux) )  )
           mederr=median(sub[indd].flerr)
           if n_elements(mad_arr) eq 0 then mad_arr=mad else mad_arr=[mad_arr, mad]
           if n_elements(mederr_arr) eq 0 then mederr_arr=mederr else mederr_arr=[mederr_arr, mederr]
           if n_elements(mag_arr) eq 0 then mag_arr=sdss[i].i else mag_arr=[mag_arr,sdss[i].i]
        endif
      endfor
    endif
  endfor

  ntot=n_elements(stdid_arr)
  result={stdid:0L, band:'', orig:'', mag:0.D, mad:0.D, mederr:0.D}
  result=replicate(result, ntot)
  result.stdid=stdid_arr & result.band=band_arr & result.orig=orig_arr
  result.mag=mag_arr & result.mad=mad_arr & result.mederr=mederr_arr
  outfile='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/lc_stat.fits'
  mwrfits, result, outfile, /create
 endif

  result=mrdfits(outfile, 1)
  band=['g', 'i']
  for i=0, 1 do begin
    figfile='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/std_lc_check_'+band[i]+'.ps'
    begplot, name=figfile, /landscape, /color
      ind=where(result.band eq band[i] and result.mederr gt 0)
      plot, result[ind].mag, result[ind].mad/result[ind].mederr, psym=3, xtitle='Magnitude', $
         ytitle='MAD / MedErr', title = 'diff img (faint sdss stars) ' + band[i], /ylog, yrange=[1d0, 1d3]
    endplot
    cgfixps, figfile
  endfor
end
