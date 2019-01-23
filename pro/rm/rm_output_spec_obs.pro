; Output the mean mjd and # of coadded spectra in each epoch

pro rm_output_spec_obs, result=result

  target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
  fibermap = mrdfits(target_file,1,/silent)
  plate=fibermap[0].plate & mjd=fibermap[0].mjd
  ind=where(plate gt 0, nnn)
  plate=plate[ind] & mjd=mjd[ind]

  result=replicate({plate:0L,mjd:0L,beg_mjd:0L, end_mjd:0L, mean_mjd:0.D,ncoadd:0L, platesn2:0.D, $
      seeing20:0., seeing50:0., seeing80:0., airmass:0.}, nnn)
  result.plate=plate & result.mjd=mjd

  for i=0, nnn - 1 do begin
     rm_expinfo, plate[i], mjd[i], expinfo=expinfo, mean_mjd=mean_mjd, /diet
     coadd=expinfo.coadded[0,*]
     tt=where(coadd eq 1, ncoadd)
     print, plate[i], mjd[i], mean_mjd, ' Ncoadd='+string(ncoadd,format='(i2)'), $
       ' platesn2=', string(expinfo[0].platesn2,format='(f4.1)'), ' seeing50=', expinfo[0].seeing50
     result[i].mean_mjd=mean_mjd & result[i].ncoadd=ncoadd
     result[i].platesn2=expinfo[0].platesn2 & result[i].seeing50=expinfo[0].seeing50
     result[i].seeing20=expinfo[0].seeing20 & result[i].seeing80=expinfo[0].seeing80
     result[i].beg_mjd=floor(min((expinfo.mjd_beg)[0,tt])) 
     result[i].end_mjd=floor(max((expinfo.mjd_end)[0,tt]))
     result[i].airmass=expinfo[0].plate_airmass
  endfor


end
