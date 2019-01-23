; plot ps1 LCs

pro plot_ps1_lc, rmid

   file = '/data3/yshen/work/PS1/yue_all.fits'
   ps1 = mrdfits(file,1)
   file = '/data3/yshen/ftp/sdssrm/collab/science_data/sample_char/sample_char_refined.fits'
   target=mrdfits(file,1)

   spherematch, target[rmid].ra, target[rmid].dec, ps1.ra, ps1.dec, 0.5/3600.D, match0, match1, dist1

   lc_mjd = ps1[match1].lc_mjd
   lc_mag = ps1[match1].lc_mag
   lc_err = ps1[match1].lc_err

   output = {mjd:lc_mjd, mag:lc_mag, err:lc_err}
   outfile = '/data3/yshen/work/PS1/rmid' + string(rmid, format='(i3.3)') + '_ps1_lc.fits'
   mwrfits, output, outfile, /create

   figfile = '/data3/yshen/work/PS1/rmid' + string(rmid, format='(i3.3)') + '.ps'
   begplot, name=figfile, /landscape
   str='PS1_' + ['g','r', 'i', 'z', 'y']

   if ~keyword_set(xrange) then xrange=[54900, 56600]
   if ~keyword_set(yrange) then yrange=[22.5, 19]
   for i=0, 4 do begin

      mjd = lc_mjd[i,*]
      mag = lc_mag[i,*]
      err = lc_err[i,*]
      ind = where(err gt 0 and mag gt 0 and mag lt 30.)
      mjd = mjd[ind]
      mag = mag[ind]
      err = err[ind]

      ploterror, mjd, mag, err, psym=6, xtitle='MJD', ytitle = str[i]+' [mag]', title='RMID'+string(rmid, format='(i3.3)'), $
        xrange=xrange, /xsty, xtickformat='(I0)',/ysty, yrange=yrange
   endfor

   endplot
   cgfixps, figfile

end
