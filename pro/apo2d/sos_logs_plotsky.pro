pro sos_logs_plotsky

   scidat = mrdfits('logfile-all.fits', 4)
   jd = 2400000.5D + scidat.tai / (24.D*3600.D)

   moonpos, jd, moon_ra, moon_dec
   sunpos, jd, sun_ra, sun_dec

   eq2hor, moon_ra, moon_dec, jd, moon_alt, moon_az, obsname='apo'
   eq2hor, sun_ra, sun_dec, jd, sun_alt, sun_az, obsname='apo'

   ; Compute the moon phase -- remap from the fraction of moon illuminated
   mphase, jd, moonfrac
   jd_tmp = 2450280.2d0 + 14.76*dindgen(1000)/1000.
   moonphase_tmp = dindgen(1000)/1000. ; spanning 0-1 for half a lunar cycle
   mphase, jd_tmp, moonfrac_tmp
   linterp, moonfrac_tmp, moonphase_tmp, moonfrac, moonphase

   moondist = djs_diff_angle(scidat.radeg, scidat.decdeg, moon_ra, moon_dec)

   dfpsplot, 'boss-sky-b2.eps', /square, /color
   qq = scidat.exptime GT 100 $
    AND scidat.mjd GE 55000+365 $
    AND scidat.camera EQ 'b2' $
    AND sun_alt LT -18
   ii = where(qq)
   csize = 1.6
   djs_plot, moon_alt[ii], scidat[ii].skypersec, psym=1, symsize=0.25, $
    xtitle='Moon altitude [deg]', ytitle='Sky continuum [counts/sec/pix]', $
    title='BOSS Sky Brightness (g-band) continuum', charsize=csize, $
    xrange=[-90,90], yrange=[0.01,1], /ylog
   i1 = where(qq AND moon_alt GT 2 AND  moonphase GT 0.25 AND moonphase LT 0.35)
   i2 = where(qq AND moon_alt GT 2 AND  moonphase GT 0.50 AND moonphase LT 0.60)
   djs_oplot, moon_alt[i1], scidat[i1].skypersec, $
    psym=1, color='red', symsize=0.5, /ylog
   djs_oplot, moon_alt[i2], scidat[i2].skypersec, $
    psym=1, color='blue', symsize=0.5, /ylog
   djs_oplot, [0,30], 0.055*[1,2], color='green', /ylog, thick=3
   xyouts, -80, 0.50, 'BLUE: Moon phase = 0.50-0.60', charsize=csize
   xyouts, -80, 0.35, 'RED:  Moon phase = 0.25-0.35', charsize=csize
   dfpsclose

   dfpsplot, 'boss-sky-r2.eps', /square, /color
   qq = scidat.exptime GT 100 $
    AND scidat.mjd GE 55000+365 $
    AND scidat.camera EQ 'r2' $
    AND sun_alt LT -18
   ii = where(qq)
   csize = 1.6
   djs_plot, moon_alt[ii], scidat[ii].skypersec, psym=1, symsize=0.25, $
    xtitle='Moon altitude [deg]', ytitle='Sky continuum [counts/sec/pix]', $
    title='BOSS Sky Brightness (i-band) continuum', charsize=csize, $
    xrange=[-90,90], yrange=[0.05,5], /ylog
   i1 = where(qq AND moon_alt GT 2 AND moonphase GT 0.25 AND moonphase LT 0.35)
   i2 = where(qq AND moon_alt GT 2 AND moonphase GT 0.50 AND moonphase LT 0.60)
   djs_oplot, moon_alt[i1], scidat[i1].skypersec, $
    psym=1, color='red', symsize=0.5, /ylog
   djs_oplot, moon_alt[i2], scidat[i2].skypersec, $
    psym=1, color='blue', symsize=0.5, /ylog
   djs_oplot, [0,45], 0.17*[1,2], color='green', /ylog, thick=3
   xyouts, -80, 2.0, 'BLUE: Moon phase = 0.50-0.60', charsize=csize
   xyouts, -80, 1.4, 'RED:  Moon phase = 0.25-0.35', charsize=csize
   dfpsclose

   dfpsplot, 'boss-moon-dist.eps', /square, /color
   djs_plot, moon_alt[ii], moonphase[ii], psym=1, symsize=0.25, $
    xtitle='Moon alt', ytitle='Moon phase', $
    title='BOSS lunar phases', charsize=csize, $
    xrange=[-90,90], yrange=[0,1]
   djs_oplot, moon_alt[i1], moonphase[i1], psym=1, color='red', symsize=0.5
   djs_oplot, moon_alt[i2], moonphase[i2], psym=1, color='blue', symsize=0.5
   dfpsclose

   return
end
