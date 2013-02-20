pro sos_logs_plotsky

   scidat = mrdfits('logfile-all.fits', 4)
   jd = 2400000.5D + scidat.tai / (24.D*3600.D)

   moonpos, jd, moon_ra, moon_dec
   sunpos, jd, sun_ra, sun_dec

   eq2hor, moon_ra, moon_dec, jd, moon_alt, moon_az, obsname='apo'
   eq2hor, sun_ra, sun_dec, jd, sun_alt, sun_az, obsname='apo'

   moondist = djs_diff_angle(scidat.ra, scidat.dec, moon_ra, moon_dec)

   igood = where(scidat.exptime GT 800 AND sun_alt LT -18 $
    AND scidat.mjd GE 55000+365 AND scidat.camera EQ 'r1' $
    AND scidat.quality EQ 'excellent')
   splot, moon_alt[igood], scidat[igood].skypersec, psym=3, $
    xtitle='Moon altitude [deg]', ytitle='Sky continuum [counts/sec]', $
    title='BOSS Sky Brightness (r1) excluding OH', charsize=1.5
; Plot different moon phases in different colors ???

stop
   return
end
