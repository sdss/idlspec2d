; make table 1 and take 3 for the technical overview paper

pro make_tech_overview_tables

;; make table 3
;file=getenv('IDLRM_DIR')+'/etc/wh-sky-mask.txt'
;fmt='(f6.4, f10.3)'
;readfmt,file,fmt, logwave, wave
;outfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/tech_sum/submission/data_tables/datafile3.txt'
;nnn=n_elements(wave)
;openw, lun, outfile,/get_lun
;for i=0L, nnn-1 do printf, lun, format='(f9.3, " ", f6.4)', wave[i], logwave[i]
;close,lun
;free_lun, lun

; make fits table for table 1
file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
tags=['RA','DEC','zfinal','sourcetype', 'PSFMAG', 'objc_type', 'release', 'plate', 'fiberid', $
  'mjd', 'med_sn']
result=mrdfits(file,1, columns=tags)

output=replicate({RMID:0L, RA:0.D, DEC:0.D, Z:0.D, sourcetype:'', psfmag:fltarr(5), $
   objctype:0L, sample:'', plate:lonarr(32), fiberid:lonarr(32), mjd:lonarr(32), medsn:dblarr(32)}, 1000L)

output.RMID=lindgen(1000L)
output.ra=result.ra & output.dec=result.dec
output.z=result.zfinal
output.sourcetype=result.sourcetype
output.psfmag=result.psfmag
output.objctype=result.objc_type
output.sample=result.release
ind=where(strmatch(output.sample, 'dr12*'))
output[ind].sample='boss'
ind=where(strmatch(output.sample, 'dr7*'))
output[ind].sample='dr7'
output.plate=result.plate
output.fiberid=result.fiberid
output.mjd=result.mjd
output.medsn=result.med_sn

outfile='/home/yshen/Research/Projects/reverberation_mapping/level2/proposal/papers/tech_sum/submission/data_tables/tiled_sample.fits'
mwrfits, output, outfile, /create

end
