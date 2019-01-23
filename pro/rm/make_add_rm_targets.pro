; make a target list of additional RM quasars for
; the MMT/Hectospec

pro make_add_rm_targets

file='/home/yshen/products/Linux/idlrm/etc/RM_targets.fits'
target=mrdfits(file,1)

; keep those not receiving a SDSS fiber
ind=where(target.tiled eq 0,nnn)
target=target[ind]

fmt='(a6, " ", f11.6, " ", f11.6, " ", f5.2, " ", f6.4)'
outfile='/home/yshen/products/Linux/idlrm/etc/add_RM_targets'
openw, lun, outfile,/get_lun
printf, lun, '## Additional RM targets for MMT/Hectospec'
printf, lun, '## i<21.7 but did not receive a SDSS fiber'
printf, lun, '## ID  RA[deg]  DEC[deg]  i-mag  redshift'

for i=0L, nnn-1 do begin

   printf, lun, format=fmt, 'RM'+string(1000L + i, format='(i4.4)'), $
       target[i].ra, target[i].dec, (target[i].psfmag)[3], target[i].z

endfor

close, lun
free_lun, lun

end
