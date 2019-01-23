; manually fix the mismatched fibers in the target_fibermap
; first run map_fiber_target to get a clean slate

pro corr_target_fibermap

;file = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits.07132017'
;file = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits.05112018'
file = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits.07192018'

tt1 = mrdfits(file,1)
tt2 = mrdfits(file,2)

plate=tt1.plate & fiber=tt1.fiberid & mjd=tt1.mjd & sn=tt1.med_sn
plate = plate[*,0] & mjd = mjd[*,0]


corrfile='/data3/yshen/ftp/sdssrm/collab/bossredux/fiber_mismatch_eboss/corr_fiber_eboss_all'
readcol,corrfile,format='l,l,l,l,l',rmid, pp, mm, old_fiber, new_fiber
nnn = n_elements(rmid)

for i=0, nnn-1 do begin

  ind1 = rmid[i]
  ind2 = where( mjd eq mm[i])

  ind3 = where( fiber[ind2, *] eq new_fiber[i])

  sn[ind2, ind1] = sn[ind2, ind3]

  print, fiber[ind2, ind1], new_fiber[i]
  fiber[ind2, ind1] = new_fiber[i]


endfor

tt1.fiberid = fiber
tt1.med_sn = sn

outfile = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
mwrfits, tt1, outfile, /create
mwrfits, tt2, outfile 

end
