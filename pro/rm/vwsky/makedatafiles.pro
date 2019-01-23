;*** make public data files for weights and espec

PRO makedatafiles

;restore, '/data/cass00a/vw/DR3/espec_sky/platenoise_a.sav' ;from platenoise.pro

;plateerror = error_plate
;continuum = cont
;plateid = plateid

;save, plateerror, continuum, plateid, file='/home/vw/idl/pro/sdss/SKY/PUBLIC/plate_weights.sav'

restore, '/data/cass00a/vw/DR4/espec_sky/platenoise_dr4.sav' ;from platenoise.pro

plateerror = error_plate
continuum = cont
plateid = plateid

save, plateerror, continuum, plateid, file='/data/cass00a/vw/DR4/espec_sky/plate_weights.sav'

stop

common sky_info, lambda,bad_sky,no_sky,minlam,maxlam,disp,espec,nbin,nbin_bad
sky_info,str='bo'

espec = espec
disp = 0.0001
lambda = findgen(nbin)*disp + minlam
pix_sky = bad_sky
pix_nosky = no_sky

save, espec,lambda,pix_sky,pix_nosky, file='/home/vw/idl/pro/sdss/SKY/PUBLIC/espec_OH.sav'

END
