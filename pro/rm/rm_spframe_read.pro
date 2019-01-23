;+
; NAME
;
; PURPOSE
;
;
;
;

pro rm_spframe_read,objname,indx,loglam=loglam,objflux=objflux,objivar=objivar $
      , plate=plate

   topdir = getenv('BOSS_SPECTRO_REDUX')
   twoddir = getenv('RUN2D')
   if not keyword_Set(plate) then plate = 7338
   platestr = string(plate,format='(i4.4)')

   filename = lookforgzip(filepath(objname, root_dir=topdir, $
      subdirectory=[twoddir,platestr]), count=ct)
   if ct eq 1 then filename=filename[0]


   if not keyword_set(adderr) then adderr = 0.03
   spframe_read,filename,indx,objflux=objflux, objivar=objivar $
    , wset=wset, loglam=loglam,dispimg=dispimg, adderr=adderr

   ; Make a map of the size of each pixel in delta-(log10-Angstroms).
   ; Re-normalize the flux to ADU/(dloglam).
   ; Re-normalize the dispersion from /(raw pixel) to /(new pixel).
   correct_dlam, objflux, objivar, wset, dlam=dloglam
   correct_dlam, dispimg, 0, wset, dlam=dloglam, /inverse


end
