pro make2dmerge, planfile

   if NOT keyword_set(planfile) then return

   yanny_read, planfile, pdata, hdr=hdr
   platefile = yanny_par(hdr,'combinefile')
   mergefile = yanny_par(hdr,'mergefile')
   yanny_free, pdata

   splog, 'Reading Plate File '+platefile
   splog, 'Output to this prefix '+mergefile

   spawn, 'mkdir -p 2dmerge'


   if NOT keyword_set(platefile) then return
   if NOT keyword_set(mergefile) then return

   flux       = mrdfits( platefile, 0, hdr)
   ivar       = mrdfits( platefile , 1)
   andmask    = mrdfits( platefile, 2)
   ormask     = mrdfits( platefile, 3)
   dispersion = mrdfits( platefile, 4)
   plugmap    = mrdfits( platefile, 5)
   snvec      = mrdfits( platefile, 6)
   mags       = mrdfits( platefile, 7)

   nfiber = (size(flux))[2]

   for i = 0, nfiber - 1 do begin

     fibername = string(format='(a,i3.3,a)',mergefile,i+1,'.fits')

     err = flux[*,i] * 0.0
     good = where(ivar[*,i] GT 0, ngood)
     if good[0] NE -1 then err[good] = 1.0/sqrt(ivar[good,i])

     thishdr = hdr

     sxaddpar, thishdr, 'OBJID', string(format='(5(i))', plugmap[i].objid), $
       before='NWORDER'
     sxaddpar, thishdr, 'MAG', string(format='(5(f11.3))', plugmap[i].mag), $
       before='NWORDER'
     sxaddpar, thishdr, 'RAOBJ', plugmap[i].ra, 'RA (deg) of object', $
       before='NWORDER'
     sxaddpar, thishdr, 'DECOBJ', plugmap[i].dec, 'DEC (deg) of object', $
       before='NWORDER'
     sxaddpar, thishdr, 'OBJTYPE', plugmap[i].objtype, $
       before='NWORDER'
     sxaddpar, thishdr, 'XFOCAL', plugmap[i].xfocal, $
       before='NWORDER'
     sxaddpar, thishdr, 'YFOCAL', plugmap[i].yfocal, $
       before='NWORDER'
     sxaddpar, thishdr, 'SPECID', plugmap[i].spectrographId, $
       before='NWORDER'
     sxaddpar, thishdr, 'PRIMTARG', plugmap[i].primtarget, $
       before='NWORDER'
     sxaddpar, thishdr, 'SECTARGE', plugmap[i].sectarget, $
       before='NWORDER'
     sxaddpar, thishdr, 'FIBERID', plugmap[i].fiberId, $
       before='NWORDER'

     sn = fltarr(3)
      m = fltarr(3)

     if keyword_set(snvec) then sn=snvec[*,i]
     if keyword_set(mags) then m=mags[*,i]

     sxaddpar, thishdr, 'NGOOD', ngood, 'Number of Good Pixels'
     sxaddpar, thishdr, 'PIXMIN', 0.0, 'Place Holding'
     sxaddpar, thishdr, 'PIXMAX', 2047.0, 'Place Holding'
     sxaddpar, thishdr, 'SN_G',  sn[0], "Median S/N ratio in g'"
     sxaddpar, thishdr, 'MAG_G',  m[0], "Synthetic magnitude in g'"
     sxaddpar, thishdr, 'SN_R',  sn[1], "Median S/N ratio in r'"
     sxaddpar, thishdr, 'MAG_R',  m[1], "Synthetic magnitude in r'"
     sxaddpar, thishdr, 'SN_I',  sn[2], "Median S/N ratio in i'"
     sxaddpar, thishdr, 'MAG_I',  m[2], "Synthetic magnitude in i'"

     sxaddpar, thishdr, 'NAXIS1', n_elements(err)
     sxaddpar, thishdr, 'NAXIS2', 3

     ; 1st HDU is flux, error, and dispersion
     mwrfits, [[flux[*,i]],[err]], fibername, thishdr, /create

     ; 2nd HDU are pixelmasks
     mwrfits, [[andmask[*,i]],[ormask[*,i]]], fibername
    
     ; 3rd HDU are pixelmasks
     mwrfits, dispersion[*,i], fibername
     
   endfor

return
end


