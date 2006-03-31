;------------------------------------------------------------------------------
pro plot_velcompare

   cspeed = 2.99792458d5

   plate = [412,412,412,412,412,412]
   mjd = [51871,51931,52235,52250,52254,52258]

   readspec, plate, mjd=mjd, zans=zans, plug=plug

   ingroup = spheregroup(zans.plug_ra, zans.plug_dec, 1./3600, $
    multgroup=multgroup, firstgroup=firstgroup, nextgroup=nextgroup)

   nobj = total(firstgroup GE 0)
   nobs = max(multgroup)
   zarr = fltarr(nobj,nobs)
   nper = lonarr(nobj)
   zmean = fltarr(nobj)
   zall = replicate(zans[0],nobj,nobs)

   for iobj=0L, nobj-1L do begin
      j = firstgroup[iobj]
      for i=0L, multgroup[iobj]-1L do begin
         zarr[iobj,i] = zans[j].z * (zans[j].zwarning EQ 0)
         zall[iobj,i] = zans[j]
         j = nextgroup[j]
      endfor
      igood = where(zarr[iobj,*] NE 0, ngood)
      if (ngood GE 1) then begin
         zmean[iobj] = median(zarr[iobj,igood])
         nper[iobj] = ngood
         zarr[iobj,igood] = zarr[iobj,igood] - zmean[iobj]
      endif
   endfor

   rmag = 22.5 - 2.5*alog10(zall.spectroflux[2]>1)
   istar = where(strmatch(zall.class,'STAR*') and zarr NE 0)
   igal = where(strmatch(zall.class,'GALAXY*') and zarr NE 0)
stop

   dfpsplot, 'velcompare.ps', /square
   dfpsclose

   return
end
;------------------------------------------------------------------------------
