pro diffplate, platenum, mjd=mjd

platenum = 306
mjd = [51637, 51690]

;platenum = 302
;mjd = [51616,51688]

   readspec, platenum, mjd=mjd[0], flux=flux1, flerr=flerr1, wave=wave1, $
    plugmap=plug1
;flux1=flux1[2000:3800,*]
;flerr1=flerr1[2000:3800,*]

   readspec, platenum, mjd=mjd[1], flux=flux2, flerr=flerr2, wave=wave2, $
    plugmap=plug2
;flux2=flux2[2000:3800,*]
;flerr2=flerr2[2000:3800,*]


   nfiber = (size(flux1,/dimens))[1]
   result = veldisp_struc(nfiber)
   res1 = result[0]

   for ifiber=0, nfiber-1 do begin
;   for ifiber=0, 50 do begin ; ???
print, 'FIBER', ifiber+1
      adist = djs_diff_angle(plug1[ifiber].ra, plug1[ifiber].dec, $
       plug2.ra, plug2.dec)
      dmin = min(adist, imin)
      if (dmin LE 1./3600 AND max(flux1[*,ifiber]) GT 0 AND $
       max(flux2[*,imin]) GT 0) then begin
         ; Match fiber IFIBER in first plate with fiber IMIN in second
         veldisp, flux1[*,ifiber], flerr1[*,ifiber], wave1[*,ifiber], $
          flux2[*,imin], flerr2[*,imin], wave2[*,ifiber], res1, /nodiff, /nobe
         copy_struct_inx, res1, result, index_from=0, index_to=ifiber
      endif
   endfor

md=djs_median(flux1,1)
stop

set_plot,'ps'
device,file='diff-306-old-red.ps'
plot,alog10(md),result.z*3e5,ps=2,xr=[0,3],yr=400*[-1,1],$
 xtitle='log(Med-Flux)', ytitle='delta-z'
oplot,[0,6],[0,0]-0
oplot,[0,6],[0,0]-25
device,/close
set_plot,'x'

i=where(md GT 10^1.8)
j=where(result.z*3e5 GT -400 AND result.z*3e5 LT 400)
print,stddev(result[i].z*3e5),stddev(result[j].z*3e5)

i=where(plug1.primtarget AND 64)
djs_oplot,alog10(md[i]),result[i].z*69.,ps=2,xr=[0,3],yr=500*[-1,1],$
 color='red'

   return
end
