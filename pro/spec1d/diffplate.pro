pro diffplate, platenum, mjd=mjd

platenum = 306
mjd = [51637, 51690]

;platenum = 302
;mjd = [51616,51688]

platenum = 406
mjd = [51817, 51869]

   readspec, platenum, mjd=mjd[0], flux=flux1, flerr=flerr1, wave=wave1, $
    plugmap=plug1
;flux1=flux1[2000:3800,*]
;flerr1=flerr1[2000:3800,*]

   readspec, platenum, mjd=mjd[1], flux=flux2, flerr=flerr2, wave=wave2, $
    plugmap=plug2
;flux2=flux2[2000:3800,*]
;flerr2=flerr2[2000:3800,*]

nfiber=640
   for ifiber=0, nfiber-1 do begin
;   for ifiber=0, 50 do begin ; ???
print, 'FIBER', ifiber+1
      adist = djs_diff_angle(plug1[ifiber].ra, plug1[ifiber].dec, $
       plug2.ra, plug2.dec)
      dmin = min(adist, imin)
      if (dmin LE 1./3600 AND max(flux1[*,ifiber]) GT 0 AND $
       max(flux2[*,imin]) GT 0) then begin
         ; Match fiber IFIBER in first plate with fiber IMIN in second
         zoffset = alog10(wave1[0,ifiber] / wave2[0,ifiber]) / 1.d-4
         res1 = veldisp(flux1[*,ifiber], flerr1[*,ifiber], $
          flux2[*,imin], flerr2[*,imin], zoffset=zoffset)
         if (NOT keyword_set(result)) then begin
            result = replicate(res1, nfiber)
         endif
         copy_struct_inx, res1, result, index_from=0, index_to=ifiber
      endif
   endfor

md=djs_median(flux1,1)
stop

set_plot,'ps'
device,file='diff-306b.ps'
djs_plot,alog10(result.zconf),result.z*70.,ps=1,syms=0.5,yr=100*[-1,1],$
 xtitle='log_{10}(zconf)', ytitle='\Delta z [km/s]', charsize=2, $
 title='Plate 306 MJD 51637 vs 51690', /ystyle
oplot,[-5,2],[0,0]-0
device,/close
set_plot,'x'

i=where(result.zconf LT 0.3)
j=where(result.z*70 GT -400 AND result.z*70 LT 400)
print,mean(result[i].z*70),mean(result[j].z*70)
print,median(result[i].z*70),median(result[j].z*70)
print,stddev(result[i].z*70),stddev(result[j].z*70)

i=where(plug1.primtarget AND 64)
djs_oplot,alog10(md[i]),result[i].z*69.,ps=2,xr=[0,3],yr=500*[-1,1],$
 color='red'

   return
end
