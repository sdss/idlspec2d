;------------------------------------------------------------------------------
pro plot_velcompare1, name, zans, qq, i1, i2, _EXTRA=EXTRA, $
 vdiff=vdiff, vchi=vchi, rmag=rmag

   cspeed = 2.99792458d5

   j = where(qq[i1] AND qq[i2])
   v1 = zans[i1[j]].z * cspeed
   v2 = zans[i2[j]].z * cspeed
   v1_err = zans[i1[j]].z_err * cspeed
   v2_err = zans[i2[j]].z_err * cspeed
   vdiff = v1 - v2
   vchi = vdiff / sqrt(v1_err^2 + v2_err^2)
   rmag = 22.5 - 2.5*alog10(zans[i1[j]].spectroflux[2]>1)
   plothist, v1-v2, title=name, _EXTRA=EXTRA
   splog, name + ': Number of unique pairs = ', n_elements(j)/2
   splog, name + ': Stddev = ', djsig(v1-v2), ' (sigma-clipped)'
   splog, name + ': Median |chi| = ', median(abs(vchi))
   splog, name + ': Median (chi^2) = ', median(vchi^2)

   return
end
;------------------------------------------------------------------------------
pro plot_velcompare

   thick = 3
   csize = 2

   dtheta = 1./3600

   ; Read only the plate 412 repeats
;   plate = [412,412,412,412,412,412]
;   mjd = [51871,51931,52235,52250,52254,52258]
;   readspec, plate, mjd=mjd, zans=zans, plug=plug

   ; Read all plates
   columns = ['CLASS','Z','Z_ERR','ZWARNING', $
    'SPECTROFLUX','PLUG_RA','PLUG_DEC', $
    'PRIMTARGET','SECTARGET']
   readpublic, columns=columns, zans=zans, plug=plug

   ingroup = spheregroup(zans.plug_ra, zans.plug_dec, dtheta, $
    multgroup=multgroup, firstgroup=firstgroup, nextgroup=nextgroup)

   nobj = long(total(firstgroup GE 0))
   nobs = max(multgroup)
   indx1 = lonarr(nobj,nobs,nobs)
   indx2 = lonarr(nobj,nobs,nobs)

   for iobj=0L, nobj-1L do begin
      jlist = lonarr(nobs)-1L
      jlist[0] = firstgroup[iobj]
      for i=0L, multgroup[iobj]-2L do jlist[i+1] = nextgroup[jlist[i]]
      indx1[iobj,*,*] = rebin(jlist,nobs,nobs)
      indx2[iobj,*,*] = transpose(indx1[iobj,*,*])
   endfor

   ipair = where(indx1 GE 0 AND indx2 GE 0 AND indx1 NE indx2)

   qstar = strmatch(zans.class,'STAR*') AND zans.zwarning EQ 0 $
    AND zans.z_err GT 0
   qgal = strmatch(zans.class,'GALAXY*') AND zans.zwarning EQ 0 $
    AND zans.z_err GT 0

   qmain = (plug.primtarget AND 2L^6 + 2L^7 + 2L^8) NE 0 $
    AND strmatch(zans.class,'GALAXY*')
   qlrg = (plug.primtarget AND 2L^5 + 2L^26) NE 0 $
    AND strmatch(zans.class,'GALAXY*') AND (qmain EQ 0)
   qfstar = (plug.sectarget AND 2L^1 + 2L^5) NE 0 $
    AND strmatch(zans.class,'STAR*')

   dfpsplot, 'velcompare.ps', /square, /color

   plot_velcompare1, 'F stars', zans, qfstar, indx1[ipair], indx2[ipair], $
    charsize=csize, cthick=thick, thick=thick, $
    xrange=[-100,100], vdiff=vdiff_fstar, vchi=vchi_fstar, rmag=rmag_fstar
   plot_velcompare1, 'Main galaxies', zans, qmain, indx1[ipair], indx2[ipair], $
    charsize=csize, cthick=thick, thick=thick, $
    xrange=[-100,100], vdiff=vdiff_main, vchi=vchi_main, rmag=rmag_main
   plot_velcompare1, 'LRGs', zans, qlrg, indx1[ipair], indx2[ipair], $
    charsize=csize, cthick=thick, thick=thick, $
    xrange=[-100,100], vdiff=vdiff_lrg, vchi=vchi_lrg, rmag=rmag_lrg

   djs_plot, rmag_main, vdiff_main, psym=3, $
    xrange=[14,22], yrange=[-200,200], $
    xtitle=textoidl('r_{spectro} [mag]'), $
    ytitle='Velocity difference [km/s]', $
    charsize=csize, cthick=thick, thick=thick
   djs_oplot, rmag_lrg, vdiff_lrg, psym=4, $
    symsize=0.25, color='red', $
    charsize=csize, cthick=thick, thick=thick
   djs_oplot, rmag_fstar, vdiff_fstar, psym=4, $
    symsize=0.25, color='blue', $
    charsize=csize, cthick=thick, thick=thick

   djs_plot, rmag_main, abs(vdiff_main), psym=3, $
    xrange=[14,22], yrange=[0.1,1000], /ylog, $
    xtitle=textoidl('r_{spectro} [mag]'), $
    ytitle='Velocity difference [km/s]', $
    charsize=csize, cthick=thick, thick=thick
   djs_oplot, rmag_lrg, abs(vdiff_lrg), psym=4, $
    symsize=0.25, color='red', $
    charsize=csize, cthick=thick, thick=thick
   djs_oplot, rmag_fstar, abs(vdiff_fstar), psym=4, $
    symsize=0.25, color='blue', $
    charsize=csize, cthick=thick, thick=thick

   djs_plot, rmag_main, abs(vchi_main), psym=3, $
    xrange=[14,22], yrange=[0.01,100], /ylog, $
    xtitle=textoidl('r_{spectro} [mag]'), $
    ytitle='\chi', $
    charsize=csize, cthick=thick, thick=thick
   djs_oplot, rmag_lrg, abs(vchi_lrg), psym=4, $
    symsize=0.25, color='red', $
    charsize=csize, cthick=thick, thick=thick
   djs_oplot, rmag_fstar, abs(vchi_fstar), psym=4, $
    symsize=0.25, color='blue', $
    charsize=csize, cthick=thick, thick=thick

   dfpsclose
save, file='velcompare.ss'

   vdiff_all = abs([vdiff_main, vdiff_lrg, vdiff_fstar])
   vchi_all = abs([vchi_main, vchi_lrg, vchi_fstar])
   rmag_all = [rmag_main, rmag_lrg, rmag_fstar]
   isort = sort(rmag_all)
   vdiff_all = vdiff_all[isort]
   vchi_all = vchi_all[isort]
   rmag_all = rmag_all[isort]
   rmin = min(rmag_all, max=rmax)
   dmag = 0.20
   nbin = fix((rmax - rmin)/dmag)
   rmag_vec = rmin + findgen(nbin+1) * 0.1
   vdiff_vec = fltarr(nbin)
   vchi_vec = fltarr(nbin)
   for j=0L, nbin-1L do begin
      ii = where(rmag_all GE rmag_vec[j] AND rmag_all LT rmag_vec[j+1], ct)
      if (ct GT 1) then begin
         jsort = sort(vdiff_all[ii])
         vdiff_vec[j] = vdiff_all[ii[jsort[0.67*ct]]]
         jsort = sort(vchi_all[ii])
         vchi_vec[j] = vchi_all[ii[jsort[0.67*ct]]]
      endif
   endfor

stop

   return
end
;------------------------------------------------------------------------------
