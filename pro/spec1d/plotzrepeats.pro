; Routine to compare redshift errors between multiple observations
; of the same plate.
pro plotzrepeats, plate, mjd

   if (NOT keyword_set(plate)) then plate = 406
   if (NOT keyword_set(mjd)) then mjd = [51817,51869,51876,51900]

   cspeed = 3.e5
   dtheta = 1.0/3600.

   plotfile = string(plate, format='("zrepeat-",i4.4,".ps")')

   ;----------
   ; Read in all MJDs for this plate

   nmjd = n_elements(mjd)
   if (nmjd LE 1) then begin
      print, 'No duplicate MJDs for this plate'
      return
   endif
   readspec, replicate(plate,nmjd), mjd=mjd, zans=zans, plug=plug

   ;----------
   ; Find all pairs of the same object

   nmatch = djs_angle_match(zans.plug_ra, zans.plug_dec, dtheta=dtheta, $
    mcount=mcount, mindx=mindx, mdist=mdist, mmax=nmjd)

   ii = where(mindx NE -1, npair)
   indx1 = long(ii / nmjd)
   indx2 = mindx[ii]

   ;----------
   ; Loop over each class of object

   !p.multi = [0,1,2]
   csize = 2.0
   dfpsplot, plotfile

   classlist = ['GALAXY', 'STAR', 'QSO']

   for iclass=0, n_elements(classlist)-1 do begin
      jj = where(strmatch(zans[indx1].class, classlist[iclass]+'*') $
       AND strmatch(zans[indx2].class, classlist[iclass]+'*') $
       AND zans[indx1].z_err GT 0 $
       AND zans[indx2].z_err GT 0 $
       AND zans[indx1].zwarning EQ 0 $
       AND zans[indx2].zwarning EQ 0, nj)

      if (classlist[iclass] EQ 'QSO') then yrange=[-2000,2000] $
       else if (classlist[iclass] EQ 'STAR') then yrange = [-200,200] $
       else yrange = [-200,200]

      if (nj GT 1) then begin

         if (nj LT 1000) then psym = 4 $
          else psym = 3

         vdiff = (zans[indx1[jj]].z - zans[indx2[jj]].z) * cspeed
         verr = sqrt(zans[indx1[jj]].z_err^2 + zans[indx2[jj]].z_err^2) $
          * cspeed
         chi = vdiff / verr
         cmed = median(abs(chi))

         print, classlist[iclass]
         abschi = abs(chi)
         chisort = abschi[sort(abschi)]
         chi67 = chisort[fix(0.67*n_elements(chisort))]
         print, 'Median(|chi|) = ', cmed
         print, 'Median(verr) = ', median(verr)
         print, '|Chi| at 67% = ', chi67

         plothist, chi, bin=0.20, xrange=[-8,8], /xstyle, $
          xtitle=textoidl('\chi = (z_1-z_2) / z_{err}'), $
          ytitle='Number', charsize=csize, $
          title=classlist[iclass]+' Redshift Errors from Repeats'

         djs_xyouts, -7, 0.9*!y.crange[1], charsize=0.75*csize, $
          'Median(cz_{err})=' + string(median(verr), format='(f6.1)') + ' km/s'
         djs_xyouts, -7, 0.8*!y.crange[1], charsize=0.75*csize, $
          '\chi at 67% =' + string(chi67, format='(f5.2)')

;         xplot = plug[indx1[jj]].sn_median
         xplot = plug[indx1[jj]].mag[2]
         plot, xplot, vdiff, $
          psym=psym, symsize=0.5, yrange=yrange, charsize=csize, $
;          xtitle='Median S/N per pixel', $
          xtitle='r-mag', xrange=[15,22], /xstyle, $
          ytitle=textoidl('(z_1-z_2) [km/s]'), $
          title=classlist[iclass]+' Redshift Errors from Repeats'
         djs_oplot, !x.crange, [0,0]
         djs_oploterr, xplot, vdiff, $
          yerr=verr
      endif
   endfor

   dfpsclose

end
