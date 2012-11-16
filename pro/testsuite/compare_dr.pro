; Compare redshifts in common between two sets of reductions
pro compare_dr_loz, sp1, sp2, indx, qgood1, qgood2
   indx = where((sp1.boss_target1 AND $
    sdss_flagval('boss_target1','GAL_LOZ')) NE 0)
   qgood1 = sp1[indx].zwarning_noqso EQ 0
   qgood2 = sp2[indx].zwarning_noqso EQ 0
   return
end
pro compare_dr_cmass, sp1, sp2, indx, qgood1, qgood2
   indx = where((sp1.boss_target1 AND $
    sdss_flagval('boss_target1','GAL_CMASS')) NE 0)
   qgood1 = sp1[indx].zwarning_noqso EQ 0
   qgood2 = sp2[indx].zwarning_noqso EQ 0
   return
end
pro compare_dr_qso, sp1, sp2, indx, qgood1, qgood2
   indx = where(strtrim(sp1.objtype,2) EQ 'QSO')
   qgood1 = sp1[indx].zwarning EQ 0
   qgood2 = sp2[indx].zwarning EQ 0
   return
end
pro compare_dr, run1=run1, run2=run2
   if (NOT keyword_set(run1)) then run1 = 'v5_4_45' ; for DR9
   if (NOT keyword_set(run2)) then run2 = 'v5_5_12' ; for DR10

   columns = ['plate','mjd','fiberid', $
    'boss_target1','boss_target2','ancillary_target1', $
    'objtype','class','subclass','class_noqso','subclass_noqso', $
    'z','zwarning','z_noqso','zwarning_noqso']

   ; Read the spAll files, and match them
   file1 = filepath('spAll-'+run1+'.fits', $
    root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run1)
   file2 = filepath('spAll-'+run2+'.fits', $
    root_dir=getenv('BOSS_SPECTRO_REDUX'), subdir=run2)
   sp1 = hogg_mrdfits(file1, 1, nchunk=20, columns=columns)
   sp2 = hogg_mrdfits(file2, 1, nchunk=20, columns=columns)
   n1 = n_elements(sp1)
   n2 = n_elements(sp2)
   iname1 = strarr(n1)
   for i=0L, n1-1L do iname1[i] = $
    string(sp1[i].plate, sp1[i].mjd, sp1[i].fiberid, format='(i6.6,i5.5,i4.4)')
   iname2 = strarr(n2)
   for i=0L, n2-1L do iname2[i] = $
    string(sp2[i].plate, sp2[i].mjd, sp2[i].fiberid, format='(i6.6,i5.5,i4.4)')
   isort1 = sort(iname1)
   sp1 = sp1[isort1]
   iname1 = iname1[isort1]
   isort2 = sort(iname2)
   sp2 = sp2[isort2]
   iname2 = iname2[isort2]
   ilist1 = lonarr(n1) - 1L
   ilist2 = lonarr(n1) - 1L
   i1 = 0L
   i2 = 0L
   while (i1 LT n1) do begin
      if (iname1[i1] EQ iname2[i2]) then begin
         ilist1[i1] = i1
         ilist2[i1] = i2
         i1++
         i2++
      endif else begin
         if (iname1[i1] LT iname2[i2]) then i1++ $
          else i2++
      endelse
   endwhile
   ikeep = where(ilist1 NE -1)
   sp1 = sp1[ilist1[ikeep]]
   sp2 = sp2[ilist2[ikeep]]

   vdiff = (sp1.z - sp2.z) * 3e5
   vdiff_noqso = (sp1.z_noqso - sp2.z_noqso) * 3e5

   compare_dr_loz, sp1, sp2, indx, qgood1, qgood2
   igood = indx[where(qgood1*qgood2)]
   splog, ''
   splog, n_elements(indx), ' LOWZ targets'
   splog, mean(qgood1*qgood2), ' fraction good'
   splog, median(vdiff_noqso[igood]), ' median vdiff'
   splog, djsig(vdiff_noqso[igood]), ' rms vdiff'
   splog, mean(abs(vdiff_noqso[igood]) GT 300), ' fraction changed redshift'
   splog, mean(qgood1 EQ 0 AND qgood2 EQ 1), ' fraction good -> bad'
   splog, mean(qgood2 EQ 0 AND qgood1 EQ 1), ' fraction bad -> good'

   compare_dr_cmass, sp1, sp2, indx, qgood1, qgood2
   igood = indx[where(qgood1*qgood2)]
   splog, ''
   splog, n_elements(indx), ' CMASS targets'
   splog, mean(qgood1*qgood2), ' fraction good'
   splog, median(vdiff_noqso[igood]), ' median vdiff'
   splog, djsig(vdiff_noqso[igood]), ' rms vdiff'
   splog, mean(abs(vdiff_noqso[igood]) GT 300), ' fraction changed redshift'
   splog, mean(qgood1 EQ 0 AND qgood2 EQ 1), ' fraction good -> bad'
   splog, mean(qgood2 EQ 0 AND qgood1 EQ 1), ' fraction bad -> good'

   compare_dr_qso, sp1, sp2, indx, qgood1, qgood2
   igood = indx[where(qgood1*qgood2)]
   splog, ''
   splog, n_elements(indx), ' QSO targets'
   splog, mean(qgood1*qgood2), ' fraction good'
   splog, median(vdiff[igood]), ' median vdiff'
   splog, djsig(vdiff[igood]), ' rms vdiff'
   splog, mean(abs(vdiff[igood]) GT 1500), ' fraction changed redshift'
   splog, mean(qgood1 EQ 0 AND qgood2 EQ 1), ' fraction good -> bad'
   splog, mean(qgood2 EQ 0 AND qgood1 EQ 1), ' fraction bad -> good'

   return
end
