; Select /RSAMPLE to random-sample in each redshift bin
;   rather than taking the 10 brightest.
pro danvb_qsolist, rsample=rsample

   danfilename = 'plateEDR.summ.fits'
   outfile = 'eigeninput_qso.dat'

   ; Make the list of plates that are unique survey-quality.
   ; Ignore plates at MJD<=51605 which had electronics problems.
   ; Include only plates at MJD >= 51789 which have better spectro-photometry.
   platelist, plist=plist
   plist = plist[where(plist.qsurvey AND plist.mjd GE 51789)]

   ; Read Van den Berk's file
   danvb = mrdfits(danfilename,1)

   ; Trim Van den Berk's file to only objects on the unique survey-quality
   ; plates
   qgood = bytarr(n_elements(danvb))
   for i=0, n_elements(danvb)-1 do $
    if ((where(danvb[i].plate EQ plist.plate))[0] NE -1) then qgood[i] = 1
   danvb = danvb[where(qgood)]

   ; Trim to manually-classified QSO's, ignoring where he set FLAG_MANUAL=1.
   ; Note that he has at least two incorrect identifications:
   ;   275/51910-353, 289/51990-313, both of which are really stars.
   ; Object 412/51931-263 is not a z=5.8441 QSO.
   ; Object 367/51997-506 may not be a z=5.298 QSO.
   indx = where(strtrim(danvb.id_manual) EQ 'QSO' $
    AND danvb.z_manual GT 0.01 AND danvb.flag_manual NE 1)
   danvb = danvb[indx]

   ; Read the P-1D outputs
   readspec, danvb.plate, danvb.fiber, mjd=danvb.mjd, zans=zans

   ; Do one of the following:
   ; (1) Random-sample this list to 10 objects in each delta-z=0.1 bin
   ; (2) Or select the 10 brightest objects (according to ZANS.SN_MEDIAN)
   ;   in each delta-z=0.1 bin.
   qgood = bytarr(n_elements(danvb))
   deltaz = 0.1
   nperbin = 10
   iseed = 1234L
   for ibin=0, long(7.0/deltaz) do begin
      zlo = deltaz * ibin
      zhi = deltaz * (ibin+1)
      ii = where(danvb.z_manual GE zlo AND danvb.z_manual LT zhi, ni)
      if (ni GE 1 and ni LE nperbin) then begin
         qgood[ii] = 1
      endif else if (ni GT nperbin) then begin
         if (keyword_set(rsample)) then begin
            rr = randomu(iseed, ni)
            qgood[ii[ (sort(rr))[0:nperbin-1] ]] = 1
         endif else begin
            thissn = zans[ii].sn_median
            qgood[ii[ (sort(-thissn))[0:nperbin-1] ]] = 1
         endelse
      endif
   endfor

   indx = where(qgood)
   danvb = danvb[indx]
   zans = zans[indx]

   outstruct = create_struct( $
    'plate'       ,  0L, $
    'mjd'         ,  0L, $
    'fiberid'     ,  0L, $
    'z'           ,  0D, $ ; Cast as double so more digits written out
    'comment'     ,  ' ' )
   outstruct = replicate(outstruct, n_elements(danvb))
   outstruct.plate = danvb.plate
   outstruct.mjd = danvb.mjd
   outstruct.fiberid = danvb.fiber
   outstruct.z = danvb.z_manual
   outstruct.comment = string(zans.sn_median, format='("S/N=",f5.2)')
   ii = where(danvb.flag_manual EQ 2)
   if (ii[0] NE -1) then $
    outstruct[ii].comment = outstruct[ii].comment + ' BAL'

   struct_print, outstruct, filename=outfile

end
