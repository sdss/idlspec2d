;+
; NAME:
;   design_sdss3test
;
; PURPOSE:
;   Design test plates for SDSS-3
;
; CALLING SEQUENCE:
;   design_sdss3test, platenum
;
; INPUTS:
;   platenum   - Plate number
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;   design_struct()
;
; REVISION HISTORY:
;   29-Oct-2006  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function design_struct, num

   result = create_struct( $
    name = 'STARDATA', $
    'ra'        ,  0.d, $
    'dec'       ,  0.d, $
    'mag'       , fltarr(5), $
    'holetype'  ,   '', $
    'objtype'   ,   '', $
    'priority'  ,   0L, $
    'primtarget',   0L, $
    'sectarget' ,   0L)
   if (keyword_set(num)) then result = replicate(result, num)

   return, result
end
;------------------------------------------------------------------------------
pro design_sdss3test, platenum

   case platenum of
   2634: begin ; Centered at plate 406
      ; Note there is a bright star near ra=35.49, dec=0.40
      racen = 35.88296
      deccen = 0.1250122
      tilenum = 9549
      runnum = [4263, 4874]
      rerun = [137, 137]
      end
   2637: begin ; Centered at plate 416
      racen =  55.49162
      deccen = 0.01375204
      tilenum = 9550
      runnum = [4136, 4145, 4874]
      rerun = [137, 137, 137]
      end
   else: begin
      print, 'Unknown plate number!'
      return
      end
   endcase

   ;----------
   ; Read all objects from pfObjc files

   objs = 0
   for irun=0, n_elements(runnum)-1 do begin
      fields = sdss_fieldlist(runnum[irun])
      for camcol=1, 6 do begin
         sdss_run2radec, runnum[irun], camcol, fields, rerun=rerun[irun], $
          ra=thisra, dec=thisdec
         indx = where(djs_diff_angle(thisra, thisdec, racen, deccen) $
          LT 1.49, ct)
         if (ct GT 0) then begin
            obj1 = sdss_readobj(runnum[irun], camcol, fields[indx], $
             rerun=rerun[irun])
            if (keyword_set(obj1)) then $
             objs = keyword_set(objs) ? [objs,obj1] : obj1
         endif
      endfor
   endfor
   objs = objs[ sdss_selectobj(objs, /trim) ]
   objs = objs[uniq(objs.thing_id,sort(objs.thing_id))]
   objs = objs[where(djs_diff_angle(objs.ra, objs.dec, racen, deccen) $
    LT 1.49)]

   fibermag = 22.5 - 2.5*alog10(objs.fiberflux>0.01)

   ;----------
   ; Select the sky objects

   isky = where(objs.objc_type EQ 8, nsky)
   skyobj = design_struct(nsky)
   skyobj.ra = objs[isky].ra
   skyobj.dec = objs[isky].dec
   skyobj.mag = 30.
   skyobj.holetype = 'OBJECT'
   skyobj.objtype = 'SKY'
   skyobj.primtarget = 0
   skyobj.sectarget = sdss_flagval('TTARGET','SKY')
   skyobj.priority = 1

   ;----------
   ; Select the guide stars, using the criteria outlines in sdss-spectro/955

   psfmag = 22.5 - 2.5*alog10(objs.psfflux>0.01) ; no extinction-correction
   psf_ugcolor = psfmag[0,*] - psfmag[1,*]
   psf_grcolor = psfmag[1,*] - psfmag[2,*]
   psf_ricolor = psfmag[2,*] - psfmag[3,*]
   psf_izcolor = psfmag[3,*] - psfmag[4,*]

   iguide = where(psfmag[1,*] GT 13.0 AND psfmag[1,*] LT 15.5 $
    AND psf_grcolor GT 0.3 AND psf_grcolor LT 1.4 $
    AND psf_ricolor GT 0.0 AND psf_ricolor LT 0.7 $
    AND psf_izcolor GT -0.4 AND psf_izcolor LT 1.0 $
    AND objs.objc_type EQ 6, nguide)
   guideobj = design_struct(nguide)
   guideobj.ra = objs[iguide].ra
   guideobj.dec = objs[iguide].dec
   guideobj.mag = fibermag[*,iguide]
   guideobj.holetype = 'GUIDE'
   guideobj.objtype = 'NA'
   guideobj.primtarget = 0
   guideobj.sectarget = sdss_flagval('TTARGET','GUIDE_STAR')
   guideobj.priority = ((100. * psf_grcolor[iguide] / 1.4) > 1) < 100

   ;----------
   ; Select the SPECTROPHOTO_STD targets,
   ; using the criteria outlines in sdss-spectro/955,
   ; but make the mag limit a little bit fainter.

   ifstar = where(psfmag[1,*] GT 16.0 AND psfmag[1,*] LT 18.5 $
    AND psf_ugcolor GT 0.6 AND psf_ugcolor LT 1.2 $
    AND psf_grcolor GT 0.0 AND psf_grcolor LT 0.6 $
    AND psf_grcolor GT 0.75*psf_ugcolor-0.45 $
    AND objs.objc_type EQ 6, nfstar)
   cdist = sqrt( (psf_ugcolor[ifstar] - 0.934)^2 $
    + (psf_grcolor[ifstar] - 0.280)^2 $
    + (psf_ricolor[ifstar] - 0.101)^2 $
    + (psf_izcolor[ifstar] - 0.013)^2 )
   fstarobj = design_struct(nfstar)
   fstarobj.ra = objs[ifstar].ra
   fstarobj.dec = objs[ifstar].dec
   fstarobj.mag = fibermag[*,ifstar]
   fstarobj.holetype = 'OBJECT'
   fstarobj.objtype = 'SPECTROPHOTO_STD'
   iredden = where(psfmag[1,ifstar] GT 17.25, nredden)
   if (nredden GT 0) then fstarobj[iredden].objtype = 'REDDEN_STD'
   fstarobj.primtarget = 0
   fstarobj.sectarget = sdss_flagval('TTARGET',fstarobj.objtype)
   fstarobj.priority = $
    (100. * (cdist - min(cdist)) / (max(cdist) - min(cdist))) > 1

   ;----------
   ; Select the QSO targets

   iqso = where(psfmag[1,*] GT 17 AND psfmag[1,*] LT 22 $
    AND ((psf_ugcolor GE 0.40 AND psf_ugcolor LE 0.70 AND psf_grcolor LT 0.5) $
     OR (psf_ugcolor GE 0.70 AND psf_grcolor LT 0.45*psf_ugcolor-0.25) $
     OR (psf_grcolor LT 0.20 AND psf_ugcolor GT 0.70)) $
    AND psf_grcolor GT -0.8 AND psf_grcolor LT 0.7 $
    AND objs.objc_type EQ 6, nqso)
   qsoobj = design_struct(nqso)
   qsoobj.ra = objs[iqso].ra
   qsoobj.dec = objs[iqso].dec
   qsoobj.mag = fibermag[*,iqso]
   qsoobj.holetype = 'OBJECT'
   qsoobj.objtype = 'QSO'
   qsoobj.primtarget = sdss_flagval('TARGET','QSO_HIZ')
   qsoobj.sectarget = 0
   qsoobj.priority = 100.
; Should prioritize by distance from the stellar locus ???

   ;----------
   ; Select the LRG targets

   modelmag = 22.5 - 2.5*alog10(objs.modelflux>0.01) - objs.extinction
   ugcolor = modelmag[0,*] - modelmag[1,*]
   grcolor = modelmag[1,*] - modelmag[2,*]
   ricolor = modelmag[2,*] - modelmag[3,*]
; ...???

   ;----------
   ; Concatenate all of the objects

   allobj = [skyobj, guideobj, fstarobj, qsoobj]

   ;----------
   ; Read the existing plates, and lower priorities for existing spectra ???

   platelist, plist=plist
   adist = djs_diff_angle(plist.ra, plist.dec, racen, deccen)
   ikeep = where(adist LT 3. $
    AND strmatch(plist.status1d,'Done*') $
    AND (strmatch(plist.platequality,'good*') $
     OR strmatch(plist.platequality,'marginal*')))
   plist = plist[ikeep]
   readspec, plist.plate, mjd=plist.mjd, zans=zans, tsobj=tsobj
   indx = where(djs_diff_angle(zans.plug_ra, zans.plug_dec, racen, deccen) $
    LT 1.49 AND zans.zwarning EQ 0)
   zans = zans[indx]
   tsobj = tsobj[indx]

   spherematch, allobj.ra, allobj.dec, zans.plug_ra, zans.plug_dec, 1./3600, $
    i1, i2, d12

   ;----------
   ; Create the plug-map file

stop
; Need to split up the targets into multiple plates ???!!!
   allobj = [skyobj, guideobj, fstarobj, qsoobj]
   design_plate, allobj, racen=racen, deccen=deccen, tilenum=tilenum, $
    platenum=platenum, nstd=16, nminsky=64

end
;------------------------------------------------------------------------------
