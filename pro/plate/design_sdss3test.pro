;+
; NAME:
;   design_sdss3test
;
; PURPOSE:
;   Design test plates for SDSS-3
;
; CALLING SEQUENCE:
;   design_sdss3test, platenum, [ nminsky=, nstd= ]
;
; INPUTS:
;   platenum   - Plate number
;
; OPTIONAL INPUTS:
;   nminsky    - Minimum number of sky fibers per plate; default to 64
;   nstd       - Number of F star standards per plate; default to 16
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   Design our test plates for the Oct 2006 drill run at Princeton:
;     setenv,'PHOTO_REDUX=/u/dss/redux'
;     setenv,'PHOTO_RESOLVE=/u/dss/redux/resolve/full_02apr06'
;     setenv,'PHOTO_CALIB=/u/dss/redux/resolve/full_02apr06/calib/default0'
;     setenv,'BOSS_SPECTRO_REDUX=/u/dss/spectro_v5'
;     design_sdss3test, 2634
;     design_sdss3test, 2638
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
;   29-Oct-2006  Written by D. Schlegel and N. Padmanabhan, LBL
;-
;------------------------------------------------------------------------------
function design_struct, num

   result = create_struct( $
    name = 'TARGETDATA', $
    'objid'     ,  lonarr(5), $
    'ra'        ,  0.d, $
    'dec'       ,  0.d, $
    'mag'       , fltarr(5), $
    'holetype'  ,   '', $
    'objtype'   ,   '', $
    'priority'  ,   1L, $
    'primtarget',   0L, $
    'sectarget' ,   0L)
   if (keyword_set(num)) then result = replicate(result, num)

   return, result
end
;------------------------------------------------------------------------------
function design_decollide, ra, dec, mindist

   if (NOT keyword_set(mindist)) then mindist = 55./3600
   ingroup = spheregroup(ra, dec, mindist, firstgroup=firstgroup)

   return, firstgroup
end
;------------------------------------------------------------------------------
function lrg_dperp, modelmag

   grcolor = modelmag[1,*]-modelmag[2,*]
   ricolor = modelmag[2,*]-modelmag[3,*]
   dperp = (ricolor) - (grcolor)/8.d0
   return, dperp

end


function lrg_cpllel, modelmag

   grcolor = modelmag[1,*]-modelmag[2,*]
   ricolor = modelmag[2,*]-modelmag[3,*]
   c= 0.7
   cpllel = c*(grcolor) + (1.d0-c)*4.0d0*((ricolor) - 0.18)
   return, cpllel

end

function lrg_cperp, modelmag

   grcolor = modelmag[1,*]-modelmag[2,*]   
   ricolor = modelmag[2,*]-modelmag[3,*]
   cperp = (ricolor) - (grcolor/4.d0)-0.18
   return, cperp
end


function lrg_select_target, calibobj, lowz=lowz, hiz=hiz, all=all,$
                            maglim=maglim, dperp0 = dperp0

     ; Get the number of galaxies
     ngal = n_elements(calibobj)

     ; Set some reasonable defaults
     if (n_elements(maglim) NE 2) then maglim=[13.6d0, 19.7d0]
     if (NOT keyword_set(dperp0)) then dperp0=0.50
     if (NOT keyword_set(dperp1)) then dperp1=0.65

     ; Check to see if calibobj has any elements
     if (ngal EQ 0) then $
       message, 'ERROR : calibobj has no elements'

     ; Extract the relevant magnitudes
     ; We really should have used cmodelmags, and will later.
     modelmag = 22.5 - 2.5*alog10(calibobj.modelflux > 0.001)
     devmag = 22.5 - 2.5*alog10(calibobj.devflux > 0.001)
     expmag = 22.5 - 2.5*alog10(calibobj.expflux > 0.001)
     cmodelmag = devmag*calibobj.fracpsf + expmag*(1.0-calibobj.fracpsf)
     fibermag = 22.5 - 2.5*alog10(calibobj.fiberflux > 0.001)
     psfmag = 22.5 - 2.5*alog10(calibobj.psfflux > 0.001)
     ; Extinction correct
     modelmag = modelmag - calibobj.extinction
     psfmag = psfmag - calibobj.extinction
     cmodelmag = cmodelmag - calibobj.extinction
     fibermag = fibermag - calibobj.extinction
     
     ; Compute cperp and cpllel
     cperp = lrg_cperp(modelmag)
     dperp = lrg_dperp(modelmag)
     cpllel = lrg_cpllel(modelmag)     

     ; First lowz cuts
     if (keyword_set(lowz) OR keyword_set(all)) then begin
         icut1 = cmodelmag[2,*] LT (maglim[0] + cpllel/0.3d0)
         icut1 = icut1 AND (abs(cperp) LT 0.2d0)
         icut1 = icut1 AND (cmodelmag[2,*] LT maglim[1]) 
         icut1 = icut1 AND (cmodelmag[3,*] LT fibermag[3,*])
         ilrg = icut1 * 2L^1
     endif

     ;	Now hiz cuts
     maglow = 18.3+2*dperp1
     if (keyword_set(hiz) OR keyword_set(all)) then begin
         icut2 = (cmodelmag[3,*] GT 18.5d0) AND (cmodelmag[3,*] LT maglow)
         icut2 = icut2 AND ((modelmag[1,*] - modelmag[2,*]) GT 1.4d0)
         icut2 = icut2 AND ((modelmag[1,*] - modelmag[2,*]) LT 3.0d0)
         icut2 = icut2 AND ((modelmag[2,*] - modelmag[3,*]) LT 2.0d0)
         icut2 = icut2 AND (dperp GT dperp0) AND (dperp LT dperp1)
         icut2 = icut2 AND (cmodelmag[3,*] LT fibermag[3,*])
         ilrg  = ilrg + icut2*2L^2
     endif


     if (keyword_set(hiz) OR keyword_set(all)) then begin
         icut3 = (cmodelmag[3,*] GT 18.5) AND (cmodelmag[3,*] LT 20.0d0)
         icut3 = icut3 AND ((modelmag[1,*] - modelmag[2,*]) GT 1.4d0)
         icut3 = icut3 AND ((modelmag[1,*] - modelmag[2,*]) LT 3.0d0)
         icut3 = icut3 AND ((modelmag[2,*] - modelmag[3,*]) LT 2.0d0)
         icut3 = icut3 AND (dperp GT dperp1)
         icut3 = icut3 AND (cmodelmag[3,*] LT fibermag[3,*])
         ilrg  = ilrg + icut3*2L^3
     endif

     ; Only work with those objects that photo calls a galaxy, and don't
     ; have processing flags thrown.
     icut3 = reform(icut2*0b)
     indx = sdss_selectobj(calibobj, objtype=3, /trim) 
     if (indx[0] GT -1) then $      
       icut3[indx] = 1
     
     return, icut3*reform(ilrg) 
end

;------------------------------------------------------------------------------
; Assign priorities evenly spaced between 1 and 100 (larger is better),
; where the priorities are assigned randomly, but giving preference to
; those that don't have close neighbors in color space.
function design_color_prioritize, colors

   ndim = size(colors,/n_dimen)
   dims = size(colors,/dimens)
   if (ndim EQ 1) then ncolor = 1 $
    else ncolor = dims[1]
   num = dims[0]
   dd = fltarr(num)
   isort = sort(randomu(1234,num)) ; First sort these randomly
   dd[isort[0]] = 0.
   for i=1L, num-1L do begin
      ; Compute distance in color space to all previously assigned objects
      thisdist = fltarr(i)
      for j=0, ncolor-1 do $
       thisdist += (colors[isort[i]] - colors[isort[0:i-1]])^2
      dd[isort[i]] += min(thisdist)
   endfor

   psort = sort(dd)
   priority = lonarr(num)
   priority[psort] = (1. + 99. * (1+findgen(num))/float(num))

   return, priority
end
;------------------------------------------------------------------------------
pro design_sdss3test, platenum, nminsky=nminsky, nstd=nstd

   if (NOT keyword_set(nminsky)) then nminsky = 64
   if (NOT keyword_set(nstd)) then nstd = 16

   case platenum of
   2634: begin ; Centered at plate 406
      ; Note there is a bright star near ra=35.49, dec=0.40
      tilenum = 9549
      racen = 35.88296
      deccen = 0.1250122
;      runnum = [4263, 4874]
;      rerun = [137, 137]
      runnum = [2709, 3384]
      rerun = [137, 137]
      end
   2638: begin ; Centered at plate 416
      tilenum = 9553
      racen =  55.49162
      deccen = 0.01375204
;      runnum = [4136, 4145, 4874]
;      rerun = [137, 137, 137]
      runnum = [2709, 3384]
      rerun = [137, 137]
      end
   else: begin
      print, 'Unknown plate number!'
      return
      end
   endcase

   ;----------
   ; Read all objects from fpObjc files

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
             rerun=rerun[irun], /silent)
            if (keyword_set(obj1)) then $
             objs = keyword_set(objs) ? [objs,obj1] : obj1
         endif
      endfor
   endfor

   ;----------
   ; Add the OBJID array for the plugmaps

   objid = replicate(create_struct('objid', lonarr(5)), n_elements(objs))
   objid.objid[0] = objs.run
   objid.objid[1] = long(objs.rerun)
   objid.objid[2] = objs.camcol
   objid.objid[3] = objs.field
   objid.objid[4] = objs.id
   objs = struct_addtags(objs, objid)

   ;----------
   ; Trim based upon the default options in the SDSS_SELECTOBJ routine,
   ; as well as the INTERP_CENTER.
   objs = objs[ sdss_selectobj(objs, /trim) ]
   indx = where( $
    (objs.objc_flags2 AND sdss_flagval('OBJECT2','INTERP_CENTER')) EQ 0)
   objs = objs[indx]
   objs = objs[uniq(objs.thing_id,sort(objs.thing_id))]
   objs = objs[where(djs_diff_angle(objs.ra, objs.dec, racen, deccen) $
    LT 1.49)]

   ; Logic to see which objects are or maybe are moving...
   qmoved = $
    (objs.objc_flags2 AND sdss_flagval('OBJECT2','DEBLENDED_AS_MOVING')) EQ 0 $
    AND (abs(objs.rowv) GE 3.0*abs(objs.rowverr) $
     OR abs(objs.colv) GE 3.0*abs(objs.colverr))
   qmaybe_moved = $
    (objs.objc_flags2 AND sdss_flagval('OBJECT2','DEBLENDED_AS_MOVING')) EQ 0 $
    AND (abs(objs.rowv) GE 3.0*abs(objs.rowverr) $
     OR abs(objs.colv) GE 3.0*abs(objs.colverr) $
     OR objs.rowverr LE 0 OR objs.colverr LE 0)

   fibermag = 22.5 - 2.5*alog10(objs.fiberflux>0.01)

   ;----------
   ; Select the sky objects

   isky = where(objs.objc_type EQ 8, nsky)
   skyobj = design_struct(nsky)
   skyobj.objid = objs[isky].objid
   skyobj.ra = objs[isky].ra
   skyobj.dec = objs[isky].dec
   skyobj.mag = 30.
   skyobj.holetype = 'OBJECT'
   skyobj.objtype = 'SKY'
   skyobj.primtarget = 0
   skyobj.sectarget = sdss_flagval('TTARGET','SKY')
   skyobj.priority = 1

   ;----------
   ; Select the guide stars, using the criteria outlines in sdss-spectro/955,
   ; but changing the magnitude limit to go 0.5 mag fainter.
   ; This is because stars can saturate at g=15.5 in good-seeing runs.

   psfmag = 22.5 - 2.5*alog10(objs.psfflux>0.01) ; no extinction-correction
   psf_ugcolor = reform(psfmag[0,*] - psfmag[1,*])
   psf_grcolor = reform(psfmag[1,*] - psfmag[2,*])
   psf_ricolor = reform(psfmag[2,*] - psfmag[3,*])
   psf_izcolor = reform(psfmag[3,*] - psfmag[4,*])

   iguide = where(psfmag[1,*] GT 13.5 AND psfmag[1,*] LT 16.0 $
    AND psf_grcolor GT 0.3 AND psf_grcolor LT 1.4 $
    AND psf_ricolor GT 0.0 AND psf_ricolor LT 0.7 $
    AND psf_izcolor GT -0.4 AND psf_izcolor LT 1.0 $
    AND objs.objc_type EQ 6, nguide)
   iguide = iguide[ sdss_selectobj(objs[iguide], ancestry='single', $
    /trim, count=nguide) ]
   guideobj = design_struct(nguide)
   guideobj.objid = objs[iguide].objid
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

   psfmag = 22.5 - 2.5*alog10(objs.psfflux>0.01) - objs.extinction ; correct!
   psf_ugcolor = reform(psfmag[0,*] - psfmag[1,*])
   psf_grcolor = reform(psfmag[1,*] - psfmag[2,*])
   psf_ricolor = reform(psfmag[2,*] - psfmag[3,*])
   psf_izcolor = reform(psfmag[3,*] - psfmag[4,*])

   ifstar = where(psfmag[1,*] GT 16.0 AND psfmag[1,*] LT 18.5 $
    AND psf_ugcolor GT 0.6 AND psf_ugcolor LT 1.2 $
    AND psf_grcolor GT 0.0 AND psf_grcolor LT 0.6 $
    AND psf_grcolor GT 0.75*psf_ugcolor-0.45 $
    AND objs.objc_type EQ 6, nfstar)
   ifstar = ifstar[ sdss_selectobj(objs[ifstar], ancestry='single', $
    /trim, count=nfstar) ]
   cdist = sqrt( (psf_ugcolor[ifstar] - 0.934)^2 $
    + (psf_grcolor[ifstar] - 0.280)^2 $
    + (psf_ricolor[ifstar] - 0.101)^2 $
    + (psf_izcolor[ifstar] - 0.013)^2 )
   fstarobj = design_struct(nfstar)
   fstarobj.objid = objs[ifstar].objid
   fstarobj.ra = objs[ifstar].ra
   fstarobj.dec = objs[ifstar].dec
   fstarobj.mag = fibermag[*,ifstar]
   fstarobj.holetype = 'OBJECT'
   fstarobj.objtype = 'SPECTROPHOTO_STD'
; Including the lines below would call some objects REDDEN_STD,
; but that may not be supported by the DESIGN_PLATE code.
   iredden = where(psfmag[1,ifstar] GT 17.25, nredden)
   if (nredden GT 0) then fstarobj[iredden].objtype = 'REDDEN_STD'
   fstarobj.primtarget = 0
   fstarobj.sectarget = sdss_flagval('TTARGET',fstarobj.objtype)
   fstarobj.priority = $
    (100. * (cdist - min(cdist)) / ((max(cdist) - min(cdist))>1)) > 1

   ;----------
   ; Select the QSO targets

   iqso = where(psfmag[1,*] GT 17 AND psfmag[1,*] LT 22 $
    AND psf_ugcolor GT 0.40 AND psf_grcolor LT 0.20 $
    AND objs.fiberflux[1] LT objs.psfflux[1] $
    AND (qmaybe_moved EQ 0) $
;    AND ((psf_ugcolor GE 0.40 AND psf_ugcolor LE 0.75 AND psf_grcolor LT 0.5) $
;     OR (psf_ugcolor GE 0.75 AND psf_grcolor LT 0.45*psf_ugcolor-0.25) $
;     OR (psf_grcolor LT 0.20 AND psf_ugcolor GT 0.75)) $
;    AND psf_grcolor GT -0.8 AND psf_grcolor LT 0.7 $
    AND objs.objc_type EQ 6, nqso)
   qsoobj = design_struct(nqso)
   qsoobj.objid = objs[iqso].objid
   qsoobj.ra = objs[iqso].ra
   qsoobj.dec = objs[iqso].dec
   qsoobj.mag = fibermag[*,iqso]
   qsoobj.holetype = 'OBJECT'
   qsoobj.objtype = 'QSO'
   qsoobj.primtarget = sdss_flagval('TARGET','QSO_HIZ')
   qsoobj.sectarget = 0
   qsoobj.priority = design_color_prioritize([ [psf_ugcolor[iqso]/3.], $
    [psf_grcolor[iqso]], [reform(psfmag[1,iqso]/50.)] ])

; Compare to the co-add catalogs
;varcat=mrdfits('/u/schlegel/varcat/varcat-ra30.fits',1)
;spherematch, qsoobj.ra, qsoobj.dec, varcat.ra, varcat.dec, 1./3600, $
; i1, i2, d12
;varmag = 22.5 - 2.5*alog10(varcat.modelflux_clip_mean>0.01) - varcat.extinction
;var_ugcolor = reform(varmag[0,*] - varmag[1,*])
;var_grcolor = reform(varmag[1,*] - varmag[2,*])
;splot, psf_ugcolor[iqso[i1]], psf_grcolor[iqso[i1]], ps=3
;soplot, var_ugcolor[i2], var_grcolor[i2], ps=3, color='red'

   ;----------
   ; Select the LRG targets, and further trim to i(fiber) > 16

   ilist = lrg_select_target(objs, /all)
   ilrg = where(ilist GT 0 AND fibermag[3,*] GT 16, nlrg)
   lrgobj = design_struct(nlrg)
   lrgobj.objid = objs[ilrg].objid
   lrgobj.ra = objs[ilrg].ra
   lrgobj.dec = objs[ilrg].dec
   lrgobj.mag = fibermag[*,ilrg]
   lrgobj.holetype = 'OBJECT'
   lrgobj.objtype = 'GALAXY'
   lrgobj.primtarget = sdss_flagval('TARGET','GALAXY_RED')
   lrgobj.sectarget = 0
   lrgobj.priority = ilist[ilrg]

   ;----------
   ; Read the existing plates, and lower priorities for existing spectra

   platelist, plist=plist
   adist = djs_diff_angle(plist.ra, plist.dec, racen, deccen)
   ikeep = where(adist LT 3. $
    AND strmatch(plist.status1d,'Done*') $
    AND (strmatch(plist.platequality,'good*') $
     OR strmatch(plist.platequality,'marginal*')))
   plist = plist[ikeep]
   readspec, plist.plate, mjd=plist.mjd, zans=zans, tsobj=tsobj, /silent
   indx = where(djs_diff_angle(zans.plug_ra, zans.plug_dec, racen, deccen) $
    LT 1.49 AND zans.zwarning EQ 0)
   zans = zans[indx]
   tsobj = tsobj[indx]

   sciobj = [qsoobj, lrgobj]
   nsci = n_elements(sciobj)
   spherematch, sciobj.ra, sciobj.dec, zans.plug_ra, zans.plug_dec, 1./3600, $
    i1, i2, d12
   qexist = bytarr(nsci)
   if (i1[0] NE -1) then qexist[i1] = 1B

   ;----------
   ; Assign science targets to each plate

   ; Decide how many plates we need based upon the number of targets
   ; not yet observed...

   maxtarget = 640 - nminsky - nstd ; Max science targets per plate
   nplate = ceil((nsci - total(qexist)) / float(maxtarget))
   splog, 'Number of required plates = ', nplate
   qassign = lonarr(nplate,nsci) ; Set to 1 wherever a target is assigned

   ; Assign science targets one at a time, starting with those not yet observed
   ; We loop over each object in a random order, then over each plate randomly
   iseed = 123456
   irandom = sort(randomu(iseed, nsci) + qexist)
   for i=0L, nsci-1L do begin
      porder = sort(randomu(iseed+i, nplate))
      qadd = 0B
      for j=0L, nplate-1L do begin
         ; Try putting object number irandom[i] on plate porder[j]
         if (total(qassign[porder[j],*]) LT maxtarget AND qadd EQ 0) then begin
            indx = where(qassign[porder[j],*], ct)
            if (ct EQ 0) then qadd = 1B $
             else qadd = min(djs_diff_angle(sciobj[irandom[i]].ra, $
              sciobj[irandom[i]].dec, sciobj[indx].ra, sciobj[indx].dec)) $
               GT 55./3600
            qassign[porder[j],irandom[i]] = qadd
         endif
      endfor
   endfor

   splog, 'Number of SKY fibers = ', nsky
   splog, 'Number of GUIDE stars = ', nguide
   splog, 'Number of F stars = ', nfstar
   splog, 'Number of QSO targets = ', nqso
   splog, 'Number of LRG targets = ', nlrg

   splog, 'Number of objects targetted = ', $
    long(total(total(qassign,1) NE 0))
   splog, 'Number of objects not targetted = ', $
    long(total(total(qassign,1) EQ 0))
   splog, 'Number of unobserved objects not targetted = ', $
    long(total(qexist EQ 0 AND total(qassign,1) EQ 0))
   splog, 'Number of observed objects not targetted = ', $
    long(total(qexist EQ 1 AND total(qassign,1) EQ 0))
   splog, 'Number of observed objects targetted = ', $
    long(total(qexist EQ 1 AND total(qassign,1) EQ 1))

   ;----------
   ; Create the plug-map files

   for iplate=0L, nplate-1L do begin
      splog, 'Generating plate number = ', platenum+iplate
      allobj = [skyobj, guideobj, fstarobj, sciobj[where(qassign[iplate,*])]]
      design_plate, allobj, racen=racen, deccen=deccen, $
       tilenum=tilenum+iplate, platenum=platenum+iplate, $
       nstd=nstd, nminsky=nminsky
;      simple_plate, allobj, racen=racen, deccen=deccen, $
;       tilenum=tilenum+iplate, platenum=platenum+iplate
   
   endfor

   return
end
;------------------------------------------------------------------------------
