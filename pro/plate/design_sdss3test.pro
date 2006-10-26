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
    name = 'TARGETDATA', $
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
     psfmag = 22.5 - 2.5*alog10(calibobj.psfflux > 0.001)
     ; Extinction correct
     modelmag = modelmag - calibobj.extinction
     psfmag = psfmag - calibobj.extinction
     
     ; Compute cperp and cpllel
     cperp = lrg_cperp(modelmag)
     dperp = lrg_dperp(modelmag)
     cpllel = lrg_cpllel(modelmag)

     ; First lowz cuts
     if (keyword_set(lowz) OR keyword_set(all)) then begin
         icut1 = modelmag[2,*] LT (maglim[0] + cpllel/0.3d0)
         icut1 = icut1 AND (abs(cperp) LT 0.2d0)
         icut1 = icut1 AND (modelmag[2,*] LT maglim[1]) 
         ilrg = icut1 * 2L^1
     endif

     ;	Now hiz cuts
     maglow = 18.3+2*dperp1
     if (keyword_set(hiz) OR keyword_set(all)) then begin
         icut2 = (modelmag[3,*] GT 18.5d0) AND (modelmag[3,*] LT maglow)
         icut2 = icut2 AND ((modelmag[1,*] - modelmag[2,*]) GT 1.4d0)
         icut2 = icut2 AND ((modelmag[1,*] - modelmag[2,*]) LT 3.0d0)
         icut2 = icut2 AND ((modelmag[2,*] - modelmag[3,*]) LT 2.0d0)
         icut2 = icut2 AND (dperp GT dperp0) AND (dperp LT dperp1)
         ilrg  = ilrg + icut2*2L^2
     endif


     if (keyword_set(hiz) OR keyword_set(all)) then begin
         icut3 = (modelmag[3,*] GT 18.5) AND (modelmag[3,*] LT 20.0d0)
         icut3 = icut3 AND ((modelmag[1,*] - modelmag[2,*]) GT 1.4d0)
         icut3 = icut3 AND ((modelmag[1,*] - modelmag[2,*]) LT 3.0d0)
         icut3 = icut3 AND ((modelmag[2,*] - modelmag[3,*]) LT 2.0d0)
         icut3 = icut3 AND (dperp GT dperp1)
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
   psf_ugcolor = reform(psfmag[0,*] - psfmag[1,*])
   psf_grcolor = reform(psfmag[1,*] - psfmag[2,*])
   psf_ricolor = reform(psfmag[2,*] - psfmag[3,*])
   psf_izcolor = reform(psfmag[3,*] - psfmag[4,*])

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
    (100. * (cdist - min(cdist)) / ((max(cdist) - min(cdist))>1)) > 1

   ;----------
   ; Select the QSO targets

   iqso = where(psfmag[1,*] GT 17 AND psfmag[1,*] LT 22 $
    AND psf_ugcolor GT 0.40 AND psf_grcolor LT 0.20 $
;    AND ((psf_ugcolor GE 0.40 AND psf_ugcolor LE 0.75 AND psf_grcolor LT 0.5) $
;     OR (psf_ugcolor GE 0.75 AND psf_grcolor LT 0.45*psf_ugcolor-0.25) $
;     OR (psf_grcolor LT 0.20 AND psf_ugcolor GT 0.75)) $
;    AND psf_grcolor GT -0.8 AND psf_grcolor LT 0.7 $
    AND objs.objc_type EQ 6, nqso)
   qsoobj = design_struct(nqso)
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
   ; Select the LRG targets

   ilist = lrg_select_target(objs, /all)
   ilrg = where(ilist GT 0)
   lrgobj = design_struct(nlrg)
   lrgobj.ra = objs[ilrg].ra
   lrgobj.dec = objs[ilrg].dec
   lrgobj.mag = fibermag[*,ilrg]
   lrgobj.holetype = 'OBJECT'
   lrgobj.objtype = 'GALAXY'
   lrgobj.primtarget = sdss_flagval('TARGET','GALAXY_RED')
   lrgobj.sectarget = 0
   lrgobj.priority = ilist[ilrg]

;   modelmag = 22.5 - 2.5*alog10(objs.modelflux>0.01) - objs.extinction
;   ugcolor = reform(modelmag[0,*] - modelmag[1,*])
;   grcolor = reform(modelmag[1,*] - modelmag[2,*])
;   ricolor = reform(modelmag[2,*] - modelmag[3,*])

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
