; Code to design the Taurus run plates
;------------------------------------------------------------------------------
function design_taurus_add, alldata, objs, $
 holetype=holetype, objtype=objtype, $
 primtarget=primtarget, sectarget=sectarget, priority=priority

   ;----------
   ; Create blank output structure

   blankdat = create_struct( $
    'OBJID'     ,    lonarr(5), $
    'HOLETYPE'  , '', $
    'RA'        , 0.d0, $
    'DEC'       , 0.d0, $
    'MAG'       , fltarr(5), $
    'OBJTYPE'   , '', $
    'PRIORITY'  , 0L, $
    'PRIMTARGET', 0L, $
    'SECTARGET' , 0L)

   appdata = replicate(blankdat, n_elements(objs))
   appdata.objid[0] = objs.run
   appdata.objid[1] = objs.rerun
   appdata.objid[2] = objs.camcol
   appdata.objid[3] = objs.field
   appdata.objid[4] = objs.id
   appdata.holetype = holetype
   appdata.mag = objs.fibercounts
   appdata.objtype = objtype
   appdata.ra = objs.ra[2]
   appdata.dec = objs.dec[2]
   appdata.priority = priority
   appdata.primtarget = primtarget
   appdata.sectarget = sectarget

   if (NOT keyword_set(alldata)) then return, appdata $
    else return, [alldata, appdata]
end

;------------------------------------------------------------------------------
pro design_taurus, designnum

   if (NOT keyword_set(designnum)) then designnum = 1

   if (designnum EQ 1) then begin
      runvec = [3511, 3557]
      racen =  [65.15, 65.15, 65.15, 65.15]
;      deccen = [26.30, 28.65, 31.00, 33.35]
      deccen = [26.80, 28.95, 31.10, 33.25]
      tileid = 9361 + lindgen(4)
      plateid = 1250 + lindgen(4)
   endif else begin
      runvec = [3512, 3559]
      racen =  [64.00, 66.50, 69.00, 71.50, 74.00]
      deccen = [25.52, 25.53, 25.50, 25.41, 25.30]
      tileid = 9365 + lindgen(5)
      plateid = 1254 + lindgen(5)
   endelse
   nplate = n_elements(racen)
   rerun = 125

   ;----------
   ; The colors of BD+17 are from Hogg in sdss-calib/845.

   bd17color = [0.934, 0.280, 0.101, 0.013]

   ;----------
   ; Read in the fpObj files -- or read a previously saved file if it exists

   columns = ['ID', 'PARENT', 'OBJC_TYPE', 'OBJC_FLAGS', 'OBJC_FLAGS2', $
    'ROWC', 'COLC', 'OBJC_ROWC', 'OBJC_COLC', 'AIRMASS', $
    'PSFCOUNTS', 'PSFCOUNTSERR', $
    'FIBERCOUNTS', 'FIBERCOUNTSERR', $
    'COUNTS_MODEL', 'COUNTS_MODELERR', $
    'RA', 'DEC' ]

   for irun=0, n_elements(runvec)-1 do begin
      for camcol=1, 6 do begin
         cachefile = string(runvec[irun], camcol, $
          format='("taurus-cache-",i6.6,"-",i1.1,".fits")')
         objs1 = mrdfits(cachefile, 1)
         if (keyword_set(objs1)) then begin
            splog, 'Read previously-cached FITS file: ' + cachefile
         endif else begin
            splog, 'Generate cache FITS file: ' + cachefile
            fstart = sdss_fieldrange(runvec[irun], fend=fend)
            fieldnum = fstart + lindgen(fend-fstart+1)
            objs1 = rdss_obj(runvec[irun], rerun, camcol, $
             fstart+lindgen(fend-fstart+1), /calib, /coord, /fieldnum, $
             columns=columns)
            indx = objc_select(objs1, /trim) ; Preliminary trimming
            objs1 = objs1[indx]
            indx = where(objs1.fibercounts[1] LT 19.5 $
             OR objs1.fibercounts[2] LT 19.5 $
             OR objs1.fibercounts[3] LT 19.5 $
             OR objs1.objc_type EQ 8) ; Bright or a blank sky position
            mwrfits_chunks, objs1, cachefile, /create, chunksize=10000L
         endelse
         if (NOT keyword_set(objs)) then objs = objs1 $
          else objs = [objs, objs1]
      endfor
   endfor

   ;----------
   ; Compute the SFD extinction

   red_fac = [5.155, 3.793, 2.751, 2.086, 1.479]

   splog, 'Reading SFD reddening maps...'
   euler, objs.ra[2], objs.dec[2], ll, bb, 1
   sfdval = red_fac # dust_getval(ll, bb, /interp, /noloop)
   sfdval = sfdval * 0.75 ; Fudge factor!!!

;   splog, 'Reading blue-edge reddening maps...'
;   junk = red_fac # rdss_blue_edge(objs.run, objs.camcol, objs.field, $
;    rerun=objs.rerun) / (red_fac[2] - red_fac[3])

   ;----------
   ; Classify each object

   splog, 'Classifying objects...'

   ; These are all the flags we may be interested in...
   bflag = djs_int2bin(ulong(objs.objc_flags), ndigit=32)
   qblend = (bflag[3,*] EQ 1 AND bflag[6,*] EQ 0)[*]
   qbright = (bflag[1,*])[*]
   qchild = (bflag[4,*])[*]
   qsingle = ((qblend EQ 0) AND (qbright EQ 0) AND (qchild EQ 0))[*]

   ; Ask whether the detection is good in each of the following bands
   qgoodmag = objs.fibercountserr LT 0.10 AND objs.fibercountserr GT 0
   qnot2bright = objs.fibercounts[1] GT 15.0 $
    AND objs.fibercounts[2] GT 15.0 $
    AND objs.fibercounts[3] GT 15.0

   qstellar = objs.objc_type EQ 6
   qfuzzy = objs.objc_type EQ 3
   starcolor = objs.psfcounts[0:3] - objs.psfcounts[1:4]
   fuzzcolor = objs.counts_model[0:3] - objs.counts_model[1:4]
   sfdcolor = sfdval[0:3,*] - sfdval[1:4,*]
   bd17diff = transpose(sqrt( $
    (starcolor[1,*] - sfdcolor[1,*] - bd17color[1])^2 $
    + (starcolor[2,*] - sfdcolor[2,*] - bd17color[2])^2 ))
   ri_edgedist = transpose(starcolor[2,*] - sfdcolor[2,*] - 0.10)

   qguide = qstellar AND qsingle AND objs.psfcounts[2] LT 16.5 $
    AND (qgoodmag[1,*])[*] AND (qgoodmag[2,*])[*] AND (qgoodmag[3,*])[*] $
    AND qnot2bright AND ri_edgedist GT -0.25 AND ri_edgedist LT 0.50
   qsky = objs.objc_type EQ 8
   qsphoto = qstellar AND qsingle AND (total(qgoodmag,1) EQ 5) $
    AND objs.psfcounts[2] LT 17.5 $
    AND objs.psfcounts[3] LT 17.5 $
    AND qnot2bright AND bd17diff LT 0.20
   qgalaxy = qfuzzy AND qnot2bright $
    AND objs.counts_model[2] LT 19.0 $
    AND objs.counts_model[3] LT 19.0
   qbluestar = qstellar AND (ri_edgedist LT 0.10) $
    AND qnot2bright AND objs.psfcounts[2] LT 17.0 $
    AND objs.psfcounts[3] LT 17.0

   ;----------
   ; Loop over each plate center

   for iplate=0, nplate-1 do begin
      splog, 'Designing plate number ', plateid[iplate]

      adist = djs_diff_angle(objs.ra[2], objs.dec[2], $
       racen[iplate], deccen[iplate])
      iguide = where(qguide AND adist LT 1.485, nguide)
      isky = where(qsky AND adist LT 1.485, nsky)
      isphoto = where(qsphoto AND adist LT 1.485, nsphoto)
      igalaxy = where(qgalaxy AND adist LT 1.485, ngalaxy)
      ibluestar = where(qbluestar AND adist LT 1.485, nbluestar)

      splog, 'Possible NSKY = ', nsky
      splog, 'Possible NGUIDE = ', nguide
      splog, 'Possible NSPHOTO = ', nsphoto
      splog, 'Possible NBLUESTAR = ', nbluestar

      alldata = 0
      randompriority = long(100 * randomu(43265L, n_elements(objs))) > 1

      alldata = design_taurus_add(alldata, objs[iguide], $
       holetype='GUIDE', objtype='NA', $
       priority=randompriority[iguide], primtarget=0, sectarget=64L)

      alldata = design_taurus_add(alldata, objs[isky], $
       holetype='OBJECT', objtype='SKY', $
       priority=randompriority[isky], primtarget=0, sectarget=16L)

      alldata = design_taurus_add(alldata, objs[isphoto], $
       holetype='OBJECT', objtype='SPECTROPHOTO_STD', $
       priority=randompriority[isphoto], primtarget=0, sectarget=32L)

;splot,objs[igalaxy].dec[2],objs[igalaxy].fibercounts[2],ps=4,syms=2 

      ; Select up to the 100 brightest galaxies (only), and given them
      ; a very high priority.
      igalaxy = igalaxy[(sort(objs[igalaxy].fibercounts[2]))[0:(ngalaxy<99)]]
      alldata = design_taurus_add(alldata, objs[igalaxy], $
       holetype='OBJECT', objtype='GALAXY', $
       priority=1000L, primtarget=2L^6, sectarget=0)

;i = where(qstellar AND objs.psfcounts[2] LT 19)
;splot, ri_edgedist[i], objs[i].psfcounts[3], ps=3,xr=[-3,3]
;soplot, ri_edgedist[ibluestar], objs[ibluestar].psfcounts[3], ps=3,color='red'
      ; Prioritize the blue stars by how blue they are
      thispriority = (((2 - ri_edgedist[ibluestar]) * 100) < 999) > 1
      alldata = design_taurus_add(alldata, objs[ibluestar], $
       holetype='OBJECT', objtype='SERENDIPITY_MANUAL', $
       priority=thispriority, primtarget=2L^24, sectarget=0)

      ; Select up to 150 sky fibers if possible
      nminsky = long(0.70*nsky) < 150
      design_plate, alldata, racen=racen[iplate], deccen=deccen[iplate], $
       tilenum=tileid[iplate], platenum=plateid[iplate], $
       nstd=32, nminsky=nminsky, nextra=50
      print, 'Press a key to continue...'
      junk = get_kbrd(1)
stop

   endfor

   return
end
;------------------------------------------------------------------------------
