;+
; NAME:
;   readplugmap
;
; PURPOSE:
;   Read plugmap file and append tags as requested
;
; CALLING SEQUENCE:
;   plugmap = readplugmap( plugfile, [ spectrographid, plugdir=, $
;    /apotags, /deredden, /calibobj, exptime=, hdr=, fibermask=, _EXTRA= ] )
;
; INPUTS:
;   plugfile  - Name of Yanny-parameter plugmap file
;
; OPTIONAL INPUTS:
;   spectrographid  - The spectrograph number, either 1 or 2;
;                     if not set (or 0), then return all object fibers
;   plugdir   - Directory for PLUGFILE
;   apotags   - If set, then add a number of tags to the output structure
;               constructed from the Yanny header.  These tags are:
;               CARTID, PLATEID, TILEID, RAPLATE, DECPLATE, REDDEN_MED.
;               Also add the tags FIBERSN[3], SYTHMAG[3] which are used
;               by the on-the-mountain reductions.
;   deredden  - If set, then deredden the MAG fields using the median
;               reddening value for the entire plate as found in the
;               Yanny header of the plugmap file; this is done for the
;               on-the-mountain reductions.
;   exptime   - Default exposure time SCI_EXPTIME to add to the output;
;               if there are multiple pointings and EXPTIME is set, then
;               the exposure time in each pointing is scaled such that
;               their sum is EXPTIME.
;   calibobj  - If set, then add a CALIBFLUX,CALIBFLUX_IVAR entries based upon
;               the calibObj (deprecated) or photoPlate files.
;               For stellar objects, this contains the
;               PSF fluxes in nMgy.  For galaxies, it contains the fiber fluxes
;               multiplied by the median (PSF/fiber) flux ratio for stars.
;               The MAG fields are left unchanged.
;               For objects with no calibObj entry, simply set these fields as:
;                 CALIBFLUX = 22.5 - 2.5*alog10(MAG), CALIBFLUX_IVAR = 0.
;               We apply the putative AB corrections to these fluxes
;               (but not to the MAG values).
;               Also add the SFD reddening values as SFD_EBV.
;   _EXTRA    - Keywords for PLUG2TSOBJ(), such as MJD,INDIR
;
; OUTPUTS:
;   plugmap   - Plugmap structure
;
; OPTIONAL OUTPUTS:
;   hdr       - Header from Yanny-formatted plugmap file
;   fibermask - Byte array with bits set for unknown fibers
;
; COMMENTS:
;   Do not use the calibObj structure if more than 10% of the non-sky
;   objects do not have fluxes.
;           
;   Reads $IDLSPEC2D_DIR/opfiles/washers.par for ZOFFSET status overrides.
;   The original plugmap files reflect what we wanted to do; the overrides
;   and the return of this function reflect what we actually did.
;
; EXAMPLES:
;
; BUGS:
;   The AB corrections are hard-wired to be the same as in the photoop
;   product as of 18 Feb 2004.
;
; PROCEDURES CALLED:
;   djs_filepath()
;   dust_getval()
;   euler
;   plug2tsobj()
;   sdss_flagval()
;   splog
;   struct_addtags()
;   yanny_par()
;   yanny_read
;
; REVISION HISTORY:
;   29-Jan-2001  Written by S. Burles, FNAL
;   07-Aug-2012  Added ZOFFSET overrides; S. Bailey, LBL
;   21-Jun-2021  Editted by S Morrison to add correction check for catalog
;-
;------------------------------------------------------------------------------
function readplugmap_sort, plugmap, hdr, fibermask=fibermask, plates=plates

   qobj = strmatch(plugmap.holetype,'OBJECT')
   indx = where(qobj, nfiber)
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber) $
   else if (n_elements(fibermask) NE nfiber) then $
    message, 'Number of elements in FIBERMASK do not match NFIBER'
    
   badstdmask = bytarr(nfiber)

    
   blankmap = plugmap[0]
   struct_assign, {junk:0}, blankmap
   plugsort = replicate(blankmap, nfiber)
   plugsort.holetype = 'OBJECT'
   plugsort.objtype = 'NA'
   plugsort.fiberid = -1
   if keyword_set(plates) then begin
      programname = (yanny_par(hdr, 'programname'))[0] 
      program_tag1 = replicate( $
        {program: programname}, nfiber)
      program_tag2 = replicate( $
        {program: programname}, n_elements(plugmap))
      
      plugsort = struct_addtags(plugsort, $
       replicate(create_struct('ORG_FIBERID', 0L), n_elements(plugsort)))

      plugmap = struct_addtags(plugmap, $
        replicate(create_struct('ORG_FIBERID', 0L), n_elements(plugmap)))
      plugmap.ORG_FIBERID=plugmap.fiberid

      plugsort = struct_addtags(plugsort, program_tag1)
      plugmap = struct_addtags(plugmap, program_tag2)
      if strmatch(programname, '*MWM*', /fold_case) eq 1 then begin
        spht = strmatch(plugmap.objtype, 'SPECTROPHOTO_STD')
        for i=0, n_elements(plugmap)-1 do begin
          if spht[i] then begin
            if plugmap[i].mag[3] ge 18.0 then begin
               ;plugmap[i].fiberid=-1
               badstdmask[i] = 1
            endif
          endif
        endfor
        ;for i=0, n_elements(plugmap)-1 do begin
        ;  if plugmap[i].mag[3] le 14.5 and plugmap[i].mag[3] ge 10.0  then begin
        ;       plugmap[i].fiberid=-1
        ;  endif
        ;endfor
      endif
      ;igood = where(qobj AND plugmap.fiberid GT 0 AND plugmap.spectrographid EQ 1, ngood); check the number -1
      ;igoodapoge = where(qobj AND plugmap.fiberid GT 0 AND plugmap.spectrographid EQ 2, napogee)
      igood = where(qobj AND plugmap.fiberid GT 0 AND plugmap.spectrographid EQ 1 AND badstdmask EQ 0, ngood); check the number -1
      iplugged = where(qobj AND plugmap.fiberid GT 0 AND plugmap.spectrographid EQ 1, nplugged); check the number -1
      igoodapoge = where(qobj AND plugmap.fiberid GT 0 AND plugmap.spectrographid EQ 2 AND badstdmask EQ 0, napogee)
   endif else begin
      igood = where(qobj AND plugmap.fiberid GT 0, ngood)
      iplugged = igood
   endelse
   if (ngood EQ 0) then $
    message, 'No fibers found in plugmap!'
   ;iplace = plugmap[igood].fiberid - 1
   ;plugsort[iplace] = plugmap[igood]
   iplace = plugmap[iplugged].fiberid - 1
   plugsort[iplace] = plugmap[iplugged]
   ; Set the appropriate fibermask bit if a fiber not found in plugmap file.
   ; Do this by first setting all bits to 1, then unsetting the good ones.
   fibermask = fibermask OR fibermask_bits('NOPLUG')
   fibermask[iplace] = fibermask[iplace] - fibermask_bits('NOPLUG')

   ; Fill in unplugged fibers with arbitrary entries, and assign
   ; them a FIBERID.  After this, plugsort.fiberid should run from 1...nfiber
   imissing = where(plugsort.fiberid LE 0, nmissing)
   splog, 'Number of missing fibers: ', nmissing-napogee
   splog, 'Number of Invalid Standards: ',total(badstdmask,/INTEGER)

   if (nmissing GT 0) then begin
      if keyword_set(plates) then begin
         ifill = where(qobj AND plugmap.fiberid LE 0 OR plugmap.spectrographid EQ 2, nfill)
      endif else begin
         ifill = where(qobj AND plugmap.fiberid LE 0, nfill)
      endelse
      plugsort[imissing] = plugmap[ifill]
      plugsort[imissing].fiberid = imissing + 1
   endif
   ;print,nmissing,nfill,ngood,nfiber
   ;print,plugsort.fiberid
   return, plugsort
end
;------------------------------------------------------------------------------
function readplugmap, plugfile, spectrographid, plugdir=plugdir, $
 apotags=apotags, deredden=deredden, exptime=exptime, calibobj=calibobj, $
 hdr=hdr, fibermask=fibermask, plates=plates, gaiaext=gaiaext, $
 MWM_fluxer=MWM_fluxer, _EXTRA=KeywordsForPhoto

   hdr = 0 ; Default return value
   if (keyword_set(fibermask)) then $
    message, 'FIBERMASK is already set!'

   ; The correction vector is here --- adjust this as necessary.
   ; These are the same numbers as in SDSSFLUX2AB in the photoop product.
   correction = [-0.042, 0.036, 0.015, 0.013, -0.002]

   ;----------
   ; Read the file

   thisfile = (findfile(djs_filepath(plugfile, root_dir=plugdir), $
    count=ct))[0]
   if (ct NE 1) then begin
      splog, 'WARNING: Cannot find plugmap file ' + plugfile
      return, 0
   endif

   yanny_read, thisfile, pstruct, hdr=hdr, stnames=stnames, /anonymous
   if (NOT keyword_set(pstruct)) then begin
      splog, 'WARNING: Invalid plugmap file ' + thisfile
      return, 0
   endif
   plugmap = *pstruct[(where(stnames EQ 'PLUGMAPOBJ'))[0]]

   plugmap.ra = (360d0 + plugmap.ra) MOD 360d0

   ;----------
   ; Trim to object fibers only, sort them, and trim to spectrographid

   plugmap = readplugmap_sort(plugmap, hdr ,fibermask=fibermask, plates=plates)
   ;print,fibermask
   ;----------
   ; Add the tags OFFSETID and SCI_EXPTIME for 

   plugmap = struct_addtags(plugmap, $
    replicate(create_struct('OFFSETID', 0L, 'SCI_EXPTIME', 0.), $
    n_elements(plugmap)))
   i = (where(stnames EQ 'PLUGMAPPOINT', ct))[0]
   if (ct GT 0) then begin
      splog, 'Using OFFSETID and SCI_EXPTIME from PLUGMAPPOINT structure'
      plugpoint = *pstruct[i]
      for j=0L, n_elements(plugpoint)-1L do begin
         k = where(abs(plugmap.xfocal - plugpoint[j].xfocal) LT 0.0001 $
          AND abs(plugmap.yfocal - plugpoint[j].yfocal) LT 0.0001, ct)
         if (ct GT 0) then begin
            plugmap[k[0]].offsetid = plugpoint[j].offsetid
            plugmap[k[0]].sci_exptime = plugpoint[j].sci_exptime
         endif
      endfor
   endif else begin
      ; Use default values
      plugmap.offsetid = 1
      sci_exptime = 1
   endelse
   if (keyword_set(exptime)) then begin
      iuniq = uniq(plugmap.offsetid, sort(plugmap.offsetid))
      exptot = total(plugmap[iuniq].sci_exptime)
      if (exptot GT 0) then begin
         splog, 'Rescaling SCI_EXPTIME values by ', exptime/exptot
         plugmap.sci_exptime = plugmap.sci_exptime * exptime/exptot
      endif
   endif

   plateid = (yanny_par(hdr, 'plateId'))[0]
   redden_med = yanny_par(hdr, 'reddeningMed')
   if (n_elements(redden_med) NE 5) then begin
      splog, 'WARNING: Wrong number of elements for reddeningMed'
      redden_med = fltarr(5)
   endif

   ;----------
   ; Append some information from the plateHoles file

   platelist_dir = getenv('PLATELIST_DIR')
   platefile = 'plateHoles-' + string(plateid,format='(i6.6)') + '.par'
   if (keyword_set(platelist_dir)) then begin
      thisfile = (findfile(djs_filepath(platefile, $
       root_dir=platelist_dir, subdir=['plates','*','*']), count=ct))[0]
      if (ct GT 0) then begin
         plateholes = yanny_readone(thisfile, /anonymous)
         iobj = where(strmatch(plateholes.holetype,'BOSS*'))
         ;plateholes = plateholes[iobj]
         isort = lonarr(n_elements(plugmap)) - 1
         for i=0L, n_elements(plugmap)-1 do $
          isort[i] = where(plateholes.xfocal EQ plugmap[i].xfocal $
           AND plateholes.yfocal EQ plugmap[i].yfocal)
         plateholes = plateholes[isort>0]
         blankhole = plateholes[0]
         struct_assign, {junk: 0}, blankhole
         ibad = where(iobj EQ -1, nbad)
         for i=0L, nbad-1 do plateholes[ibad[i]] = blankhole
         if keyword_set(plates) then begin
           htags = ['SOURCETYPE','LAMBDA_EFF','ZOFFSET','BLUEFIBER', $
            'BOSS_TARGET*','ANCILLARY_TARGET*', 'EBOSS_TARGET*', $
            'CATALOGID','SDSSV_BOSS_TARGET*','FIRSTCARTON', $
            'GAIA_G','GAIA_BP','GAIA_RP',  $
            'RUN','RERUN','CAMCOL','FIELD','ID', 'THING_ID_TARGETING']
         endif else begin
           htags = ['SOURCETYPE','LAMBDA_EFF','ZOFFSET','BLUEFIBER', $
            'BOSS_TARGET*','ANCILLARY_TARGET*', 'EBOSS_TARGET*', $
            'RUN','RERUN','CAMCOL','FIELD','ID', 'THING_ID_TARGETING']
         endelse
         plugmap = struct_addtags(plugmap, $
          struct_selecttags(plateholes, select_tags=htags))
          
         ;- We never used washers < 175 microns
         ii = where(plugmap.zoffset lt 175)
         plugmap[ii].zoffset = 0
      
         ;- Check opfiles/washers.par for overrides to ZOFFSET
         ;;; print, "Reading washers.par"
         washers_file = getenv('IDLSPEC2D_DIR') + '/opfiles/washers.par'
         washers = yanny_readone(washers_file)

         ; extract plugmap file name to match header keyword NAME
         ; plPlugMapM-5317-56000-01.par -> 5317-56000-01
         tmp = strsplit(file_basename(plugfile), '-.', /extract)
         plugname = strjoin(tmp[1:3], '-')
         mjd = long(tmp[2])
         
         ii = where(washers.plugname eq plugname, n)
         if (n gt 1) then $
             message, "ERROR: multiple washers.par entries for " + plugname
         if (n eq 1) then begin
             status = washers[ii[0]].status
             splog, "INFO: washer ZOFFSET override ", plugname, " ", status
             if (status ne 'Y') then begin
               if (status eq 'N') then plugmap.zoffset = 0.0
               if (status eq 'L') then plugmap.zoffset = (plugmap.zoffset ne 0) * 300.0
               if (status eq 'T') then plugmap.zoffset = (plugmap.zoffset eq 300) * 300.0
               if (status eq 'X') then begin
                 splog, "WARNING: We know that we don't know ZOFFSET washer status for ", plugname
                 splog, "WARNING: setting washer ZOFFSET to default 0.0 for ", plugname
                 plugmap.zoffset = 0.0
               endif
             endif  ; status ne 'Y'
         endif else begin
             ; No explicit override; check mjd before washers were available
             ; don't print info for current plates to keep SoS quiet
             if (mjd lt 56200) then splog, "INFO: no washers.par entry for ", plugname
             if (mjd lt 55442) then begin
                 splog, "INFO: setting ZOFFSET=0 for MJD", mjd, " < 55442"
                 plugmap.zoffset = 0.0
             endif
             ; No more washers after plate 7184
             if (plateid gt 7184) then plugmap.zoffset = 0.0
         endelse
         
         if keyword_set(plates) then begin
            ;Check if plate is in list of plates to be corrected
            catfile_dir = getenv('IDLSPEC2D_DIR')+ '/catfiles/'
            catfile=catfile_dir+'Corrected_values_plate'+plateid+'_design*.fits'

            if ~FILE_TEST(catfile_dir) then $
            splog, 'WARNING: No Catalogid Correction Files found in '+catfile_dir


            if FILE_TEST(catfile) then begin
                ;loads corrected catalog file
                splog, 'Correcting Catalogid, Carton, and SDSS Magnitudes for ',plateid
                catdata=mrdfits (catfile,1)
           
                addtags = replicate(create_struct( $
                    'SDSSV_APOGEE_TARGET0', 0L, $
                    'GRI_GAIA_TRANSFORM', 0L, $
                    'org_catid',0LL), n_elements(plugmap))
                plugmap = struct_addtags(plugmap, addtags)
                plugmap.org_catid = plugmap.catalogid
                ;;; File Check
                matchlength=2.0/3600 ;Fiber Diameter in Degrees
                for j=0L, n_elements(plugmap)-1 do begin
                    spherematch, catdata.RA, catdata.DEC, plugmap[j].RA, plugmap[j].DEC, matchlength, id, match2, distance12
                    if id eq -1 then begin
                        splog, "WARNING: No Match for Correction of fiber ",j," @RA:",plugmap[j].RA," ",plugmap[j].DEC
                    endif else begin
                        if plugmap[j].catalogid ne catdata[id].Original_CatalogID then begin
                            splog, "WARNING: Catalogid Missmatch ",j," @RA DEC:",plugmap[j].RA," ",plugmap[j].DEC
                        endif else begin
                            splog, "Correcting Fiber ",STRTRIM(j,1)," Info @RA DEC:",plugmap[j].RA," ",plugmap[j].DEC
                 
                            plugmap[j].catalogid=catdata[id].Final_CatalogID
                            plugmap[j].mag[1]=catdata[id].GMAG
                            plugmap[j].mag[2]=catdata[id].RMAG
                            plugmap[j].mag[3]=catdata[id].IMAG
                            plugmap[j].mag[4]=catdata[id].ZMAG
                            plugmap[j].FIRSTCARTON=catdata[id].First_carton
                            plugmap[j].SDSSV_APOGEE_TARGET0=catdata[id].APOGEE_FLAG
                            plugmap[j].sdssv_boss_target0=catdata[id].BOSS_FLAG
                            plugmap[j].GRI_GAIA_TRANSFORM=catdata[id].TRANSFORMATION_FLAG
                            ;Instrument ;Fiber_Type
                        endelse
                    endelse
                endfor
            endif
         endif
      endif  ; ct gt 0
   endif  ; platelist_dir set
   ;----------
   ; Optionally add tags for SOS

   if (keyword_set(apotags)) then begin
      addtags = { $
       cartid   : long((yanny_par(hdr, 'cartridgeId'))[0]), $
       plateid  : long(plateid), $
       tileid   : long((yanny_par(hdr, 'tileId'))[0]), $
       raplate  : float((yanny_par(hdr, 'raCen'))[0]), $
       decplate : float((yanny_par(hdr, 'decCen'))[0]), $
       redden_med : float(redden_med), $
       fibersn    : fltarr(3), $
       synthmag   : fltarr(3) }
      plugmap = struct_addtags(plugmap, replicate(addtags, n_elements(plugmap)))
      healpix_now=0; There is no need for healpix info in the SOS 
   endif else begin
      healpix_now=0; HJIM: I decided to pass this part as an afterburner for the spAll file
   endelse
   
   if keyword_set(plates) then begin
      programname = (yanny_par(hdr, 'programname'))[0] 
      if strmatch(programname, '*eFEDS*', /fold_case) eq 1 then begin
        psffibercor = [0.7978, 0.8138, 0.8230, 0.8235]
        spht = strmatch(plugmap.FIRSTCARTON, '*bhm_spiders_clusters-efeds*')
        for i=0, n_elements(plugmap)-1 do begin
          if spht[i] then begin
            plugmap[i].mag=plugmap[i].mag-psffibercor
          endif
        endfor
      endif
      splog, 'Adding healpix info'
      addtags = replicate(create_struct( $
       'HEALPIX', 0L, $
       'HEALPIXGRP', 0, $
       'HEALPIX_DIR', ' '), n_elements(plugmap))
      plugmap = struct_addtags(plugmap, addtags)
      
      if keyword_set(healpix_now) then begin
        mwm_root='$MWM_HEALPIX';getenv('MWM_ROOT')
        run2d=getenv('RUN2D')
        ;healpix_t=plugmap.healpix
        ;healpixgrp_t=plugmap.healpixgrp
        healpix_dir_t=plugmap.healpix_dir
        ;for fid = 0L, n_elements(plugmap)-1 do begin
          ;healp=coords_to_healpix(string(plugmap[fid].ra), string(plugmap[fid].dec))
          ;healpix_t[fid]=healp.healpix
          ;healpixgrp_t[fid]=healp.healpixgrp
        healp=coords_to_healpix(string(plugmap.ra), string(plugmap.dec))
        plugmap.healpix=healp.healpix
        plugmap.healpixgrp=healp.healpixgrp
        for fid = 0L, n_elements(plugmap)-1 do begin  
          if not keyword_set(legacy) then begin
            healpix_dir_t[fid]=mwm_root + $
            strtrim(string(healp[fid].healpixgrp),2) + '/' + $
            strtrim(string(healp[fid].healpix),2) + '/boss/' + $
            strtrim(run2d, 2)+ '/' + $
            'spec-' + string(plateid,format='(i6.6)') + '-XXXX' + $;strtrim(string(mjd),2) + $
            '-' + string(plugmap[fid].catalogid,format='(i11.11)')+'.fits'
          endif  
        endfor
        ;plugmap.healpix=healpix_t
        ;plugmap.healpixgrp=healpixgrp_t
        plugmap.healpix_dir=healpix_dir_t
      endif
      
      for istd=0, n_elements(plugmap)-1 do begin
      endfor
   endif
   ;if (keyword_set(plates)) then begin
   ;    addtags = { $
   ;    targetid : long(plugmap.catalogid[0])}
   ;    plugmap = struct_addtags(plugmap, replicate(addtags, n_elements(plugmap)))
   ;    plugmap.targetid=plugmap.catalogid
   ;endif
   ;----------
   ; Read calibObj or photoPlate photometry data
   if (keyword_set(calibobj)) then begin
      splog, 'Adding fields from calibObj file'
      addtags = replicate(create_struct( $
       'CALIBFLUX', fltarr(5), $
       'CALIBFLUX_IVAR', fltarr(5), $
       'CALIB_STATUS', lonarr(5), $
       'SFD_EBV', 0., $
       'WISE_MAG', fltarr(4), $
       'TWOMASS_MAG', fltarr(3), $ 
       'GUVCAT_MAG', fltarr(2), $
       'GAIA_PARALLAX', 0.0, $
       'GAIA_PMRA', 0.0, $
       'GAIA_PMDEC', 0.0), n_elements(plugmap))
      plugmap = struct_addtags(plugmap, addtags)
      ra_temp=plugmap.ra
      dec_temp=plugmap.dec
      ;----------
      ;Obtain the WISE, TWOMASS, GUVCAT and GAIA pm
      wise_temp=fltarr(4,n_elements(plugmap))
      two_temp=fltarr(3,n_elements(plugmap))
      guv_temp=fltarr(2,n_elements(plugmap))
      parallax_temp=fltarr(n_elements(plugmap))
      pmra_temp=fltarr(n_elements(plugmap))
      pmdec_temp=fltarr(n_elements(plugmap))
      
      splog, "Obtaing the WISE, TWOMASS, GUVCAT and GAIA parallax and pm"
      openw,lun,'catalog.inp',/get_lun
      for istd=0, n_elements(ra_temp)-1 do begin
          ;cmd = "catalogdb_photo "+strtrim(string(ra_temp[istd]),2)+" "+strtrim(string(dec_temp[istd]),2)
          ;spawn, cmd, dat
          printf,lun,strtrim(string(ra_temp[istd]),2)+" "+strtrim(string(dec_temp[istd]),2)
      endfor
      free_lun, lun
      cmd = "catalogdb_photo_file catalog.inp"
      spawn, cmd, alldat
      for istd=0, n_elements(ra_temp)-1 do begin
        dat=alldat[istd]
        splog, "Photometric Data for fiber "+string(istd+1)+": "+dat
        tp=strsplit(dat,' ',/extract)
        wise_temp[0,istd]=float(tp[0])
        wise_temp[1,istd]=float(tp[1])
        wise_temp[2,istd]=float(tp[2])
        wise_temp[3,istd]=float(tp[3])
        two_temp[0,istd]=float(tp[4])
        two_temp[1,istd]=float(tp[5])
        two_temp[2,istd]=float(tp[6])
        guv_temp[0,istd]=float(tp[7])
        guv_temp[1,istd]=float(tp[8])
        parallax_temp[istd]=float(tp[9])
        pmra_temp[istd]=float(tp[10])
        pmdec_temp[istd]=float(tp[11])
      endfor
      if FILE_TEST('catalog.inp') then FILE_DELETE, 'catalog.inp', /QUIET
      plugmap.wise_mag=wise_temp
      plugmap.twomass_mag=two_temp
      plugmap.guvcat_mag=guv_temp
      plugmap.gaia_parallax=parallax_temp
      plugmap.gaia_pmra=pmra_temp
      plugmap.gaia_pmdec=pmdec_temp
      ;----------
      ; Read the SFD dust maps
      euler, plugmap.ra, plugmap.dec, ll, bb, 1
      plugmap.sfd_ebv = dust_getval(ll, bb, /interp)
      ;---------
      ;Redefine the Extintion using the RJCE extintion method, see Majewski, Zasowski & Nidever (2011) and Zasowski et al. (2013)
      rjce_extintion=0
      if keyword_set(rjce_extintion) then begin
        if keyword_set(plates) then begin
          programname = (yanny_par(hdr, 'programname'))[0] 
          if strmatch(programname, '*MWM*', /fold_case) eq 1 then begin
            spht = strmatch(plugmap.objtype, 'SPECTROPHOTO_STD')
            ispht = where(spht, nspht)
            stsph=plugmap(ispht)
            catid=stsph.catalogid
            ebv_std=stsph.sfd_ebv
            fib_std=stsph.fiberId
            splog, "RJCE extintion"
            for istd=0, n_elements(catid)-1 do begin
              cmd = "catalogdb_ev "+strtrim(string(catid[istd]),2)
              spawn, cmd, dat
              dat=double(dat)
              if (finite(dat) ne 0) and (dat gt 0) and (dat le 1.2*ebv_std[istd]) then begin
                splog,"change SFD E(B-V) "+strtrim(string(ebv_std[istd]),2)+" by RJCE E(B-V) " + $
                strtrim(string(dat),2)+" on SPECTROPHOTO_STD fiber "+strtrim(string(fib_std[istd]),2) + $
                " with CATALOGID "+strtrim(string(catid[istd]),2)
                ebv_std[istd]=dat
              endif
            endfor
            stsph.sfd_ebv=ebv_std
            plugmap(ispht)=stsph
            splog, "done with RJCE extintion"
            gaiaext=0
          endif
        endif
      endif
      if keyword_set(MWM_fluxer) then begin
        if keyword_set(plates) then begin
          programname = (yanny_par(hdr, 'programname'))[0]
          if ((strmatch(programname, '*MWM*', /fold_case) eq 1) $
           || (strmatch(programname, '*OFFSET*', /fold_case) eq 1)) then begin
           gaiaext = 1
          endif
        endif
      endif
      if keyword_set(gaiaext) then begin
        ;---------
        ;Redefine the Extintion using the Bayestar 3D dust extintion maps
        spht = strmatch(plugmap.objtype, 'SPECTROPHOTO_STD')
        ispht = where(spht, nspht)
        stsph=plugmap(ispht)
        ll_std=ll[ispht]
        bb_std=bb[ispht]
        ra_plate=float(yanny_par(hdr, 'raCen'))
        dec_plate=float(yanny_par(hdr, 'decCen'))
        rm_read_gaia, ra_plate,dec_plate,stsph,dist_std=dist_std
        ;ra_std=stsph.ra
        ;dec_std=stsph.dec
        ebv_std=stsph.sfd_ebv
        fib_std=stsph.fiberId
        splog, "running  dust_3d_map"
        ll_std_str='['+strjoin(strtrim(ll_std,1),',')+']'
        bb_std_str='['+strjoin(strtrim(bb_std,1),',')+']'
        dist_std_str='['+strjoin(strtrim(dist_std,1),',')+']'
        cmd = "dust_3d_map_arr.py "+ll_std_str+" "+bb_std_str+" "+dist_std_str
        spawn, cmd, dat_arr
        remchar, dat_arr, '['
        remchar, dat_arr, ']'
        dat_arr=double(STRSPLIT(strjoin(dat_arr,' '),/EXTRACT))
        for istd=0, n_elements(dist_std)-1 do begin
            dat=dat_arr[istd]
            if (finite(dat) ne 0) and (dat le ebv_std[istd]) then begin
                splog,"change E(B-V) "+strtrim(string(ebv_std[istd]),2)+" by " + $
                 strtrim(string(dat),2)+" on SPECTROPHOTO_STD fiber "+strtrim(string(fib_std[istd]),2) + $
                 "; Gaia DR2 parallax distance: "+strtrim(string(dist_std[istd]),2)+" pc"
                ebv_std[istd]=dat
            endif
        endfor
;        for istd=0, n_elements(dist_std)-1 do begin
;          cmd = "dust_3d_map.py "+strtrim(string(ll_std[istd]),2)+" "+strtrim(string(bb_std[istd]),2)+" "+strtrim(string(dist_std[istd]),2)
;          ;print, cmd
;          ;spawn, cmd
;          spawn, cmd, dat
;          dat=double(dat)
;          ;print,dust_getval(ll_std[istd], bb_std[istd], /interp)
;          ;print,ll_std[istd], bb_std[istd]
;          if (finite(dat) ne 0) and (dat le ebv_std[istd]) then begin
;            splog,"change E(B-V) "+strtrim(string(ebv_std[istd]),2)+" by " + $
;            strtrim(string(dat),2)+" on SPECTROPHOTO_STD fiber "+strtrim(string(fib_std[istd]),2) + $
;             "; Gaia DR2 parallax distance: "+strtrim(string(dist_std[istd]),2)+" pc"
;            ebv_std[istd]=dat
;          endif
;        endfor
        stsph.sfd_ebv=ebv_std
        plugmap(ispht)=stsph
        splog, "done with 3D dust matches"
      endif
      
      ;----------
      ; Attempt to read the calibObj photometry data
      
      tsobj = plug2tsobj(plateid, /plates, /legacy, _EXTRA=KeywordsForPhoto)

      ; Do not use the calibObj structure if more than 20% of the non-sky
      ; objects do not have fluxes.
      if (keyword_set(tsobj)) then begin
         qexist = tsobj.psfflux[2] NE 0
         qsky = strmatch(plugmap.objtype,'SKY*')
         splog, 'Matched ', fix(total(qsky EQ 0 AND qexist)), $
          ' of ', fix(total(qsky EQ 0)), ' non-SKY objects'
         if (total((qsky EQ 0) AND qexist) LT 0.80*total(qsky EQ 0)) then begin
            splog, 'Discarding calibObj structure because < 80% matches'
            tsobj = 0
         endif
      endif

      if (keyword_set(tsobj)) then begin

         ; Propagate CALIB_STATUS information:
         if tag_exist(tsobj, 'CALIB_STATUS') then $
            plugmap.calib_status = tsobj.calib_status

         ; Assume that all objects not called a 'GALAXY' are stellar objects
         qstar = strmatch(plugmap.objtype, 'GALAXY*') EQ 0
         istar = where(qstar AND qexist, nstar)
         igal = where(qstar EQ 0 AND qexist, ngal)
         if (tag_exist(tsobj,'FIBER2FLUX')) then begin
            fiberflux = transpose(tsobj.fiber2flux)
            fiberflux_ivar = transpose(tsobj.fiber2flux_ivar)
            pratio = [2.085, 2.085, 2.116, 2.134, 2.135]
         endif else begin
            fiberflux = transpose(tsobj.fiberflux)
            fiberflux_ivar = transpose(tsobj.fiberflux_ivar)
            pratio = [1.343, 1.336, 1.354, 1.363, 1.367]
         endelse
         if (nstar GT 0) then begin
            plugmap[istar].calibflux = tsobj[istar].psfflux
            plugmap[istar].calibflux_ivar = tsobj[istar].psfflux_ivar
;            ; Compute the ratio of PSF/FIBER flux for stars in each filter,
;            ; using only stars that are brighter than 30 nMgy (= 18.8 mag).
;            ; If no such stars, then this ratio is set to unity.
;            for ifilt=0, 4 do begin
;               v1 = tsobj[istar].psfflux[ifilt]
;               v2 = fiberflux[istar,ifilt]
;               jj = where(v1 GT 30 AND v2 GT 30, ct)
;               if (ct GT 0) then pratio[ifilt] = median([ v1[jj] / v2[jj] ])
;            endfor
         endif
         splog, 'PSF/fiber flux ratios = ', pratio
         if (ngal GT 0) then begin
            for ifilt=0, 4 do begin
               plugmap[igal].calibflux[ifilt] = $
                fiberflux[igal,ifilt] * pratio[ifilt]
               plugmap[igal].calibflux_ivar[ifilt] = $
                fiberflux_ivar[igal,ifilt] / (pratio[ifilt])^2
            endfor
         endif

         ; Reject any fluxes based upon suspect PHOTO measurements,
         ; as indicated by the PHOTO flags.
         badbits2 = sdss_flagval('OBJECT2','SATUR_CENTER') $
          OR sdss_flagval('OBJECT2','INTERP_CENTER') $
          OR sdss_flagval('OBJECT2','PSF_FLUX_INTERP')
         qgoodphot = (tsobj.flags2 AND badbits2) EQ 0
         plugmap.calibflux = plugmap.calibflux * qgoodphot
         plugmap.calibflux_ivar = plugmap.calibflux_ivar * qgoodphot
      endif else begin
         splog, 'WARNING: No calibObj structure found for plate ', plateid
      endelse

      ;----------
      ; For any objects that do not have photometry from the calibObj
      ; structure, simply translate the flux from the plugmap MAG values
      ; (as long as those values are in the range 0 < MAG < +50).

      for ifilt=0, 4 do begin
         pratio = [2.085, 2.085, 2.116, 2.134, 2.135];ratio of fiber2flux
         splog, 'PSF/fiber flux ratios = ', pratio
         ibad = where(plugmap.calibflux[ifilt] EQ 0 $
          AND plugmap.mag[ifilt] GT 0 $
          AND plugmap.mag[ifilt] LT 50, nbad)
         if (nbad GT 0) then begin
            splog, 'Using plug-map fluxes for ', nbad, $
             ' values in filter ', ifilt
            plugmap[ibad].calibflux[ifilt] = $
             10.^((22.5 - plugmap[ibad].mag[ifilt]) / 2.5)*pratio[ifilt]
            plugmap[ibad].calibflux_ivar[ifilt] = 0
         endif
      endfor

      ;----------
      ; Apply AB corrections to the CALIBFLUX values (but not to MAG)

      factor = exp(-correction/2.5 * alog(10))
      for j=0,4 do plugmap.calibflux[j] = plugmap.calibflux[j] * factor[j]
      for j=0,4 do $
       plugmap.calibflux_ivar[j] = plugmap.calibflux_ivar[j] / factor[j]^2
   endif

   ;-----
   ; The following lines are modified to apply the new extinctiion coefficients.
   ; They are just the ratios of new/old,  
   redden_corr=[0.822308,0.870815,0.830607,0.813998,0.853955]	
   if (keyword_set(deredden)) then begin
      splog, 'Applying reddening vector ', redden_med
      for ifilt=0, 4 do begin
         if (plateid ge 7572) then begin
            plugmap.mag[ifilt] = plugmap.mag[ifilt] - redden_med[ifilt]
         endif else begin
            plugmap.mag[ifilt] = plugmap.mag[ifilt] - redden_med[ifilt]*redden_corr[ifilt]
            splog, 'modified extinction co-efficients: ',redden_med[ifilt]*redden_corr[ifilt]
         endelse
      endfor	
   endif


   ; Optionally trim to selected spectrograph
   if keyword_set(plates) then begin
      nfiber = 1000;n_elements(plugmap)
      plugmap.spectrographid=1
   endif else begin
      nfiber = n_elements(plugmap)
   endelse
   if (keyword_set(spectrographid)) then begin
      print,spectrographid,nfiber
      indx = (spectrographid-1)*nfiber/2 + lindgen(nfiber/2)
      plugmap = plugmap[indx]
      fibermask = fibermask[indx]
   endif
   ;print,nfiber
   ;print,plugmap.fiberid
   ;print,fibermask
   ;print,(size(plugmap, /dimens))[0]
   ;exit
   
   if tag_exist(plugmap,'ORG_FIBERID') || tag_exist(plugmap,'org_catid') then begin
;     struct_print, plugmap, filename='plugmap_'+STRTRIM(mjd,1)+'_correcting.html', /html
     plugmap = struct_trimtags(plugmap,except_tags=['ORG_FIBERID','org_catid'])
   endif
;   struct_print, plugmap, filename='plugmap_'+STRTRIM(mjd,1)+'.html', /html

   return, plugmap
end
;------------------------------------------------------------------------------
