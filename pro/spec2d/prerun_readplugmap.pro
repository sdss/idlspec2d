;+
; NAME:
;   prerun_readplugmap
;
; PURPOSE:
;   Read plugmap file prior to spreduce2d and preforms all crossmatch calcuations, writing correction values to file
;
; CALLING SEQUENCE:
;   plugmap = prerun_readplugmap( plugfile, [ plugdir=, $
;    /apotags, exptime=, hdr=, fibermask=, _EXTRA= ] )
;
; INPUTS:
;   plugfile  - Name of Yanny-parameter plugmap file
;
; OPTIONAL INPUTS:
;   plugdir   - Directory for PLUGFILE
;   apotags   - If set, then add a number of tags to the output structure
;               constructed from the Yanny header.  These tags are:
;               CARTID, PLATEID, TILEID, RAPLATE, DECPLATE, REDDEN_MED.
;               Also add the tags FIBERSN[3], SYTHMAG[3] which are used
;               by the on-the-mountain reductions.
;   exptime   - Default exposure time SCI_EXPTIME to add to the output;
;               if there are multiple pointings and EXPTIME is set, then
;               the exposure time in each pointing is scaled such that
;               their sum is EXPTIME.
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
; INTERNAL FUNCTIONS CALLED:
;   calibrobj
;   ApogeeToBOSSRobomap
;   BossToApogeeRobomap
;   readFPSobsSummary
;   readPlateplugmap_sort
;   readPlateplugMap
;
; REVISION HISTORY:
;   29-Jan-2001  Written by S. Burles, FNAL
;   07-Aug-2012  Added ZOFFSET overrides; S. Bailey, LBL
;   21-Jun-2021  Editted by S Morrison to add correction check for catalog
;   15-Oct-2021  Modified by S Morrison to prepare for FPS:
;                       -adding readFPSobsSummary, ApogeeToBOSSRobomap, BossToApogeeRobomap
;                       -breaking readPlateplugMap in to a new function to preseve for FPS
;-
;------------------------------------------------------------------------------
function calibrobj, plugfile, fibermap, fieldid, rafield, decfield, programname=programname,$
                    plates=plates, legacy=legacy, fps=fps, MWM_fluxer=MWM_fluxer,$
                    KeywordsForPhoto=KeywordsForPhoto

    ; The correction vector is here --- adjust this as necessary.
    ; These are the same numbers as in SDSSFLUX2AB in the photoop product.
    correction = [-0.042, 0.036, 0.015, 0.013, -0.002]

    splog, 'Adding fields from calibObj file'
    addtags = replicate(create_struct( $
            'CALIBFLUX', fltarr(5), $
            'CALIBFLUX_IVAR', fltarr(5), $
            'CALIB_STATUS', lonarr(5), $
            'SFD_EBV', 0., $
            'SFD_EBV_gaia', 0., $
            'SFD_EBV_RJCE', 0., $
            'WISE_MAG', fltarr(4), $
            'TWOMASS_MAG', fltarr(3), $
            'GUVCAT_MAG', fltarr(2), $
            'GAIA_PARALLAX', 0.0, $
            'GAIA_PMRA', 0.0, $
            'GAIA_PMDEC', 0.0), n_elements(fibermap))

    fibermap = struct_addtags(fibermap, addtags)


    splog, "Running 'run_plugmap_supplements.py' to retreive supplementary info:"
    flags=" --mags --rjce --gaia"
    if keyword_set(logfile) then begin
        flags = flags+" --log "+ logfile
    endif
    catfile=FILE_BASENAME(plugfile,'.par')+'.inp'
    supfile=FILE_BASENAME(plugfile,'.par')+'_supp.fits'
    
    ra_temp=fibermap.ra
    dec_temp=fibermap.dec
    catid_temp=fibermap.catalogid
    if keyword_set(fps) then begin
        spht = strmatch(fibermap.program, 'ops_std', /FOLD_CASE)
    endif else begin
        spht = strmatch(fibermap.objtype, 'SPECTROPHOTO_STD', /FOLD_CASE)
    endelse
    stdflag=BYTARR(n_elements(ra_temp))
    ispht = where(spht, nspht)
    stdflag[ispht]=1
    euler, fibermap.ra, fibermap.dec, ll, bb, 1
    stsph=fibermap(ispht)
    dist=DBLARR(n_elements(ra_temp))
    stsph.RA[where(abs(stsph.RA) gt 90)] = stsph.RAcat[where(abs(stsph.RA) gt 90)]
    stsph.dec[where(abs(stsph.dec) gt 90)] = stsph.deccat[where(abs(stsph.dec) gt 90)]
    if abs(stsph.DEC)gt 90 then stsph.DEC = 0.0
    if abs(stsph.RA)gt 90 then stsph.RA = 0.0
    if abs(decfield)gt 90 then decfield = 0.0
    if abs(rafield)gt 90 then rafield = 0.0
    rm_read_gaia, rafield,decfield,stsph,dist_std=dist_std
    dist[ispht]=dist_std

    openw,lun1,catfile,/get_lun
    for istd=0, n_elements(ra_temp)-1 do begin
        printf,lun1,strtrim(string(ra_temp[istd]),2)+$
            " "+strtrim(string(dec_temp[istd]),2)+$
            " "+string(catid_temp[istd])+" "+string(stdflag[istd],/PRINT)+$
            " "+string(ll[istd])+" "+string(bb[istd])+" "+string(dist[istd])
    endfor
    free_lun, lun1
    
    ;Obtain the WISE, TWOMASS, GUVCAT and GAIA pm
    wise_temp=fltarr(4,n_elements(fibermap))
    two_temp=fltarr(3,n_elements(fibermap))
    guv_temp=fltarr(2,n_elements(fibermap))
    parallax_temp=fltarr(n_elements(fibermap))
    pmra_temp=fltarr(n_elements(fibermap))
    pmdec_temp=fltarr(n_elements(fibermap))

    cmd = "run_plugmap_supplements.py " + catfile + flags
    FILE_DELETE, supfile, /ALLOW_NONEXISTENT

    
    while not FILE_TEST(supfile) DO spawn, cmd, dat
    if not keyword_set(logfile) then $
        foreach row, dat do splog, row
    supplements = mrdfits(supfile,1)
    wise_temp[0,*]=supplements.w1mpro
    wise_temp[1,*]=supplements.w2mpro
    wise_temp[2,*]=supplements.w3mpro
    wise_temp[3,*]=supplements.w4mpro
    two_temp[0,*]=supplements.j2mass
    two_temp[1,*]=supplements.h2mass
    two_temp[2,*]=supplements.k2mass
    guv_temp[0,*]=supplements.fuv
    guv_temp[1,*]=supplements.nuv
    parallax_temp=supplements.parallax
    pmra_temp=supplements.pmra
    pmdec_temp=supplements.pmdec

    fibermap.wise_mag=wise_temp
    fibermap.twomass_mag=two_temp
    fibermap.guvcat_mag=guv_temp
    fibermap.gaia_parallax=parallax_temp
    fibermap.gaia_pmra=pmra_temp
    fibermap.gaia_pmdec=pmdec_temp

    ; Read the SFD dust maps
    euler, fibermap.ra, fibermap.dec, ll, bb, 1
    fibermap.sfd_ebv = dust_getval(ll, bb, /interp)
    sfd_ebv_RJCE = fibermap.sfd_ebv
    sfd_ebv_gaia = fibermap.sfd_ebv

    ;---------
    ;Redefine the Extintion using the RJCE extintion method, see Majewski, Zasowski & Nidever (2011) and Zasowski et al. (2013)
    sfd_ebv_RJCE[ispht] = supplements[ispht].EB_RJCE
    fibermap.sfd_ebv_RJCE=sfd_ebv_RJCE

    ;Redefine the Extintion using the Bayestar 3D dust extintion maps
    sfd_ebv_gaia[ispht] = supplements[ispht].REDDENING_GAIA
    fibermap.sfd_ebv_gaia=sfd_ebv_gaia
    
    ;----------
    ; Attempt to read the calibObj photometry data
    if keyword_set(legacy) then begin
      tsobj = plug2tsobj(fieldid, plates=plates, legacy=legacy, _EXTRA=KeywordsForPhoto)
      ; Do not use the calibObj structure if more than 20% of the non-sky
      ; objects do not have fluxes.
      if (keyword_set(tsobj)) then begin
          qexist = tsobj.psfflux[2] NE 0
          qsky = strmatch(fibermap.objtype,'SKY*')
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
              fibermap.calib_status = tsobj.calib_status
          ; Assume that all objects not called a 'GALAXY' are stellar objects
          qstar = strmatch(fibermap.objtype, 'GALAXY*') EQ 0
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
              fibermap[istar].calibflux = tsobj[istar].psfflux
              fibermap[istar].calibflux_ivar = tsobj[istar].psfflux_ivar
          endif
          splog, 'PSF/fiber flux ratios = ', pratio
          if (ngal GT 0) then begin
              for ifilt=0, 4 do begin
                  fibermap[igal].calibflux[ifilt] = $
                      fiberflux[igal,ifilt] * pratio[ifilt]
                  fibermap[igal].calibflux_ivar[ifilt] = $
                      fiberflux_ivar[igal,ifilt] / (pratio[ifilt])^2
              endfor
          endif
          ; Reject any fluxes based upon suspect PHOTO measurements,
          ; as indicated by the PHOTO flags.
          badbits2 = sdss_flagval('OBJECT2','SATUR_CENTER') $
                  OR sdss_flagval('OBJECT2','INTERP_CENTER') $
                  OR sdss_flagval('OBJECT2','PSF_FLUX_INTERP')
          qgoodphot = (tsobj.flags2 AND badbits2) EQ 0
          fibermap.calibflux = fibermap.calibflux * qgoodphot
          fibermap.calibflux_ivar = fibermap.calibflux_ivar * qgoodphot
      endif else begin
          if not keyword_set(fps) then begin
              splog, 'WARNING: No calibObj structure found for plate ', fieldid
          endif else splog, 'WARNING: No calibObj structure found for field ', fieldid
      endelse
    endif
    ;----------
    ; For any objects that do not have photometry from the calibObj
    ; structure, simply translate the flux from the fibermap MAG values
    ; (as long as those values are in the range 0 < MAG < +50).

    for ifilt=0, 4 do begin
       pratio = [2.085, 2.085, 2.116, 2.134, 2.135];ratio of fiber2flux
       splog, 'PSF/fiber flux ratios = ', pratio
       ibad = where(fibermap.calibflux[ifilt] EQ 0 $
                      AND fibermap.mag[ifilt] GT 0 $
                      AND fibermap.mag[ifilt] LT 50, nbad)
       if (nbad GT 0) then begin
          splog, 'Using plug-map fluxes for ', nbad, $
                 ' values in filter ', ifilt
          fibermap[ibad].calibflux[ifilt] = $
              10.^((22.5 - fibermap[ibad].mag[ifilt]) / 2.5)*pratio[ifilt]
          fibermap[ibad].calibflux_ivar[ifilt] = 0
       endif
    endfor

    ;----------
    ; Apply AB corrections to the CALIBFLUX values (but not to MAG)
    factor = exp(-correction/2.5 * alog(10))
    for j=0,4 do fibermap.calibflux[j] = fibermap.calibflux[j] * factor[j]
    for j=0,4 do $
      fibermap.calibflux_ivar[j] = fibermap.calibflux_ivar[j] / factor[j]^2
      
    FILE_DELETE, catfile, /ALLOW_NONEXISTENT
    FILE_DELETE, supfile, /ALLOW_NONEXISTENT
  return, fibermap
end

;------------------------------------------------------------------------------

function rename_tags, struct, oldtags, newtags
  tags = tag_names(struct)
  FOR i=0L, n_elements(tags)-1 DO BEGIN
    w=where(oldtags EQ tags[i], nw)
    IF nw NE 0 THEN taguse = newtags[w[0]] ELSE taguse = tags[i]
    IF i EQ 0 THEN newst = create_struct(taguse, struct[0].(i)) $
    ELSE newst = create_struct(newst, taguse, struct[0].(i))
  ENDFOR
  newstruct = replicate(newst, n_elements(struct))
  FOR i=0L, n_elements(tags)-1 DO newstruct.(i) = struct.(i)
  return, newstruct
END

function ApogeeToBOSSRobomap, robomap
  iassigned_apogee = where(robomap.spectrographId eq 0 AND robomap.Assigned EQ 1)
;  iassigned_apogee = where(robomap.spectrographId eq 2 AND robomap.Assigned EQ 1)
  uassigned_boss = where(robomap.spectrographId eq 1 AND robomap.Assigned EQ 0)

  boss = where(robomap.spectrographID eq 1)

  Assigned_Apogee_fiberIDs = robomap[iassigned_apogee].POSITIONERID
  paired_BOSS_fiberIDs = Assigned_Apogee_fiberIDs

  match, robomap[boss].POSITIONERID, paired_BOSS_fiberIDs, matched_fiber, paired
  catid_apogee = robomap[iassigned_apogee[paired]].CATALOGID
  catid_apogee='u'+(strtrim(catid_apogee,1))

  ;robomap[boss][matched_fiber].CATALOGID= string(robomap[[matched_fiber].CATALOGID)

  robomap_boss=robomap[boss]
  robomap_boss[matched_fiber].CATALOGID = catid_apogee
  robomap[boss]=robomap_boss


  return, robomap
end

function BossToApogeeRobomap, robomap
  iassigned_boss = where(robomap.spectrographId eq 1 AND robomap.Assigned EQ 1)
  uassigned_apogee = where(robomap.spectrographId eq 0 AND robomap.Assigned EQ 0)
;  uassigned_apogee = where(robomap.spectrographId eq 2 AND robomap.Assigned EQ 0)

  apogee = where(robomap.spectrographID eq 0)
;  apogee = where(robomap.spectrographID eq 2)


  Assigned_BOSS_fiberIDs = robomap[iassigned_boss].POSITIONERID
  paired_Apogee_fiberIDs = Assigned_BOSS_fiberIDs

  match, robomap[apogee].POSITIONERID, paired_Apogee_fiberIDs, matched_fiber, paired

  catid_BOSS = robomap[iassigned_boss[paired]].CATALOGID
  catid_BOSS='u'+(strtrim(catid_BOSS,1))

  robomap_apogee=robomap[apogee]
  robomap_apogee[matched_fiber].CATALOGID = catid_BOSS
  robomap[apogee]=robomap_apogee

  return, robomap
end

function NoAssignRobomap, robomap

  sign = MAKE_ARRAY(n_elements(robomap.DEC), /STRING, VALUE = '+')
  sign[where(robomap.DEC lt 0)] = '-'

  radec, robomap.RA, robomap.DEC, ihr, imin, xsec, ideg, imn, xsc
  xsc_str=string(xsc,format='(f4.1)')
  xsc_str[where(xsc lt 10)]='0'+string(xsc_str[where(xsc lt 10)],format='(f3.1)')

  xsec_str=string(xsec,format='(f4.1)')
  xsec_str[where(xsec lt 10)]='0'+string(xsc_str[where(xsec lt 10)],format='(f3.1)')

  dummy_catid = 'u'+strtrim(string(ihr,format='(I2.2)'),2)+$
    string(imin,format='(I2.2)')+xsec_str+$
    sign+strtrim(string(abs(ideg),format='(I3.3)'),2)+$
    string(imn,format='(I2.2)')+xsc_str

  All_apogee_fibers = where(robomap.spectrographId eq 0)
;  All_apogee_fibers = where(robomap.spectrographId eq 2)

  iunassigned_boss = where(robomap.spectrographId eq 1 AND robomap.Assigned EQ 0)
  iunassigned_apogee = where(robomap.spectrographId eq 0 AND robomap.Assigned EQ 0)
;  iunassigned_apogee = where(robomap.spectrographId eq 2 AND robomap.Assigned EQ 0)

  match, robomap[iunassigned_boss].POSITIONERID, robomap[iunassigned_apogee].POSITIONERID, matched_fiber, paired
  robomap_boss=robomap[iunassigned_boss]
  boss_dummpycat=dummy_catid[iunassigned_boss]
  robomap_boss[matched_fiber].CATALOGID = boss_dummpycat[matched_fiber]
  robomap[iunassigned_boss]=robomap_boss

  match, robomap[iunassigned_apogee].POSITIONERID, robomap[iunassigned_boss].POSITIONERID, matched_fiber, paired
  robomap_apogee=robomap[iunassigned_apogee]
  apogee_dummpycat=dummy_catid[iunassigned_apogee]
  robomap_apogee[matched_fiber].CATALOGID = apogee_dummpycat[matched_fiber]
  robomap[iunassigned_apogee]=robomap_apogee

  foreach boss_fiber, iunassigned_boss do begin
    ind = where(boss_fiber eq All_apogee_fibers, count)
    if count eq 0 then robomap[boss_fiber].CATALOGID = dummy_catid[boss_fiber]
  endforeach

  return, robomap
end

;------------------------------------------------------------------------------

function readFPSobsSummary, plugfile, robomap, stnames, mjd, hdr=hdr, $
            apotags=apotags, exptime=exptime, fibermask=fibermask, $
            _EXTRA=KeywordsForPhoto

   ; The correction vector is here --- adjust this as necessary.
   ; These are the same numbers as in SDSSFLUX2AB in the photoop product.
   correction = [-0.042, 0.036, 0.015, 0.013, -0.002]

   if not TAG_EXIST(robomap,'fiberid') then begin
        robomap = struct_addtags(robomap, $
            replicate(create_struct('fiberid', 0L), n_elements(robomap)))
        robomap.fiberid=robomap.POSITIONERID
   endif

   nfiber=n_elements(robomap)

   robomap = rename_tags(robomap, 'CATALOGID','iCATALOGID')
   robomap = struct_addtags(robomap, $
     replicate(create_struct('CATALOGID', ''), n_elements(robomap)))
   robomap.CATALOGID=strtrim(robomap.iCATALOGID,1)
   robomap = ApogeeToBOSSRobomap(robomap)
   robomap = BossToApogeeRobomap(robomap)
   robomap = NoAssignRobomap(robomap)

   if (NOT keyword_set(fibermask)) then fibermask = bytarr(nfiber) $
   else if (n_elements(fibermask) NE nfiber) then $
    message, 'Number of elements in FIBERMASK do not match NFIBER'


   robomap = struct_addtags(robomap, $
        replicate(create_struct('BADSTDMASK', 0L), n_elements(robomap)))

   robomap = struct_addtags(robomap, $
        replicate(create_struct('fibermask', 0L), n_elements(robomap)))

   robomap = struct_addtags(robomap, $
        replicate(create_struct('OBJTYPE', ''), n_elements(robomap)))
   robomap.OBJTYPE = robomap.category

   robomap = rename_tags(robomap, 'ra','ra_obs')
   robomap = rename_tags(robomap, 'dec','dec_obs')

   robomap = rename_tags(robomap, 'racat','ra')
   robomap = rename_tags(robomap, 'deccat','dec')

   iunAssigned = where(robomap.Assigned EQ 0, nAssigned)

   fibermask = fibermask OR fibermask_bits('NOPLUG')
   fibermask[iunAssigned] = fibermask[iunAssigned] - fibermask_bits('NOPLUG') ; TODO: NEED TO UPDATE with new fibermask bit value


   if (keyword_set(apotags)) then begin
        addtags = { configuration_id    :   long((yanny_par(hdr, 'configuration_id'))[0]), $
                    targeting_vers      :   long((yanny_par(hdr, 'targeting_version'))[0]), $
                    observation_id      :   long((yanny_par(hdr, 'observation_id'))[0]), $
                    fieldid             :   long((yanny_par(hdr, 'fieldid'))[0]), $
                    MJD                 :   long((yanny_par(hdr, 'MJD'))[0]), $
                    rafield             :   float((yanny_par(hdr, 'RACEN'))[0]), $
                    decfield            :   float((yanny_par(hdr, 'DECCEN'))[0]), $
                    redden_med          :   float((yanny_par(hdr, 'reddeningMed'))[0]), $
                    fibersn             :   fltarr(3), $
                    synthmag            :   fltarr(3), $
                    hrmed               :   float((yanny_par(hdr, 'haMed'))[0])$
                  }
      robomap = struct_addtags(robomap, replicate(addtags, n_elements(robomap)))
      healpix_now=0; There is no need for healpix info in the SOS
   endif else begin
      healpix_now=0; HJIM: I decided to pass this part as an afterburner for the spAll file
   endelse

   mjd=long((yanny_par(hdr, 'MJD'))[0])
   ;----------
   if not keyword_set(apotags) then begin
        ; Read calibObj or photoPlate photometry data
        fieldid = (yanny_par(hdr, 'field_id'))[0]
        ra_field=float(yanny_par(hdr, 'raCen'))
        dec_field=float(yanny_par(hdr, 'decCen'))
        robomap=calibrobj(plugfile,robomap, fieldid, ra_field, dec_field, /fps, $
                          KeywordsForPhoto=KeywordsForPhoto)
   endif
   ;-----
    
   ;robomap = struct_addtags(robomap, replicate({orig_objtype:''}, n_elements(robomap)))
   ;robomap.orig_objtype=robomap.objtype
   ;print,robomap.orig_objtype
   robomap[where(strtrim(robomap.objtype,2) EQ 'sky_boss')].objtype = 'SKY'
   robomap.fibermask=fibermask
   return, robomap
end

;------------------------------------------------------------------------------


function readPlateplugmap_sort, plugmap, hdr, fibermask=fibermask, plates=plates

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

   plugsort = struct_addtags(plugsort, $
   replicate(create_struct('BADSTDMASK', 0L), n_elements(plugsort)))
   plugmap = struct_addtags(plugmap, $
   replicate(create_struct('BADSTDMASK', 0L), n_elements(plugmap)))

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
               badstdmask[i] = 1
               plugmap[i].badstdmask = 1
            endif
          endif
        endfor
      endif
      igood = where(qobj AND plugmap.fiberid GT 0 AND plugmap.spectrographid EQ 1 AND badstdmask EQ 0, ngood); check the number -1
      iplugged = where(qobj AND plugmap.fiberid GT 0 AND plugmap.spectrographid EQ 1, nplugged); check the number -1
      igoodapoge = where(qobj AND plugmap.fiberid GT 0 AND plugmap.spectrographid EQ 2 AND badstdmask EQ 0, napogee)
   endif else begin
      igood = where(qobj AND plugmap.fiberid GT 0, ngood)
      iplugged = igood
   endelse
   if (ngood EQ 0) then $
    message, 'No fibers found in plugmap!'
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
      plugsort[imissing].spectrographid = 1
   endif
   ;print,nmissing,nfill,ngood,nfiber
   ;print,plugsort.fiberid
   return, plugsort
end

;------------------------------------------------------------------------------

function readPlateplugMap, plugfile, plugmap,stnames, spectrographid, mjd, $
    apotags=apotags, exptime=exptime, $
    hdr=hdr, fibermask=fibermask, plates=plates, $
    legacy=legacy, _EXTRA=KeywordsForPhoto

   ; The correction vector is here --- adjust this as necessary.
   ; These are the same numbers as in SDSSFLUX2AB in the photoop product.
   correction = [-0.042, 0.036, 0.015, 0.013, -0.002]

   ;----------
   ; Trim to object fibers only, sort them, and trim to spectrographid

   plugmap = readPlateplugmap_sort(plugmap, hdr ,fibermask=fibermask, plates=plates)
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
      endif  ; ct gt 0
   endif  ; platelist_dir set

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


   ;----------
   ; Optionally add tags for SOS

   if (keyword_set(apotags)) then begin
      addtags = { cartid   : strtrim((yanny_par(hdr, 'cartridgeId'))[0],2), $
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
      addtags = replicate(create_struct('HEALPIX', 0L, $
                                        'HEALPIXGRP', 0, $
                                        'HEALPIX_DIR', ' '), n_elements(plugmap))
      plugmap = struct_addtags(plugmap, addtags)

      if keyword_set(healpix_now) then begin
            mwm_root='$MWM_HEALPIX'
            run2d=getenv('RUN2D')
            healpix_dir_t=plugmap.healpix_dir
            healp=coords_to_healpix(string(plugmap.ra), string(plugmap.dec))
            plugmap.healpix=healp.healpix
            plugmap.healpixgrp=healp.healpixgrp
            for fid = 0L, n_elements(plugmap)-1 do begin
                if not keyword_set(legacy) then begin
                    healpix_dir_t[fid]=mwm_root + $
                        strtrim(string(healp[fid].healpixgrp),2) + '/' + $
                        strtrim(string(healp[fid].healpix),2) + '/boss/' + $
                        strtrim(run2d, 2)+ '/' + 'spec-' + $
                        string(plateid,format='(i6.6)') + '-XXXX' + '-' + $
                        string(plugmap[fid].catalogid,format='(i11.11)')+'.fits'
                endif
            endfor
            plugmap.healpix_dir=healpix_dir_t
      endif
   endif
   ;----------
   if not keyword_set(apotags) then begin
        ; Read calibObj or photoPlate photometry data
        programname = (yanny_par(hdr, 'programname'))[0]
        ra_plate=float(yanny_par(hdr, 'raCen'))
        dec_plate=float(yanny_par(hdr, 'decCen'))
        plugmap=calibrobj(plugfile,plugmap, plateid, ra_plate,dec_plate, plates=plates,$
                            legacy=legacy, KeywordsForPhoto=KeywordsForPhoto)
   endif
   ;-----

   if tag_exist(plugmap,'ORG_FIBERID') || tag_exist(plugmap,'org_catid') then begin
     plugmap = struct_trimtags(plugmap,except_tags=['ORG_FIBERID','org_catid'])
   endif

   return, plugmap
end


;------------------------------------------------------------------------------
function prerun_readplugmap, plugfile, outfile, plugdir=plugdir, apotags=apotags, logfile=logfile, $
    exptime=exptime, hdr=hdr, fibermask=fibermask, plates=plates, legacy=legacy, $
    nfiles=nfiles, _EXTRA=KeywordsForPhoto
    
    if keyword_set(logfile) then begin
      cpbackup, logfile
      splog, filename=logfile
      splog, 'Log file ' + logfile + ' opened ' + systime()
   endif

   if n_elements(plugfile) GT 1 then begin
      for i=0, n_elements(plugfile)-1 do begin
        plugmap1 = readplugmap(plugfile[i], plugdir=plugdir, $
            apotags=apotags, deredden=deredden, exptime=exptime,  $
            hdr=hdr1, fibermask=fibermask1, plates=plates, legacy=legacy, $
            /nfiles, _EXTRA=KeywordsForPhoto)
        if keyword_set(plugmap) then plugmap = [plugmap ,plugmap1] else plugmap = plugmap1
        if keyword_set(hdr) then hdr = [hdr,'cut',hdr1] else hdr = ['cut',hdr1]
        if keyword_set(fibermask) then fibermask = [fibermask,-100, fibermask1] $
        else fibermask=[-100,fibermask1]
        Undefine, fibermask1
      endfor
      hdr = [hdr,'cut', '  ']
      fibermask=[fibermask,-100]
      
      if not keyword_set(apotags) then begin
        MWRFITS, junk, outfile, hdr, /create, /silent
            
        sxaddpar, plugmap_val, 'PLUGFILE', plugfile, ' plugfiles'
        sxaddpar, plugmap_val, 'EXTNAME', 'FIBERMAP', ' Complete Plugmap/confSummary'
        MWRFITS, plugmap, outfile, plugmap_val, Status=Status
        sxdelpar, plugmap_val, 'COMMENT'

        sxaddpar, fibermask_val, 'EXTNAME', 'FIBERMASK', ' Fibermask'
        MWRFITS, fibermask, outfile, fibermask_val ,Status=Status
        sxdelpar, fibermask_val, 'COMMENT'
        
      endif

      return, plugmap
   endif
    hdr = 0 ; Default return value
    if (keyword_set(fibermask)) then begin
        message, 'FIBERMASK is already set!'
    endif

    if keyword_set(plates) or keyword_set(legacy) then begin
        PMobjName='PLUGMAPOBJ'
        maptype='PlugMap'
        thisfile = (findfile(djs_filepath(plugfile, root_dir=plugdir), count=ct))[0]
    endif else begin
        PMobjName='FIBERMAP'
        maptype='confSummary'
        if keyword_set(plugdir) then begin
            thisfile = (findfile(djs_filepath(plugfile, root_dir=plugdir, subdir='*'), count=ct))[0]
        endif else begin
            thisfile = plugfile
            ct=1
        endelse
    endelse
    ;----------
    ; Read the file
    if (ct NE 1) then begin
        splog, 'WARNING: Cannot find '+PMobjName+' file ' + plugfile
        return, 0
    endif
    
            

    yanny_read, thisfile, pstruct, hdr=hdr, stnames=stnames, /anonymous
    if (NOT keyword_set(pstruct)) then begin
        splog, 'WARNING: Invalid '+PMobjName+' file ' + thisfile
        return, 0
    endif

    plugmap = *pstruct[(where(stnames EQ PMobjName))[0]]

    plugmap.ra = (360d0 + plugmap.ra) MOD 360d0


    if keyword_set(plates) or keyword_set(legacy) then begin
        plugmap = readPlateplugMap(plugfile,plugmap,stnames, mjd=mjd,$
                                   apotags=apotags, deredden=deredden, $
                                   exptime=exptime,  $
                                   hdr=hdr, fibermask=fibermask, $
                                   plates=plates,legacy=legacy, $
                                   _EXTRA=KeywordsForPhoto)
    endif else begin
        plugmap = readFPSobsSummary(plugfile,plugmap,stnames, mjd=mjd, $
                                   apotags=apotags,deredden=deredden, $
                                   exptime=exptime,  $
                                   hdr=hdr,fibermask=fibermask,  $
                                   _EXTRA=KeywordsForPhoto)
       if not keyword_set(nfiles) then begin
            hdr= ['cut',hdr,'cut', ' ']
            fibermask = [-100,fibermask,-100]
       endif
    endelse
    ; struct_print, plugmap, filename=repstr(plugfile,'.par','.html'), /html

    if (not keyword_set(apotags)) AND (not keyword_set(nfiles)) then  begin
        MWRFITS, junk, outfile, hdr, /create, /silent
    
        sxaddpar, plugmap_val, 'PLUGFILE', plugfile, ' plugfiles'
        sxaddpar, plugmap_val, 'EXTNAME', 'FIBERMAP', ' Complete Plugmap/confSummary'
        MWRFITS, plugmap, outfile, plugmap_val, Status=Status
        sxdelpar, plugmap_val, 'COMMENT'

        sxaddpar, fibermask_val, 'EXTNAME', 'FIBERMASK', ' Fibermask'
        MWRFITS, fibermask, outfile, fibermask_val ,Status=Status
        sxdelpar, fibermask_val, 'COMMENT'
        
    endif
    if keyword_set(logfile) then splog, /close
    return, plugmap
end
