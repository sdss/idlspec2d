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
function readplugmap, plugfile, spectrographid, plugdir=plugdir, $
    apotags=apotags, deredden=deredden, exptime=exptime, calibobj=calibobj, $
    hdr=hdr, fibermask=fibermask, plates=plates, legacy=legacy, gaiaext=gaiaext, $
    MWM_fluxer=MWM_fluxer, nfiles=nfiles, _EXTRA=KeywordsForPhoto

    if not FILE_TEST(repstr(plugfile,'.par','.fits')) then begin
        splog, 'No fits fiber map exists: '+repstr(plugfile,'.par','.fits')
        splog, 'Creating New fits fiber map from: '+plugfile
        plugmap = prerun_readplugmap(plugfile, plugdir=plugdir, apotags=apotags, $
                                    exptime=exptime, hdr=hdr, fibermask=fibermask, $
                                    plates=plates, legacy=legacy, _EXTRA=KeywordsForPhoto)
    endif else begin
        splog, 'Reading fiber map from '+repstr(plugfile,'.par','.fits')
        ;rdfits_struct, repstr(plugfile,'.par','.fits'), plugmap, exten=[0]
        junk = mrdfits(repstr(plugfile,'.par','.fits'), 0, hdr)
        plugmap = mrdfits(repstr(plugfile,'.par','.fits'), 1)
        fibermask = mrdfits(repstr(plugfile,'.par','.fits'), 2)
    endelse

    fieldid = (yanny_par(hdr, 'field_id'))[0]
    ra_field=float(yanny_par(hdr, 'raCen'))
    dec_field=float(yanny_par(hdr, 'decCen'))

    if keyword_set(plates) and keyword_set(programname) then $
              stdtype = 'SPECTROPHOTO_STD' else stdtype = 'standard_boss'

    if keyword_set(calibobj) then begin
        rjce_extintion = 0
        if keyword_set(rjce_extintion) then begin
            if keyword_set(plates) and keyword_set(programname) then begin
                if ((strmatch(programname, '*MWM*', /fold_case) eq 1) $
                    || (strmatch(programname, '*OFFSET*', /fold_case) eq 1)) then begin
                    splog, "Using RJCE extintion"
                    spht = strmatch(fibermap.objtype, stdtype, /FOLD_CASE)
                    ispht = where(spht, nspht)
                    for i=0, n_elements(plugmap[ispht])-1 do begin
                        dat=(plugmap[ispht])[i].sfd_ebv_rjce
                        if (finite(dat) ne 0) and (dat gt 0) and (dat le 1.2*(plugmap[ispht])[i].sfd_ebv) then $
                            (plugmap[ispht])[i].sfd_ebv=(plugmap[ispht])[i].sfd_ebv_rjce
                        ;plugmap[ispht].sfd_ebv=plugmap[ispht].sfd_ebv_rjce
                    endfor
                    gaiaext = 0
                endif
            endif
        endif

        if keyword_set(MWM_fluxer) then begin
            if keyword_set(plates) then begin
                if keyword_set(programname) then begin
                    if ((strmatch(programname, '*MWM*', /fold_case) eq 1) $
                        || (strmatch(programname, '*OFFSET*', /fold_case) eq 1)) then begin
                        gaiaext = 1
                    endif
                endif
            endif else begin
                euler, ra_field, dec_field, ll_field, bb_field, 1
                if abs(bb_field) lt 15 then gaiaext = 1
            endelse
        endif

        if keyword_set(gaiaext) then begin
            splog, "Using dust_3d_map"
            spht = strmatch(fibermap.objtype, stdtype, /FOLD_CASE)
            ispht = where(spht, nspht)
            for i=0, n_elements(plugmap[ispht])-1 do begin
                dat=(plugmap[ispht])[i].sfd_ebv_gaia
                if (finite(dat) ne 0) and (dat le (plugmap[ispht])[i].sfd_ebv) then $
                    (plugmap[ispht])[i].sfd_ebv=(plugmap[ispht])[i].sfd_ebv_gaia
            endfor
        endif
    endif

    if keyword_set(plates) or keyword_set(legacy) then begin
        if keyword_set(deredden) then begin
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
        endif
        if keyword_set(plates) then begin
            nfiber = 1000
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
    endif else begin
        if (keyword_set(spectrographid)) then begin
            indx = where(plugmap.spectrographid eq spectrographid)
            plugmap = plugmap[indx]
            fibermask = fibermask[indx]
        endif
    endelse
    
    
    
    return, plugmap
end
