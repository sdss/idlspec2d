;+
; NAME:
;   ztweak_star
;
; PURPOSE:
;   Find the best-fit Elodie spectrum to a set of spectra.
;
; CALLING SEQUENCE:
;   ztweak_star, [ filename, zmin=, zmax=, /overwrite ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   filename   - Yanny parameter file with at least the entries PLATE,
;                MJD, FIBERID, and optionally CZ.
;                Default to "$IDLSPEC2D_DIR/templates/eigeninput_star.par".
;   zmin       - Minimum redshift to consider; default to -0.00333
;                (-1000 km/sec).
;   zmax       - Minimum redshift to consider; default to +0.00333
;                (+1000 km/sec).
;   overwrite  - If set, then overwrite the input file with CZ replaced
;                with the best-fit value.
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
;   No attempt is made to precisely match the instrumental dispersion
;   of the SDSS spectra and the Elodie spectra.  The Elodie spectra are
;   smoothed to an instrumental dispersion of 70 km/sec.
;
; PROCEDURES CALLED:
;   elodie_best()
;   readspec
;   struct_addtags()
;   yanny_free
;   yanny_read
;
; REVISION HISTORY:
;   03-Apr-2002  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro ztweak_star, filename, zmin=zmin, zmax=zmax, overwrite=overwrite

   snmax = 100
   if (NOT keyword_set(zmin)) then zmin = -0.00333
   if (NOT keyword_set(zmax)) then zmax = 0.00333
   cspeed = 2.99792458d5

   ;----------
   ; Read the input spectra

   if (NOT keyword_set(filename)) then $
    filename = filepath('eigeninput_star.par', $
     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='templates')
   yanny_read, filename, pdat, hdr=hdr, enums=enums, stnames=stnames
   slist = *pdat[0]
   yanny_free, pdat
   nobj = n_elements(slist)

   ;----------
   ; If CZ doesn't exist in the structure, then add it as all zeros.

   if ((where(tag_names(slist) EQ 'CZ'))[0] EQ -1) then $
    slist = struct_addtags(slist, replicate( {cz: 0.0}, nobj) )

; ???
;ii=where(strtrim(slist.class,2) EQ 'K')
;slist=slist[ii]
;slist=slist[0:9]
;slist.plate = 406
;slist.mjd = 51869
;slist.fiberid = [9,18,21,34,39,53,62,64,80,102]
;nobj = n_elements(slist)
   readspec, slist.plate, slist.fiberid, mjd=slist.mjd, $
    flux=objflux, invvar=objivar, $
    andmask=andmask, ormask=ormask, plugmap=plugmap, loglam=objloglam, /align

   ;----------
   ; Insist that all of the requested spectra exist

   imissing = where(plugmap.fiberid EQ 0, nmissing)
   if (nmissing GT 0) then begin
      for i=0, nmissing-1 do $
       print, 'Missing plate=', slist[imissing[i]].plate, $
        ' mjd=', slist[imissing[i]].mjd, $
        ' fiber=', slist[imissing[i]].fiberid
      message, string(nmissing) + ' missing object(s)'
   endif

   ;----------
   ; Do not fit where the spectrum may be dominated by sky-sub residuals.

   objivar = skymask(objivar, andmask, ormask)
andmask = 0 ; Free memory
ormask = 0 ; Free memory

   if (keyword_set(snmax)) then begin
      ifix = where(objflux^2 * objivar GT snmax^2)
      if (ifix[0] NE -1) then objivar[ifix] = (snmax/objflux[ifix])^2
   endif

   ;----------
   ; Find the best-fit Elodie star for each spectrum

   objdloglam = objloglam[1] - objloglam[0]
   res = elodie_best(objflux, objivar, $
    objloglam0=objloglam[0], objdloglam=objdloglam, zmin=zmin, zmax=zmax)

   ;----------
   ; Print the differences between the input velocities and best-fit ones

   cz_in = slist.cz
   cz_out = res.elodie_z * cspeed
   cz_diff = cz_out - cz_in

   splog, file='ztweak_star.log'
   splog, 'PLATE  MJD   FIBER CZ_IN   CZ_OUT  CZ_DIFF SPTYPE  '
   splog, '-----  ----- ----- ------- ------- ------- --------'
   for iobj=0, nobj-1 do $
    splog, slist[iobj].plate, slist[iobj].mjd, slist[iobj].fiberid, $
     cz_in[iobj], cz_out[iobj], cz_diff[iobj], res[iobj].elodie_sptype, $
     format='(i5,i7,i6,3f8.1," ",a8)'
   splog, /close

   ;----------
   ; Optionally overwrite the input file.
   ; Do not use the input STRUCTS from YANNY_READ, since we may have
   ; modified the structure.

   if (keyword_set(overwrite)) then begin
      slist.cz = cz_out
      yanny_write, filename, ptr_new(slist), hdr=hdr, enums=enums, $
       stnames=stnames
   endif

stop

   return
end
;------------------------------------------------------------------------------
