;+
; NAME:
;   elodie_best
;
; PURPOSE:
;   Find the best-fit Elodie spectrum to a set of spectra.
;
; CALLING SEQUENCE:
;   res = elodie_best(objflux, objivar, $
;    hdr=, objloglam0=, objdloglam=, zmin=, zmax=, fitindx= ])
;
; INPUTS:
;   objflux    - Flux for spectra [NPIX,NOBJECT]
;   objivar    - Inverse variance of flux [NPIX,NOBJECT]
;
; OPTIONAL INPUTS:
;   hdr        - FITS header for objects, used to construct the wavelengths
;                from the following keywords: COEFF0, COEFF1.
;                Must be specified if OBJLOGLAM0,OBJDLOGLAM are not set.
;   objloglam0 - Zero-pint of log-10(Angstrom) wavelength mapping of OBJFLUX.
;   objdloglam - Wavelength spacing for OBJFLUX in log-10(Angstrom).
;   zmin       - Minimum redshift to consider; default to -0.00333
;                (-1000 km/sec).
;   zmax       - Minimum redshift to consider; default to +0.00333
;                (+1000 km/sec).
;   fitindx    - If set, then only fit for the objects specified by
;                these indices.  Set to -1 to not fit any objects (but
;                return empty output structures); default to fitting all.
;
; OUTPUTS:
;   res        - Output structure with result for each object [NOBJECT].
;                The following elements are from the FITS header of the
;                best-fit Elodie spectrum:
;                  ELODIE_FILENAME   - Filename
;                  ELODIE_OBJECT     - Object name
;                  ELODIE_SPTYPE     - Spectral type
;                  ELODIE_BV         - (B-V) color
;                  ELODIE_TEFF       - T_effective
;                  ELODIE_LOGG       - Log10(gravity)
;                  ELODIE_FEH        - [Fe/H]
;                  ELODIE_Z_MODELERR - The standard deviation in redshift
;                                      amongst the 12 best-fit stars
;                The following elements are from the ZFIND() function:
;                  ELODIE_Z          - Redshift
;                  ELODIE_Z_ERR      - Redshift error
;                  ELODIE_RCHI2      - Reduced chi^2
;                  ELODIE_DOF        - Degrees of freedom for fit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Any stars labelled as possible binary or triple stars are ignored.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;   read_elodie()
;   splog
;   sxpar()
;   zfind()
;
; INTERNAL SUPPORT ROUTINES:
;   elodie_struct()
;
; REVISION HISTORY:
;   11-Mar-2002  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
function elodie_struct

   result = create_struct( $
    name = 'ZELODIE', $
    'elodie_filename'  , ' ', $
    'elodie_object'    , ' ', $
    'elodie_sptype'    , ' ', $
    'elodie_bv'        , 0.0, $
    'elodie_teff'      , 0.0, $
    'elodie_logg'      , 0.0, $
    'elodie_feh'       , 0.0, $
    'elodie_z'         , 0.0, $
    'elodie_z_err'     , 0.0, $
    'elodie_z_modelerr', 0.0, $
    'elodie_rchi2'     , 0.0, $
    'elodie_dof'       ,  0L  $
   )

   return, result
end

;------------------------------------------------------------------------------
function elodie_best, objflux, objivar, $
 hdr=objhdr, objloglam0=objloglam0, objdloglam=objdloglam, $
 zmin=zmin, zmax=zmax, fitindx=fitindx

   if (n_params() LT 2) then begin
      print, 'Syntax - res = elodie_best(objflux, objivar, [ hdr=, $'
      print, ' objloglam0=, objdloglam=, zmin=, zmax=, fitindx= ]'
      return, 0
   endif
   if (NOT keyword_set(objhdr)) then begin
      if (NOT keyword_set(objloglam0) OR NOT keyword_set(objdloglam)) then begin
         print, 'Either HDR or OBJLOGLAM0,OBJDLOGLAM must be set'
         return, 0
      endif
      objhdr = ['COEFF0  = ' + string(objloglam0), $
                'COEFF1  = ' + string(objdloglam) ]
   endif
   if (NOT keyword_set(zmin)) then zmin = -0.00333
   if (NOT keyword_set(zmax)) then zmax = 0.00333

   stime0 = systime(1)

   ;----------
   ; Create the output structure, and return empty data if FITINDX=[-1].

   ndim = size(objflux, /n_dimen)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = (size(objflux, /dimens))[1]

   if (n_elements(fitindx) GT 0) then thisindx = fitindx $
    else thisindx = lindgen(nobj)

   res_best = replicate(elodie_struct(), nobj)
   if (thisindx[0] EQ -1) then return, res_best

   ;----------
   ; Read all the Elodie spectra

   allfiles = elodie_filelist(minwave=minwave)
   nstar = n_elements(allfiles)

   t0 = systime(1)
   starhdr = replicate(ptr_new(), nstar)
   for istar=0L, nstar-1 do begin
      splog, 'Reading file ', istar+1, ' of ', nstar
      thisflux = read_elodie(allfiles[istar], loglam=starloglam, hdr=thishdr, minloglamclip=alog10(minwave))
      if (NOT keyword_set(starflux)) then starflux = thisflux $
       else starflux = [[starflux],[thisflux]]
      starhdr[istar] = ptr_new(thishdr)
   endfor
   npix = n_elements(starloglam)
   splog, 'Time to read all files = ', systime(1) - t0

   ;----------
   ; Trim to stars not labelled as binary or triple stars

   qkeep = bytarr(nstar) + 1B
   for istar=0L, nstar-1 do begin
      q_vr = sxpar(*starhdr[istar],'Q_VR')
      if (strmatch(q_vr,'*double*') $
       OR strmatch(q_vr,'*triple*') $
       OR strmatch(q_vr,'*binary*')) then begin
         qkeep[istar] = 0
         ptr_free, starhdr[istar]
      endif
   endfor
   ikeep = where(qkeep, nstar)
   if (nstar EQ 0) then $
    message, 'All stars seem to be labelled as binaries or triples!'
   starflux = starflux[*,ikeep]
   starhdr = starhdr[ikeep]

   ;----------
   ; Trim wavelengths to those covered by the majority of the objects

   fracgpix = total(starflux NE 0, 2) / nstar
   igood = where(fracgpix GT 0.95)
   i1 = igood[0]
   i2 = (reverse(igood))[0]
   starloglam = starloglam[i1:i2]
   starflux = starflux[i1:i2,*]

   ;----------
   ; Interpolate over bad data, of which there is very little

   starflux = djs_maskinterp(starflux, starflux EQ 0, /const, iaxis=0)

   ;----------
   ; Compute the best-fit between each object and each star

   for istar=0, nstar-1 do begin
      splog, 'Star number ', istar, ' of ', nstar

      if (n_elements(fitindx) GT 0) then $
       res1 = zfind(objflux[*,thisindx], objivar[*,thisindx], $
        hdr=objhdr, starflux=starflux[*,istar], $
        starloglam0=starloglam[0], npoly=3, zmin=zmin, zmax=zmax) $
      else $
       res1 = zfind(objflux, objivar, $
        hdr=objhdr, starflux=starflux[*,istar], $
        starloglam0=starloglam[0], npoly=3, zmin=zmin, zmax=zmax)

      if (istar EQ 0) then res_all = replicate(res1[0], nobj, nstar)
      res_all[thisindx,istar] = res1[*]
   endfor

   ;----------
   ; For each object, select the best-fit star

   for ifit=0, n_elements(thisindx)-1 do begin
      iobj = thisindx[ifit]
      junk = min(res_all[iobj,*].rchi2, imin)
      res_best[iobj].elodie_filename = sxpar((*starhdr[imin]), 'FILENAME')
      res_best[iobj].elodie_object = sxpar((*starhdr[imin]), 'OBJECT')
      res_best[iobj].elodie_sptype = sxpar((*starhdr[imin]), 'SPTYPE')
      res_best[iobj].elodie_bv = sxpar((*starhdr[imin]), 'B-V')
      res_best[iobj].elodie_teff = sxpar((*starhdr[imin]), 'TEFF')
      res_best[iobj].elodie_logg = sxpar((*starhdr[imin]), 'LOGG')
      ipar = (where(strmatch(*starhdr[imin], 'HIERARCH \[Fe/H\]*')))[0]
      if (ipar NE -1) then $
       res_best[iobj].elodie_feh = float( (strsplit((*starhdr[imin])[ipar], $
        "HIERARCH [Fe/H] = '", /extract))[0] )
      res_best[iobj].elodie_z = res_all[iobj,imin].z
      res_best[iobj].elodie_z_err = res_all[iobj,imin].z_err
      res_best[iobj].elodie_rchi2 = res_all[iobj,imin].rchi2
      res_best[iobj].elodie_dof = res_all[iobj,imin].dof
      isort = sort(res_all[iobj,*].rchi2)
      res_best[iobj].elodie_z_modelerr = $
       stddev(res_all[iobj,isort[0:11]].z, /double)
   endfor

   splog, 'Total time for ELODIE_BEST = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'

   ; Free memory
   for istar=0L, nstar-1 do $
    ptr_free, starhdr[istar]

; Test???
;plot,zans.z*3e5,res_best.elodie_z*3e5,ps=7
;for i=0,nstar-1 do oplot,zans.z*3e5,res_all[*,i].z*3e5,ps=3
;oplot,[-200,200],[-200,200]
;
;junk=poly_array(2172,3)
;synflux=fltarr(2172,nobj)
;for i=0,26 do synflux[*,i]=junk#res_best[i].theta[1:3]
;plot,djs_median(abs(synflux),1)/djs_median(objflux,1)

; Also return the best-fit Elodie spectrum???
; Use same procedure as in SYNTHSPEC for sub-pixel shifts???

   return, res_best
end
;------------------------------------------------------------------------------
