;+
; NAME:
;   zfind
;
; PURPOSE:
;   Find possible redshift matches for a set of spectra using a set of
;   eigen-templates.
;
; CALLING SEQUENCE:
;   result = zfind( objflux, objivar, [ hdr=, objloglam0=, objdloglam=, $
;    zmin=, zmax= ])
;
; INPUTS:
;   objflux    - Flux for spectra [NPIX,NOBJECT]
;   objivar    - Inverse variance of flux [NPIX,NOBJECT]
;
; OPTIONAL INPUTS:
;   hdr        - FITS header for objects, used to construct the wavelengths
;                from the following keywords: COEFF0, COEFF1.
;                Must be specified if OBJLOGLAM0,OBJDLOGLAM are not set.
;   objloglam0 - Zero-pint of log-10(Angstrom) wavelength mapping ofOBJFLUX.
;   objdloglam - Wavelength spacing for OBJFLUX in log-10(Angstrom).
;   zmin       - Minimum redshift to consider; default to -0.002
;                (-600 km/sec).
;   zmax       - Minimum redshift to consider; default to +0.002
;                (+600 km/sec).
;
; OUTPUTS:
;   res        - Output structure with result for each object [NOBJECT].
;                The following elements are from the FITS header of the
;                best-fit Elodie spectrum:
;                  ELODIE_FILENAME  - Filename
;                  ELODIE_OBJECT    - Object name
;                  ELODIE_SPTYPE    - Spectral type
;                  ELODIE_BV        - (B-V) color
;                  ELODIE_TEFF      - T_effective
;                  ELODIE_LOGG      - Log10(gravity)
;                  ELODIE_FEH       - [Fe/H]
;                The following elements are from the ZFIND() function:
;                  ELODIE_Z         - Redshift
;                  ELODIE_Z_ERR     - Redshift error
;                  ELODIE_RCHI2     - Reduced chi^2
;                  ELODIE_DOF       - Degrees of freedom for fit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
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
    'elodie_filename' ,  '', $
    'elodie_object'   ,  '', $
    'elodie_sptype'   ,  '', $
    'elodie_bv'       , 0.0, $
    'elodie_teff'     , 0.0, $
    'elodie_logg'     , 0.0, $
    'elodie_feh'      , 0.0, $
    'elodie_z'        , 0.0, $
    'elodie_z_err'    , 0.0, $
    'elodie_rchi2'    , 0.0, $
    'elodie_dof'      ,  0L  $
   )

   return, result
end

;------------------------------------------------------------------------------
function elodie_best, objflux, objivar, $
 hdr=objhdr, objloglam0=objloglam0, objdloglam=objdloglam, $
 zmin=zmin, zmax=zmax

   if (n_params() LT 2) then begin
      print, 'Syntax - res = elodie_best(objflux, objivar, [ hdr=, $'
      print, ' objloglam0=, objdloglam=, zmin=, zmax= ]'
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
   if (NOT keyword_set(zmin)) then zmin = -0.002
   if (NOT keyword_set(zmax)) then zmax = 0.002

   stime0 = systime(1)

   ;----------
   ; Read all the Elodie spectra

   elodie_path = getenv('ELODIE_PATH')
   allfiles = findfile(filepath('00*', root_dir=elodie_path), count=nstar)
; ???
nstar=50
allfiles = allfiles[0:nstar-1]
   t0 = systime(1)
   starhdr = replicate(ptr_new(), nstar)
   for istar=0, nstar-1 do begin
      splog, 'Reading file ', istar+1, ' of ', nstar
      thisflux = read_elodie(allfiles[istar], loglam=starloglam, hdr=thishdr)
      if (NOT keyword_set(starflux)) then starflux = thisflux $
       else starflux = [[starflux],[thisflux]]
      starhdr[istar] = ptr_new(thishdr)
   endfor
   npix = n_elements(starloglam)
   splog, 'Time to read all files = ', systime(1) - t0

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

   ndim = size(objflux, /n_dimen)
   if (ndim EQ 1) then nobj = 1 $
    else nobj = (size(objflux, /dimens))[1]

   for istar=0, nstar-1 do begin
splog, 'Star number ', istar, ' of ', nstar
      res1 = zfind(objflux, objivar, hdr=objhdr, starflux=starflux[*,istar], $
       starloglam0=starloglam[0], npoly=3, zmin=zmin, zmax=zmax)
      if (istar EQ 0) then res_all = replicate(res1[0], nobj, nstar)
      res_all[*,istar] = res1[*]
   endfor

   ;----------
   ; For each object, select the best-fit star

   res_best = replicate(elodie_struct(), nobj)
   for iobj=0, nobj-1 do begin
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
   endfor

   splog, 'Total time for ELODIE_BEST = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'

; ???
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
