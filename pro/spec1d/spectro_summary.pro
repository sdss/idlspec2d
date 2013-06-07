;+
;
; NAME:
;  spectro_summary
;
; PURPOSE:
;  Summarize spectro numbers for BOSS
;
; USAGE:
;  struc = spectro_summary(topdir=topdir, run2d=run2d, run1d=run1d)
;
; ARGUMENTS:
;  topdir: optional override for $BOSS_SPECTRO_REDUX
;  run2d: optional override for $RUN2D
;  run1d: optional override for $RUN1D
;
; RETURNS:
;  Structure with the following tags:
;
;   RUN2D          : run2d string used
;   RUN1D          : run1d string used
;   N_PLUGGING     : number of plate pluggings
;   N_PLATE        : number of unique plates
;   N_TILE         : number of unique tiles
;   N_SPECTRA      : number of spectra
;   N_SPECTRA_EFF  : effetive number of spectra
;   N_SPECTRA_UNIQ : number of unique spectro lines of sight
;   N_CMASS        : number of CMASS spectra
;   N_CMASS_UNIQ   : number of unique CMASS spectra
;   N_CMASS_ZGOOD  : number of previous with confident class & z
;   N_CMASS_ISGAL  : number of previous that are galaxies
;   N_LOZ          : number of LOZ spectra
;   N_LOZ_UNIQ     : number of unique LOZ spectra
;   N_LOZ_ZGOOD    : number of previous with confident class & z
;   N_LOZ_ISGAL    : number of previous that are galaxies
;   N_CMLOZ        : number of spectra that are both CMASS and LOZ
;   N_CMLOZ_UNIQ   : number of unique CMASS && LOZ spectra
;   N_CMLOZ_ZGOOD  : number of previous with confident class & z
;   N_CMLOZ_ISGAL  : number of previous that are galaxies
;   N_QSO          : number of QSO target sample spectra
;   N_QSO_UNIQ     : number of unique QSO target spectra
;   N_QSO_ZGOOD    : number of previous with confident class & z
;   N_QSO_ISQSO    : number of previous that are QSOs (spectro class)
;   N_QSO_LYA      : number of previous that have 2.15 < z < 3.5
;   N_ANCIL        : number of ancillary target spectra
;   N_ANCIL_UNIQ   : number of unique ancillary targets
;   N_SKY          : number of sky spectra
;   N_SKY_UNIQ     : number of unique sky-spectrum lines of sight
;   N_STD          : number of spectrophoto standard spectra
;   N_STD_UNIQ     : number of unique spectrophoto standards
;   N_OTHER        : number of other spectra
;   N_OTHER_UNIQ   : number of unique other objects
;   N_UNKNOWN      : number of spectra with uncertain classification
;   N_UNKNOWN_UNIQ : number of unique spectra with uncertain classification
;        
;  Removed for DR10:
;   N_QSO_SCANNED  : number of FPG-scanned unique QSO sample targets
;   N_LYA_INCOMP   : confident 2.1 < z < 3.0 FPG QSOs missed by pipeline
;   N_LYA_IMPURE   : pipeline 2.1 < z < 3.0 QSOs w/ conflicting FPG info
;
;
; WRITTEN BY:
;  bolton@utah 2012apr--jun
;        
;  UPDATES:
;     Decoupled FPG QSO scans from spAll for DR10;
;     removed incompleteness and impurity statistics here.
;     Could be re-added based upon final DRxQ.fits input when available.
;        -- Stephen Bailey, Oct 2012
;
;     Also exclude BAD_TARGETs from effective number of spectra
;        -- bolton@utah 2012dec
;
;-

function spectro_summary, topdir=topdir, run2d=run2d, run1d=run1d, printtable=printtable

if (not keyword_set(topdir)) then topdir = getenv('BOSS_SPECTRO_REDUX')
if (not keyword_set(run2d)) then run2d = getenv('RUN2D')
if (not keyword_set(run1d)) then run1d = getenv('RUN1D')

cols = [ $
       'OBJTYPE', $
       'TILE', $
       'PLATE', $
       'MJD', $
       'FIBERID', $
       'PLUG_RA', $
       'PLUG_DEC', $
       'THING_ID', $
       'PLATEQUALITY', $
       'Z', $
       'Z_ERR', $
       'ZWARNING', $
       'CLASS', $
       'Z_NOQSO', $
       'Z_ERR_NOQSO', $
       'ZWARNING_NOQSO', $
       'CLASS_NOQSO', $
       'RCHI2DIFF_NOQSO', $
       'BOSS_TARGET1', $
       'ANCILLARY_TARGET1', $
       'ANCILLARY_TARGET2', $
       'FIBER2FLUX', $
       'FIBER2MAG', $
       'SPECPRIMARY', $
       'SN_MEDIAN' ]
       
spf = topdir + '/' + run2d + '/spAll-' + run1d + '.fits'
print, 'Reading spAll file:'
print, spf
spall = hogg_mrdfits(spf,1,columns=cols,/silent)

ostruc = {run2d: run2d, $
          run1d: run1d, $
          n_plugging: 0L, $
          n_plate: 0L, $
          n_tile: 0L, $
          n_spectra: 0L, $
          n_spectra_eff: 0L, $
          n_spectra_uniq: 0L, $
          n_gal: 0L, $
          n_gal_uniq: 0L, $
          n_cmass: 0L, $
          n_cmass_uniq: 0L, $
          n_cmass_zgood: 0L, $
          n_cmass_isgal: 0L, $
          n_cmass_totgal: 0L, $
          n_loz: 0L, $
          n_loz_uniq: 0L, $
          n_loz_zgood: 0L, $
          n_loz_isgal: 0L, $
          n_loz_totgal: 0L, $
          n_cmloz: 0L, $
          n_cmloz_uniq: 0L, $
          n_cmloz_zgood: 0L, $
          n_cmloz_isgal: 0L, $
          n_qso: 0L, $
          n_qso_uniq: 0L, $
          n_qso_zgood: 0L, $
          n_qso_isqso: 0L, $
          n_qso_lya: 0L, $
          n_qso_totisqso: 0L, $
          n_qso_totlya: 0L, $
          n_qso_allspec: 0L, $
          n_qso_alluniq: 0L, $
          ;; n_qso_scanned: 0L, $
          ;; n_lya_incomp: 0L, $
          ;; n_lya_impure: 0L, $
          n_ancil: 0L, $
          n_ancil_uniq: 0L, $
          n_star: 0L, $
          n_star_uniq: 0L, $
          n_sky: 0L, $
          n_sky_uniq: 0L, $
          n_std: 0L, $
          n_std_uniq: 0L, $
          n_unknown: 0L, $
          n_unknown_uniq: 0L, $
          n_other: 0L, $
          n_other_uniq: 0L}

print, 'Computing summary statistics...'

; Number of pluggings:
ostruc.n_plugging = n_elements(spall)/1000L

; Number of unique plates:
uplate = spall.plate
uplate = uplate[sort(uplate)]
uplate = uplate[uniq(uplate)]
ostruc.n_plate = n_elements(uplate)

; Number of tiles:
utile = spall.tile
utile = utile[sort(utile)]
utile = utile[uniq(utile)]
ostruc.n_tile = n_elements(utile)

; Number of spectra:
ostruc.n_spectra = n_elements(spall)

; Flag out unplugged, little coverage, and bad target:
with_data = (spall.zwarning and (2^1 + 2^7 + 2^8)) eq 0
ostruc.n_spectra_eff = total(with_data)

; Unique spectra:
is_uniq = (spall.specprimary gt 0)
ostruc.n_spectra_uniq = total(is_uniq * with_data)

; Get the CMASS galaxy subsest:
;;;;is_cmass = ((spall.boss_target1 AND 2L^1+2L^2+2L^3+2L^7) NE 0) $
is_cmass = ((spall.boss_target1 AND sdss_flagval('BOSS_TARGET1', 'GAL_CMASS')) NE 0) $
 and (spall.fiber2mag[3] lt 21.5) and with_data

; How many CMASS and uniq CMASS targets?
ostruc.n_cmass = total(is_cmass)
ostruc.n_cmass_uniq = total(is_cmass * is_uniq)

; How many CMASS with good redshifts?
ostruc.n_cmass_zgood = total(is_cmass * is_uniq * (spall.zwarning_noqso eq 0))

; How many are known to be galaxies?
ostruc.n_cmass_isgal = total(is_cmass * is_uniq * (spall.zwarning_noqso eq 0) $
                             * (strtrim(spall.class_noqso, 2) eq 'GALAXY'))

; Total CMASS spectra known to be galaxies
ostruc.n_cmass_totgal = total(is_cmass * (spall.zwarning_noqso eq 0) $
                              * (strtrim(spall.class_noqso, 2) eq 'GALAXY'))

; Get the LOZ galaxy subset:
is_loz = ((spall.boss_target1 AND sdss_flagval('BOSS_TARGET1', 'GAL_LOZ')) NE 0) $
 and (spall.fiber2mag[3] lt 21.5) and with_data

; How many LOZ and uniq LOZ targets?
ostruc.n_loz = total(is_loz)
ostruc.n_loz_uniq = total(is_loz * is_uniq)

; How many LOZ with good redshifts?
ostruc.n_loz_zgood = total(is_loz * is_uniq * (spall.zwarning_noqso eq 0))

; How many are known to be galaxies?
ostruc.n_loz_isgal = total(is_loz * is_uniq * (spall.zwarning_noqso eq 0) $
                           * (strtrim(spall.class_noqso, 2) eq 'GALAXY'))

; Total LOZ spectra known to be galaxies
ostruc.n_loz_totgal = total(is_loz * (spall.zwarning_noqso eq 0) $
                           * (strtrim(spall.class_noqso, 2) eq 'GALAXY'))

; How many are both CMASS and LOZ?
ostruc.n_cmloz = total(is_loz * is_cmass)
ostruc.n_cmloz_uniq = total(is_loz * is_cmass * is_uniq)

; How many CMASS && LOZ with good redshifts?
ostruc.n_cmloz_zgood = total(is_loz * is_cmass * (spall.specprimary gt 0.) * (spall.zwarning_noqso eq 0))

; How many are known to be galaxies?
ostruc.n_cmloz_isgal = total(is_loz * is_cmass * is_uniq * (spall.zwarning_noqso eq 0) $
                             * (strtrim(spall.class_noqso, 2) eq 'GALAXY'))

; All galaxies, including non-CMASS and non-LOZ galaxies
;   If LOZ or CMASS, use class_noqso and zwarning_noqso
;   Otherwise, use class and zwarning
is_gal_cm_or_loz = (is_cmass OR is_loz) AND (strtrim(spall.class_noqso, 2) eq 'GALAXY') AND (spall.zwarning_noqso eq 0)
is_gal_not_cm_or_loz = ((is_cmass eq 0) AND (is_loz eq 0) AND (strtrim(spall.class, 2) eq 'GALAXY') AND (spall.zwarning eq 0))
ostruc.n_gal = total(is_gal_cm_or_loz) + total(is_gal_not_cm_or_loz)
ostruc.n_gal_uniq = total(is_gal_cm_or_loz * is_uniq) + total(is_gal_not_cm_or_loz * is_uniq)

; Get the Main QSO sample:
qso_targbits = 0LL
qso_targbits += sdss_flagval('BOSS_TARGET1', 'QSO_KNOWN_MIDZ')
qso_targbits += sdss_flagval('BOSS_TARGET1', 'QSO_FIRST_BOSS')
qso_targbits += sdss_flagval('BOSS_TARGET1', 'QSO_CORE_MAIN')
qso_targbits += sdss_flagval('BOSS_TARGET1', 'QSO_BONUS_MAIN')

is_qsotarg = ((spall.boss_target1 AND qso_targbits) NE 0) * with_data

; Classified as QSO and not targeted as either CMASS or LOZ
is_qso = (spall.zwarning eq 0) * (strtrim(spall.class, 2) eq 'QSO') $
         * (is_loz eq 0) * (is_cmass eq 0)

ostruc.n_qso = total(is_qsotarg)
ostruc.n_qso_uniq = total(is_qsotarg * is_uniq)
ostruc.n_qso_zgood = total(is_qsotarg * is_uniq * (spall.zwarning eq 0))
ostruc.n_qso_isqso = total(is_qsotarg * is_uniq * is_qso)
ostruc.n_qso_lya = total(is_qsotarg * is_uniq * is_qso $
                         * (spall.z ge 2.15) * (spall.z le 3.5))

;; Same thing without specprimary uniqueness cut
ostruc.n_qso_totisqso = total(is_qsotarg * is_qso)
ostruc.n_qso_totlya = total(is_qsotarg * is_qso $
                       * (spall.z ge 2.15) * (spall.z le 3.5))

;; QSOs regardless of target type *except* excluding LOZ and CMASS
;; --> DR9 paper included LOZ and CMASS targets but shouldn't have
ostruc.n_qso_allspec = total(is_qso)
ostruc.n_qso_alluniq = total(is_qso * is_uniq)

; Ancillary targets:
is_ancil = ((spall.ancillary_target1 gt 0) or (spall.ancillary_target2 gt 0)) * with_data
ostruc.n_ancil = total(is_ancil)
ostruc.n_ancil_uniq = total(is_ancil * (spall.specprimary gt 0))

; Stars
is_star = (strtrim(spall.class, 2) eq 'STAR') AND (spall.zwarning eq 0) * with_data
ostruc.n_star = total(is_star)
ostruc.n_star_uniq = total(is_star * (spall.specprimary gt 0))

; Spectrophotometric standards:
is_std = strmatch(spall.objtype, '*SPECTROPHOTO_STD*') * with_data
ostruc.n_std = total(is_std)
ostruc.n_std_uniq = total(is_std * (spall.specprimary gt 0))

; Sky fibers:
is_sky = strmatch(spall.objtype, '*SKY*') * with_data
ostruc.n_sky = total(is_sky)
ostruc.n_sky_uniq = total(is_sky * (spall.specprimary gt 0))

; Everything else:
is_other = (is_cmass eq 0) * (is_loz eq 0) * (is_qsotarg eq 0) * (is_ancil eq 0) * $
           (is_sky eq 0) * (is_std eq 0) * with_data
ostruc.n_other = total(is_other)
ostruc.n_other_uniq = total(is_other * (spall.specprimary gt 0))

; Unknown
is_unknown = (((is_cmass OR is_loz) AND (spall.zwarning_noqso gt 0)) OR $
             (spall.zwarning gt 0)) AND with_data
is_unknown *= is_unknown and (is_sky eq 0)
ostruc.n_unknown = total(is_unknown)
ostruc.n_unknown_uniq = total(is_unknown * is_uniq)

if keyword_set(printtable) then begin
   print, "Plates                 ", ostruc.n_plugging, ostruc.n_plate
   print, "Spectra                ", ostruc.n_spectra_eff, ostruc.n_spectra_uniq
   print, "All Galaxies           ", ostruc.n_gal, ostruc.n_gal_uniq
   print, "CMASS Galaxies         ", ostruc.n_cmass_totgal, ostruc.n_cmass_isgal
   print, "LOWZ Galaxies          ", ostruc.n_loz_totgal, ostruc.n_loz_isgal
   print, "All Quasars            ", ostruc.n_qso_allspec, ostruc.n_qso_alluniq
   print, "Main Quasars           ", ostruc.n_qso_totisqso, ostruc.n_qso_isqso
   print, "  with 2.15 <= z <= 3.5", ostruc.n_qso_totlya, ostruc.n_qso_lya
   print, "Ancillary spectra      ", ostruc.n_ancil, ostruc.n_ancil_uniq
   print, "Stars                  ", ostruc.n_star, ostruc.n_star_uniq
   print, "Standard Stars         ", ostruc.n_std, ostruc.n_std_uniq
   print, "Sky                    ", ostruc.n_sky, ostruc.n_sky_uniq
   print, "Unknown                ", ostruc.n_unknown, ostruc.n_unknown_uniq
end

;;-------------------------------------------------------------------
;; START FPG QSO SCAN COMPARISON
;;  ; spAll subset for objects with FPG inspections:
;;  wh_fpg = where(spall.z_conf_person gt 0)
;;  spall_fpg = spall[wh_fpg]
;;  
;;  ; spAll subset for QSO target specprimary:
;;  wh_qso = where(is_qsotarg and (spall.specprimary gt 0))
;;  spall_qso = spall[wh_qso]
;;  
;;  ;help, where(spall_qso.z_conf_person eq 0) 
;;  ;<Expression>    LONG      = Array[299]
;;  
;;  ; See if all the FPG stuff matches internally:
;;  ;spherematch, spall_fpg.plug_ra, spall_fpg.plug_dec, spall_fpg.plug_ra, spall_fpg.plug_dec, $
;;  ;             1.0/3600., m1, m2, d12, maxmatch=0
;;  ;wh_cross = where(m1 ne m2)
;;  ;m1 = m1[wh_cross]
;;  ;m2 = m2[wh_cross]
;;  ;splot, spall_fpg[m1].z_person, spall_fpg[m2].z_person, ps=1, color=2
;;  ;print, minmax(spall_fpg[m1].z_person - spall_fpg[m2].z_person)
;;  ;print, minmax(spall_fpg[m1].z_conf_person - spall_fpg[m2].z_conf_person)
;;  ;print, minmax(spall_fpg[m1].class_person - spall_fpg[m2].class_person)
;;  ; There are a couple strange mis-matches, but all seems to be in
;;  ; reasonably good agreements...
;;  
;;  spherematch, spall_qso.plug_ra, spall_qso.plug_dec, spall_fpg.plug_ra, spall_fpg.plug_dec, $
;;              1./3600., m_qso, m_fpg, d12, maxmatch=0
;;  
;;  ; Fill in:
;;  spall_qso[m_qso].z_person = spall_fpg[m_fpg].z_person
;;  spall_qso[m_qso].z_conf_person = spall_fpg[m_fpg].z_conf_person
;;  spall_qso[m_qso].class_person = spall_fpg[m_fpg].class_person
;;  
;;  ; Comparison computations:
;;  dzmax = 0.05
;;  pipe_zgood = (spall_qso.zwarning eq 0) and (strtrim(spall_qso.class, 2) eq 'QSO')
;;  fpg_zgood = (spall_qso.z_conf_person ge 3) and (spall_qso.class_person eq 3)
;;  pipe_zconf = (spall_qso.zwarning eq 0)
;;  fpg_zconf = (spall_qso.z_conf_person ge 3)
;;  in_agreement = abs(spall_qso.z - spall_qso.z_person) le dzmax
;;  is_scanned = spall_qso.z_conf_person gt 0
;;  in_range_pipe = (spall_qso.z ge 2.1) and (spall_qso.z le 3.0)
;;  in_range_fpg = (spall_qso.z_person ge 2.1) and (spall_qso.z_person le 3.0)
;;  
;;  
;;  
;;  ; Various other trial comparisons, none of which appear to be super-informative:
;;  ;whx = where(is_scanned * pipe_zgood * (spall_qso.z_conf_person ge 3) * $
;;  ;            (spall_qso.class_person eq 4) * (spall_qso.z ge 2.0))
;;  
;;  ;whx = where(is_scanned * pipe_zgood * (spall_qso.z_conf_person ge 3) * (spall_qso.z ge 2.0) * $
;;  ;            (spall_qso.class_person eq 1))
;;  
;;  
;;  ; N.B.: If it's a GALAXY, FPG doesn't usually provide a redshift.
;;  
;;  ; Number scanned by FPG:
;;  ostruc.n_qso_scanned = total(is_scanned)
;;  
;;  ; "Impurity": both confident, but in disagreement:
;;  is_impure = (in_agreement eq 0) * pipe_zgood * in_range_pipe * fpg_zconf
;;  
;;  ;wh_junk = where(is_impure)
;;  ;uuplotspec, spall_qso[wh_junk].plate, spall_qso[wh_junk].fiberid, mjd=spall_qso[wh_junk].mjd
;;  ; A few absorption issues, a few line mis-IDs.
;;  
;;  ; "Incompleteness": FPG confident, pipeline not:
;;  is_incomplete = fpg_zgood * (pipe_zgood eq 0) * in_range_fpg
;;  ;wh_junk = where(is_incomplete)
;;  ;uuplotspec, spall_qso[wh_junk].plate, spall_qso[wh_junk].fiberid, mjd=spall_qso[wh_junk].mjd
;;  ; Low SNR, negative emission, small-delta-chi2 between two reasonably
;;  ; close z's.  Those are the main explanations.
;;  
;;  ostruc.n_lya_incomp = total(is_incomplete)
;;  ostruc.n_lya_impure = total(is_impure)
;;  
;; END FPG QSO SCAN COMPARISON
;;-------------------------------------------------------------------

return, ostruc

end
