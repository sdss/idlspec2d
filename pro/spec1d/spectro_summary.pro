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
;   N_QSO_SCANNED  : number of FPG-scanned unique QSO sample targets
;   N_QSO_INCOMP   : confident FPG QSOs w/o confident pipeline class & z
;   N_QSO_IMPURE   : FPG & pipeline confident QSOs, with different z's
;   N_SKY          : number of sky spectra
;   N_SKY_UNIQ     : number of unique sky-spectrum lines of sight
;   N_STD          : number of spectrophoto standard spectra
;   N_STD_UNIQ     : number of unique spectrophoto standards
;   N_OTHER        : number of other spectra
;   N_OTHER_UNIQ   : number of unique other objects
;
;
; WRITTEN BY:
;  bolton@utah 2012apr--jun
;
;-

function spectro_summary, topdir=topdir, run2d=run2d, run1d=run1d

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
       'FIBER2FLUX', $
       'FIBER2MAG', $
       'SPECPRIMARY', $
       'SN_MEDIAN', $
       'Z_PERSON', $
       'CLASS_PERSON', $
       'Z_CONF_PERSON']

ostruc = {run2d: run2d, $
          run1d: run1d, $
          n_plugging: 0L, $
          n_plate: 0L, $
          n_tile: 0L, $
          n_spectra: 0L, $
          n_spectra_uniq: 0L, $
          n_cmass: 0L, $
          n_cmass_uniq: 0L, $
          n_cmass_zgood: 0L, $
          n_cmass_isgal: 0L, $
          n_loz: 0L, $
          n_loz_uniq: 0L, $
          n_loz_zgood: 0L, $
          n_loz_isgal: 0L, $
          n_cmloz: 0L, $
          n_cmloz_uniq: 0L, $
          n_cmloz_zgood: 0L, $
          n_cmloz_isgal: 0L, $
          n_qso: 0L, $
          n_qso_uniq: 0L, $
          n_qso_zgood: 0L, $
          n_qso_isqso: 0L, $
          n_qso_scanned: 0L, $
          n_qso_incomp: 0L, $
          n_qso_impure: 0L, $
          n_sky: 0L, $
          n_sky_uniq: 0L, $
          n_std: 0L, $
          n_std_uniq: 0L, $
          n_other: 0L, $
          n_other_uniq: 0L}

spf = topdir + '/' + run2d + '/spAll-' + run1d + '.fits'
print, 'Reading spAll file:'
print, spf
spall = hogg_mrdfits(spf,1,columns=cols,/silent)

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
ostruc.n_spectra_uniq = total(spall.specprimary gt 0)

; Get the CMASS galaxy subsest:
;;;;is_cmass = ((spall.boss_target1 AND 2L^1+2L^2+2L^3+2L^7) NE 0) $
is_cmass = ((spall.boss_target1 AND 2L^1) NE 0) $
 and (spall.fiber2mag[3] lt 21.5)

; How many CMASS?
ostruc.n_cmass = total(is_cmass)
ostruc.n_cmass_uniq = total(is_cmass * (spall.specprimary gt 0))

; How many CMASS with good redshifts?
ostruc.n_cmass_zgood = total(is_cmass * (spall.specprimary gt 0.) * (spall.zwarning_noqso eq 0))

; How many are known to be galaxies?
ostruc.n_cmass_isgal = total(is_cmass * (spall.specprimary gt 0.) * (spall.zwarning_noqso eq 0) $
                             * (strtrim(spall.class_noqso, 2) eq 'GALAXY'))

; Get the LOZ galaxy subset:
is_loz = ((spall.boss_target1 AND 2L^0) NE 0) $
 and (spall.fiber2mag[3] lt 21.5)

; How many LOZ?
ostruc.n_loz = total(is_loz)
ostruc.n_loz_uniq = total(is_loz * (spall.specprimary gt 0))

; How many LOZ with good redshifts?
ostruc.n_loz_zgood = total(is_loz * (spall.specprimary gt 0.) * (spall.zwarning_noqso eq 0))

; How many are known to be galaxies?
ostruc.n_loz_isgal = total(is_loz * (spall.specprimary gt 0.) * (spall.zwarning_noqso eq 0) $
                           * (strtrim(spall.class_noqso, 2) eq 'GALAXY'))

; How many are both CMASS and LOZ?
ostruc.n_cmloz = total(is_loz * is_cmass)
ostruc.n_cmloz_uniq = total(is_loz * is_cmass * (spall.specprimary gt 0))

; How many CMASS && LOZ with good redshifts?
ostruc.n_cmloz_zgood = total(is_loz * is_cmass * (spall.specprimary gt 0.) * (spall.zwarning_noqso eq 0))

; How many are known to be galaxies?
ostruc.n_cmloz_isgal = total(is_loz * is_cmass * (spall.specprimary gt 0.) * (spall.zwarning_noqso eq 0) $
                             * (strtrim(spall.class_noqso, 2) eq 'GALAXY'))

; Get the QSO sample:
is_qsotarg = ((spall.boss_target1 AND 3298535930880LL) NE 0)

ostruc.n_qso = total(is_qsotarg)
ostruc.n_qso_uniq = total(is_qsotarg * (spall.specprimary gt 0))
ostruc.n_qso_zgood = total(is_qsotarg * (spall.specprimary gt 0) * (spall.zwarning eq 0))
ostruc.n_qso_isqso = total(is_qsotarg * (spall.specprimary gt 0) * (spall.zwarning eq 0) $
                           * (strtrim(spall.class, 2) eq 'QSO'))

; spAll subset for objects with FPG inspections:
wh_fpg = where(spall.z_conf_person gt 0)
spall_fpg = spall[wh_fpg]

; spAll subset for QSO target specprimary:
wh_qso = where(is_qsotarg and (spall.specprimary gt 0))
spall_qso = spall[wh_qso]

;help, where(spall_qso.z_conf_person eq 0) 
;<Expression>    LONG      = Array[299]

; See if all the FPG stuff matches internally:
;spherematch, spall_fpg.plug_ra, spall_fpg.plug_dec, spall_fpg.plug_ra, spall_fpg.plug_dec, $
;             1.0/3600., m1, m2, d12, maxmatch=0
;wh_cross = where(m1 ne m2)
;m1 = m1[wh_cross]
;m2 = m2[wh_cross]
;splot, spall_fpg[m1].z_person, spall_fpg[m2].z_person, ps=1, color=2
;print, minmax(spall_fpg[m1].z_person - spall_fpg[m2].z_person)
;print, minmax(spall_fpg[m1].z_conf_person - spall_fpg[m2].z_conf_person)
;print, minmax(spall_fpg[m1].class_person - spall_fpg[m2].class_person)
; There are a couple strange mis-matches, but all seems to be in
; reasonably good agreements...

spherematch, spall_qso.plug_ra, spall_qso.plug_dec, spall_fpg.plug_ra, spall_fpg.plug_dec, $
             1./3600., m_qso, m_fpg, d12, maxmatch=0

; Fill in:
spall_qso[m_qso].z_person = spall_fpg[m_fpg].z_person
spall_qso[m_qso].z_conf_person = spall_fpg[m_fpg].z_conf_person
spall_qso[m_qso].class_person = spall_fpg[m_fpg].class_person

; Comparison computations:
dzmax = 0.05
pipe_zgood = (spall_qso.zwarning eq 0) and (strtrim(spall_qso.class, 2) eq 'QSO')
fpg_zgood = (spall_qso.z_conf_person ge 3) and (spall_qso.class_person eq 3)
in_agreement = abs(spall_qso.z - spall_qso.z_person) le dzmax
is_scanned = spall_qso.z_conf_person gt 0

; N.B.: If it's a GALAXY, FPG doesn't bother with the redshift.

; Number scanned by FPG:
ostruc.n_qso_scanned = total(is_scanned)

; "Impurity": both confident, but in disagreement:
is_impure = is_scanned * (in_agreement eq 0) * pipe_zgood * fpg_zgood

; "Incompleteness": FPG confident, pipeline not:
is_incomplete = is_scanned * fpg_zgood * (pipe_zgood eq 0)

ostruc.n_qso_incomp = total(is_incomplete)
ostruc.n_qso_impure = total(is_impure)

; Spectrophotometric standards:
is_std = strmatch(spall.objtype, '*SPECTROPHOTO_STD*')
ostruc.n_std = total(is_std)
ostruc.n_std_uniq = total(is_std * (spall.specprimary gt 0))

; Sky fibers:
is_sky = strmatch(spall.objtype, '*SKY*')
ostruc.n_sky = total(is_sky)
ostruc.n_sky_uniq = total(is_sky * (spall.specprimary gt 0))

; Everything else:
is_other = (is_cmass eq 0) * (is_loz eq 0) * (is_qsotarg eq 0) * (is_sky eq 0) * (is_std eq 0)
ostruc.n_other = total(is_other)
ostruc.n_other_uniq = total(is_other * (spall.specprimary gt 0))

return, ostruc

end
