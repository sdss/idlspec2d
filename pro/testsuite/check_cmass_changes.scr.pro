;
; check_cmass_changes.scr.pro
;
; script to check changed CMASS redshifts.
;
; bolton@utah nov 2012
;

if (not keyword_set(v_old)) then v_old = 'v5_4_45'
if (not keyword_set(v_new)) then v_new = 'v5_5_12'

; The two spAll files:
f_old = getenv('BOSS_SPECTRO_REDUX') + '/' $
        + v_old + '/spAll-' + v_old + '.fits'
f_new = getenv('BOSS_SPECTRO_REDUX') + '/' $
        + v_new + '/spAll-' + v_new + '.fits'

; columns:
cols = ['plate', 'mjd', 'fiberid', $
        'plug_ra', 'plug_dec', 'boss_target1', $
        'run2d', 'run1d', 'class', 'z', 'zwarning', $
        'class_noqso', 'z_noqso', 'zwarning_noqso']

sp_old = hogg_mrdfits(f_old, 1, columns=cols)
sp_new = hogg_mrdfits(f_new, 1, columns=cols)

match_spall, sp_old, sp_new, m_old, m_new

print, minmax(sp_new[m_new].plate - sp_old[m_old].plate)
;           0           0
print, minmax(sp_new[m_new].mjd - sp_old[m_old].mjd)
;           0           0
print, minmax(sp_new[m_new].fiberid - sp_old[m_old].fiberid)
;           0           0

sp_old = sp_old[m_old]
sp_new = sp_new[m_new]

cmass_val = sdss_flagval('boss_target1', 'gal_cmass')

dz = 0.001
is_cmass = ((sp_old.boss_target1 and cmass_val) gt 0) and $
           ((sp_new.boss_target1 and cmass_val) gt 0)
is_goodboth = (sp_old.zwarning_noqso eq 0) and $
              (sp_new.zwarning_noqso eq 0)
is_zdiff = (abs(sp_new.z_noqso - sp_old.z_noqso) gt dz)

print, total(is_cmass * is_goodboth)
;      373490.
print, total(is_cmass * is_goodboth * is_zdiff) / total(is_cmass * is_goodboth)
;  0.000198131

