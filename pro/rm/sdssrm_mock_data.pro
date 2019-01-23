; make mock LC data for the sdss-rm sample, e.g., same cadence and L-z range,
; but random LCs from DRW models

pro sdssrm_mock_data


; get the mock sample
infile = '/data2/yshen/work/RM/TDSS/level2/sim_qso_w_LC_tf_width_tenth.fits'
sim_LC_all = mrdfits(infile, 1)  ; all 100,000 objects

; get the actual numbers in each z-Mi grid
z_grid = 0.1 + findgen(25)*0.2
imag_grid = 15.5 + findgen(7)*1.0
file = '/data2/yshen/work/RM/TDSS/level2/true_nqso_HRH_LF.fits'
ttt = mrdfits(file, 1, /silent)
true_nqso = ttt.true_nqso ; can be decimal numbers


; downgrade to a single plate
subind = downsample_mock_qso(nqso = true_nqso, sim_LC_all = sim_LC_all)
sim_LC = sim_LC_all[subind]
; keep only i<21.7
ind = where(sim_LC.imag lt 21.7)
sim_LC = sim_LC[ind]
ntot = n_elements(sim_LC)


; get the actual cadence from SDSS-RM
file = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
target=mrdfits(file,1)
mjd = target[0].mjd
mjd=mjd[0:31]
imjd=mjd - min(mjd)
nmjd=n_elements(mjd)

remove_tags, sim_LC, ['CONTI_LC','LINE_LC'], result
tmp = {conti_lc:dblarr(nmjd), line_lc:dblarr(nmjd), mjd:mjd}
tmp = replicate(tmp, ntot)
conti_lc = sim_LC.conti_lc & line_lc = sim_LC.line_lc

result = struct_addtags(result, tmp)
result.conti_lc = conti_lc[imjd,*]
result.line_lc = line_lc[imjd,*]

outfile = '/data3/yshen/work/composite_lag/mock_sdssrm_lc.fits'

mwrfits, result, outfile, /create

end
