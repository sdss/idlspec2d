; coadd the spectra taken during each year by calling rm_coaddspec
; this does not incorporate the fix of the mismatched fibers

pro script_coaddspec

target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
fibermap = mrdfits(target_file,1,/silent)

mjdall = fibermap[0].mjd

; 2016 coadd
;ind = where(mjdall ge 57428 and mjdall le 57576)
;outfile = '/data3/yshen/spectro/bossredux/v5_7_1/0000/wh_skysub/spPlate-0000-57428-57576.fits'
;rm_coaddspec, mjd_coadd=mjdall[ind], outfile=outfile

; 2017 coadd
;ind = where(mjdall ge 57781 and mjdall le 57934)
;outfile = '/data3/yshen/spectro/bossredux/v5_7_1/0000/wh_skysub/spPlate-0000-57781-57934.fits'
;rm_coaddspec, mjd_coadd=mjdall[ind], outfile=outfile


end

