; output a fits file of all 1214 known quasars in the SDSS-RM field
; only 849 of them were tiled

pro output_all_qso_rmfield

file = '/home/yshen/products/Linux/idlrm/etc/targets_final.fits'
columns = ['RA', 'DEC', 'Z', 'PLATE', 'FIBERID', 'MJD', 'PSFMAG', 'OBJC_TYPE', 'RELEASE']
tt = mrdfits(file,1, columns=columns)

outfile = '/data3/yshen/work/sdssrm_sample_char/ancil_data/allqso_sdssrm.fits'

mwrfits, tt, outfile, /create

end
