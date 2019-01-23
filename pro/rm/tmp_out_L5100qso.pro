; make a table of rmid and logL5100_qso (from Shen++15b) for Jennifer Li

pro tmp_out_L5100qso

file='/data3/yshen/work/agn_host/decomp_final.fits'
decomp=mrdfits(file,1)
nobj=n_elements(decomp)

output=replicate({rmid:0L, logL5100_qso:0.D, logL5100_qso_err:-1.D}, nobj)
output.rmid=decomp.fiber - 1L
output.logL5100_qso=decomp.logL5100_qso
output.logL5100_qso_err=decomp.logL5100_qso_err

outfile='/data3/yshen/work/composite_lag/rmid_logLqso.fits'
mwrfits, output, outfile,/create

end
