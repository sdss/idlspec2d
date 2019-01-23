; output the qsofits for the BOSS-XMM sample (Mentzel et al.)

pro output_xmm_objfits

outdir_in='/data1/quasar/yshen/DR12_QSO_fits/'

file1=outdir_in+'qso_prop-xmm_boss.fits'
qso1=mrdfits(file1,1)
file2=outdir_in+'qso_prop-xmm_boss_add.fits'
qso2=mrdfits(file2,1)


linename = ['Halpha_br', 'Halpha_na', 'NII6549', $
  'NII6585', 'SII6718','SII6732','Hbeta_br', 'Hbeta_na', $
  'OIII4959', 'OIII5007','MgII_br', 'MgII_na', 'CIII_br', 'SIIII1892', $
  'ALIII1857', 'CIV_br', 'HEII1640', 'OIII1663']


outfile1='qso_prop-xmm_boss_w_optFe.fits'
rm_compile_qsofit, qso1.plate, qso1.mjd, $
     input_fiber=qso1.fiberid,rm_plate=0, $
     outdir_in=outdir_in,outfile=outfile1, linename=linename


outfile2='qso_prop-xmm_boss_add_w_optFe.fits'
rm_compile_qsofit, qso2.plate, qso2.mjd, $
     input_fiber=qso2.fiberid,rm_plate=0, $
     outdir_in=outdir_in,outfile=outfile2, linename=linename


end
