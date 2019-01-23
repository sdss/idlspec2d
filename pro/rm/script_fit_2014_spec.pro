; this is a script to fit the 2014 epochs (32) with a detailed linelist


pro script_fit_2014_spec, range=range, silent=silent

file='/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
target=mrdfits(file,1)

plate=target[0].plate
mjd=target[0].mjd

plate = plate[0:31]
mjd = mjd[0:31]
nep = n_elements(mjd)

if ~keyword_set(range) then range = [0,nep-1]

plate = plate[range[0]:range[1]]
mjd = mjd[range[0]:range[1]]
nep = n_elements(mjd)


ntrial = 25
emparfile=getenv('IDLRM_DIR')+'/etc/qsoline_all.par'
;emparfile=getenv('IDLRM_DIR')+'/etc/qsoline.par'
linename=['SII6718','Halpha', 'Halpha_br', 'Hbeta', 'Hbeta_br', 'HeII4687','HeII4687_BR', 'OIII5007', $
                  'OIII5007c','MgII', 'MgII_br', 'CIII','CIII_br', 'SiIII1892','AlIII1857', 'NIII1750','CIV', 'CIV_br', $
                  'HeII1640','HeII1640_br', 'SiIV_OIV', 'OI1304', 'Lya', 'NV1240']

for i1=0, nep-1 do begin

  outdir = getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
       + '/' + string(plate[i1],format='(i4.4)') + '/wh_skysub/qsofit/'

  rm_fitplate,plate[i1],mjd[i1],ntrial=ntrial,emparfile=emparfile,outdir=outdir,$
       linename=linename,silent=silent

endfor


end
