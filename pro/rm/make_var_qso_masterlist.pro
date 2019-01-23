; make a master fits file for all the SDSS quasars with photometric
; variability from PS1 or Stripe 82

pro make_var_qso_masterlist

outdir = '/data3/quasar/yshen/work/md_qso/'

struct={sample:'',ra:0.D,dec:0.D,plate:0L,fiber:0L,mjd:0L,z_pipe:0.D,z_vi:0.D}

; MD01pilotQSOs.csv
;file = outdir + 'MD01pilotQSOs.csv'
;readcol,file,format='l,l,l,d,d,d',plate,mjd,fiber,ra,dec,z,delimiter=','
;nnn=n_elements(plate)
;result=replicate(struct, nnn)
;result.sample='MD01'
;result.plate=plate & result.fiber=fiber & result.mjd=mjd
;result.ra=ra & result.dec=dec & result.z_pipe=z
;result_all = result
; updated by P. Green; 
file = outdir + '6369-56217V_viQSO0gd_uwold2.csv'
readcol,file,format='l,x,l,l,x,d,d,x,x,d', delimiter=',',plate,mjd,fiber,ra,dec,z
nnn=n_elements(plate)
result=replicate(struct, nnn)
result.sample='MD01'
result.plate=plate & result.fiber=fiber & result.mjd=mjd
result.ra=ra & result.dec=dec & result.z_pipe=z
result_all = result


; MD03pilotQSOs.csv
;file = outdir + 'MD03pilotQSOs.csv'
;readcol,file,format='l,l,l,d,d,d',plate,mjd,fiber,ra,dec,z,delimiter=','
;nnn=n_elements(plate)
;result=replicate(struct, nnn)
;result.sample='MD03'
;result.plate=plate & result.fiber=fiber & result.mjd=mjd
;result.ra=ra & result.dec=dec & result.z_pipe=z
;result_all = [result_all, result]
; updated by P. Green:
file = outdir + '6783-56284V_viQSO0gd_uwold.csv'
readcol,file,format='l,x,l,l,x,d,d,x,x,d', delimiter=',',plate,mjd,fiber,ra,dec,z
nnn=n_elements(plate)
result=replicate(struct, nnn)
result.sample='MD03'
result.plate=plate & result.fiber=fiber & result.mjd=mjd
result.ra=ra & result.dec=dec & result.z_pipe=z
result_all = [result_all, result]


; Other MD covered in DR10
file = outdir + 'DR10specMDSqsos.csv'
readcol,file,format='x,l,l,l,d,x,x,x,x,d,x,x,x,d,d',plate,mjd,fiber,z,z_vi,ra,dec,delimiter=','
nnn=n_elements(plate)
result=replicate(struct,nnn)
result.sample='MD_DR10'
result.plate=plate & result.fiber=fiber & result.mjd=mjd
result.ra=ra & result.dec=dec & result.z_pipe=z
result.z_vi = z_vi
result_all = [result_all, result]

; Stripe 82 quasars in DR7+DR10
file = outdir + 'DR10specS82qsos.csv'
readcol,file,format='x,l,l,l,d,x,x,x,x,x,d,d',plate,mjd,fiber,z,ra,dec,delimiter=','
nnn=n_elements(plate)
result=replicate(struct,nnn)
result.sample='S82_DR10'
result.plate=plate & result.fiber=fiber & result.mjd=mjd
result.ra=ra & result.dec=dec & result.z_pipe=z
result_all = [result_all, result]

; NOW, get z_VI from spall
;file = '/data1/quasar/yshen/data/sdss3/boss/v5_7_1_dr12/spAll-v5_7_0.fits'
;spall=hogg_mrdfits(file,1,columns=)

outfile = outdir + 'var_qso.fits'
mwrfits, result_all, outfile, /create
end
