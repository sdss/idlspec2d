; output the qsofits for the SDSS-RM objects that Yanfei finds interesting

pro output_yanfei_objfits

; get the master file
file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
target=mrdfits(file,1)

; read in the list of objects to be compiled
file='/data3/quasar/yshen/work/photo_lag/GoodList.txt'
readcol,file,format='l,d,d',id1,ra1,dec1
file='/data3/quasar/yshen/work/photo_lag/CheckList.txt'
readcol,file,format='l,d,d',id2,ra2,dec2
file='/data3/quasar/yshen/work/photo_lag/Z_list.txt'
readcol,file,format='x,d,d', ra3, dec3
id3=replicate(0, n_elements(ra3))

id=[id1,id2,id3]
ra=[ra1,ra2,ra3] & dec=[dec1,dec2,dec3]

; Yanfei provided very crude coordinates
spherematch, ra,dec, target.ra, target.dec, 5./3600.D, match0, match1, distance

; find the list of objects
ind=sort(match0)
target=target[match1[ind]]

plate=(target.plate)[0:31,*]
fiber=(target.fiberid)[0:31,*]
mjd=(target.mjd)[0:31,*]

nnn=32

finaldir='/data3/quasar/yshen/work/photo_lag/jiang/'
if file_test(finaldir,/dir) eq 0 then spawn, 'mkdir ' + finaldir
linename = ['Halpha', 'Hbeta', 'OIII5007','MgII', 'CIII', 'CIV']
for i=0L, nnn - 1 do begin

  outdir_in=getenv('BOSS_SPECTRO_REDUX') + '/' + getenv('RUN2D') $
       + '/' + string(plate[i,0],format='(i4.4)') + '/' $
       + 'qsofit/'
  platestr=string(plate[i,0],format='(i4.4)')
  mjdstr=string(mjd[i,0],format='(i5.5)')

  outfile='jiang-'+platestr+'-'+mjdstr+'.fits'

  rm_compile_qsofit, reform(plate[i,*]), reform(mjd[i,*]), $
     input_fiber=reform(fiber[i,*]),rm_plate=0,rm_ID=id, $
     outdir_in=outdir_in,outfile=outfile, linename=linename

  spawn, 'mv ' + outdir_in+outfile + ' ' + finaldir
endfor


end
