; set variables
path= '/data/spectro/veldisp/brgtemplate'
mjd= 51609
plate= 304
suffix= '-'+string(mjd,format='(I5)')+'-0'+string(plate,format='(I3)')+'.fits'

; read and write
veldisp= mrdfits(path+'/spVelDisp'+suffix,1)
readspec, plate,tsobj=tsobj
help, tsobj
fiber= veldisp.fiber
tsobjtrim= tsobj[(fiber-1)]
mwrfits, tsobjtrim,path+'/tsObj-trim'+suffix,/create

; done
;exit
