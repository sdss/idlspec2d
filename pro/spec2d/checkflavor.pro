function checkflavor, image, flavor, camname
;
;	Purpose is to take global statistics of frame and decide
;	

	totalpix = n_elements(image)
	if (totalpix EQ 0) then return, 'unk'

	imagemed = median(image)

	binsize = 100
	hist = histogram(image,binsize=binsize,min=imagemed)
	totalhist = total(hist)
	nhist = n_elements(hist)

	if (totalhist EQ 0) then return, 'unk'
	
	over = fltarr(5)
	if (nhist GT 1) then over[0] = total(hist[1:*])/totalhist
	if (nhist GT 5) then over[1] = total(hist[5:*])/totalhist
	if (nhist GT 10) then over[2] = total(hist[10:*])/totalhist
	if (nhist GT 100) then over[3] = total(hist[100:*])/totalhist
	if (nhist GT 300) then over[4] = total(hist[300:*])/totalhist

;
;	guess flavors
;
	newflavor = 'unk'
	if (over[0] GT 0.02 AND over[1] GT 0.0005 AND $
             over[2] GT 0.0002) then newflavor = flavor
	if (over[0] GT 0.2 AND over[2] GT 0.04) then newflavor = 'arc'
	if (over[2] GT 0.5) then newflavor ='flat'
	if (over[4] GT 0.2) then newflavor ='satur'
	print, over, flavor, newflavor, camname

	return, flavor
end
