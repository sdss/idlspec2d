sdssproc,'sdR-01-00000384.fit',image,invvar
 medbin, image, binim, 8
xcen = trace_crude(binim, yset=ycen, ystart=10)

; Fixing up some bad traces
xcen[*,17] = xcen[*,19] - 12.3
xcen[*,15] = xcen[*,19] - 18.4
xcen[*,274] = xcen[*,275] - 6.1

x2 = fltarr(128,317)
y2 = fltarr(128,317)
x2[*,0:15] = xcen[*,0:15]
x2[*,16:173] = xcen[*,17:174]
x2[*,174:316] = xcen[*,176:318]
y2[*,0:15] = ycen[*,0:15]
y2[*,16:173] = ycen[*,17:174]
y2[*,174:316] = ycen[*,176:318]


;  Extraction of binned image
bininv = 1.0/(abs(binim) +5.0)
extract_image,binim,bininv,x2,1.2,fluxbin,errorbin

xy2traceset,y2,x2,tset
tset.xmax=2047
traceset2xy,tset,yfull,xfull

;  Full extraction
extract_image,image,invvar,xfull,1.2,flux,error

;  Fancier extraction
;  Takes about 4 minutes to run
extract_image,image,invvar,xfull,1.2,flux,error,ymodel=ymodel,wfixed=[1,1,1]

; Look at residuals in fitting

atv,(image-ymodel)*(image-ymodel)*invvar
