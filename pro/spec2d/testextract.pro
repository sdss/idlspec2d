;---------------
; Run in directory with 15may99 data

; Read image and make crude traces
sdssproc,'sdR-01-00000384.fit',image,invvar
xcen = trace_crude(image, yset=ycen, nmed=5, nave=21)

; Display image with traces
ntrace = (size(xcen,/dim))[1]
atv,image
for i=0,ntrace-1 do atvplot,xcen[*,i],ycen[*,i]

; Fix traces and overplot
xnew = trace_fix(xcen)
for i=0,ntrace-1 do atvplot,xnew[*,i],ycen[*,i], color='green'

xy2traceset,ycen,xnew,tset
traceset2xy,tset,ycen,xsol

; Optimal extraction
extract_image,image,invvar,xsol,1.0,flux,error,ymodel=ymodel
;---------------

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

createfluxstruct, tt, flux1136, finv1136, tset, invset, plugmap, 1
sky = where(tt.plugmap.objtype EQ 'SKY')
fitarcimage, flux1148, 'blue', lamplist, xpeak, ypeak, lambda, tset, $ 
   invset, ans=ans


for i=0,319 do skyeach[*,i] = skycombine(skystruct, tt[i].coeff, skyspline=skyspline)

;;;;;;;;;;;;;;;;;;;;  Extracting blue1 ;;;;;;;;;;;;;;;;;;;;;;;;;

sdssproc,'sdR-01-00001151.fit',image11511, var11511
sdssproc,'sdR-01-00001136.fit',image11361,var11361,hdr=hdr11361
sdssproc,'sdR-01-00001148.fit',image11481, var11481
satval = 62000.0
var11361 = var11361 * (image11361 LT satval)
xcen11511 = trace320crude(image11511, yset=ycen11511)
xy2traceset,ycen11511,xcen11511,tset11511,ncoeff=4
traceset2xy,tset11511,ycen11511,xsol11511

extract_image,image11361,var11361,xsol11511,1.0,flux11361,finv11361,ymodel=ymodel11361,ansimage=ansimage11361,proftype=4,wfixed=[1,1,1],highrej=15

extract_image,image11361,var11361,xsol11511,1.0,flux11362,finv11362,ymodel=,ymodel=ymodel11362,ansimage=ansimage11362,proftype=2

extract_image,image11481,var11481,xsol11511,1.0,flux11481,finv11481,proftype=2

readcol,getenv('EVIL_SRC')+'spec2d/lamphgcdne.dat',lampwave,lampinten
.com lampfit
lamplist = [[lampwave],[lampinten]]
fitarcimage, flux11481, 'blue', lamplist, xpeak, ypeak, wset11481, invset11481, ans=ans

yanny_read, 'plPlugMapM-0198-51433-01.par', pdata, hdr=hdrplug
plugmap = *pdata[0]
createfluxstruct, tt11361, flux11361, finv11361, wset11481, invset11481, plugmap, 1

ttsub11361 = skysubtract(tt11361)

writespectra, tt[0:1], hdr11361,

restore,'test1.dat'
row = 708
aa = extract_row(image11361[*,row], var11361[*,row], xsol11511[row,*], 1.0, ymodel=ymodel, proftype=5, wfixed=[1,1,1,1,1], highrej=30)
plot,image11361[*,row],xr=[860,920],yr=[0,1000],ps=10
djs_oplot,ymodel,color='red',ps=1
print,aa[*,138]
plot,(image11361[*,row]-ymodel6)*sqrt(var11361[*,row]),xr=[880,900]
plot,image11361[*,708],xr=[880,900],ps=10
plot,image11361[*,708],xr=[880,900],ps=10,yr=[0,1000]
djs_oplot,ymodel6,color='green',ps=1  

yrow=lindgen(205)*10+5
extract_image,image11361,var11361,xsol11511,1.0,f,finv,ymodel=ymodel, proftype=5, wfixed=[1,1,1,1,1], highrej=30,lowrej=30,ansimage=ansimage,yrow=yrow

extract_image,image11361,var11361,xsol11511,1.0,f,finv,ymodel=ymodel3, proftype=4, wfixed=[1,0,0,0], highrej=15,lowrej=30,fitans=fitans
extract_image,image11361,var11361,xsol11511,1.0,f,finv,proftype=4, wfixed=[1,1,1,1]
aa = extract_row(image11361[*,row], var11361[*,row], xsol11511[row,*], 1.0, ymodel=ymodel, proftype=4, wfixed=[1,0,0,0], highrej=30,/squashprofile,inputans=fitans[*,*,row])


;
;	Call extraction routine
;

doaframe, 'sdR-04-00001151.fit','sdR-04-00001148.fit',['sdR-04-00001136.fit'], 'plPlugMapM-0198-51433-01.par',outputDir='out'

tt = doaframe('sdR-04-00001151.fit','sdR-04-00001148.fit',['sdR-04-00001136.fit'], 'plPlugMapM-0198-51433-01.par',outputDir='/ide_disk/51433/out', inputDir='/ide_disk/51433')


plot,10^(poly(pixarray,ans)),spec,xr=[5800,7000]
djs_oplot,linelist[*,0],linelist[*,1]/50,ps=1,color='red'
djs_oplot,10^lambda,lambda*0.0+400,ps=4
djs_oplot,10^(poly((2.0*xcen[row,*]-2047.0)/2047.0,ans)),xcen[row,*]*0.0+500,ps=6


IDL> 
tt = doframe('sdR-04-00001151.fit','sdR-04-00001148.fit',$
     ['sdR-04-00001136.fit'],'plPlugMapM-0198-51433-01.par',$
       inputDir='/s1/data/SDSS/51433',outputDir='/s1/data/SDSS/51433/out')
objectnums = ['-00001136','-00001138']
objectnums = ['-00001435','-00001437','-00001439'] ; skyflats
objectnums = ['-00001459','-00001461','-00001463','-00001465'] ;plate 202
flatnum = '-00001431'
flatnum = '-00001473'
arcnum = '-00001427'
arcnum = '-00001467'

cam = '04'
flatnum = '-00001451'
flat = 'sdR-'+cam+flatnum+'.fit'
arcnum = '-00001449'
arc = 'sdR-'+cam+arcnum+'.fit'
objectnums = ['-00001441','-00001443','-00001445','-00001447'] ;plate 214
objects = 'sdR-'+cam+objectnums+'.fit'
plugMapFile = 'plPlugMapM-0202-51434-01.par'
inputDir='/s1/data/SDSS/51441'
outputDir = '/s1/data/SDSS/51441/out'
plugMapDir = inputDir+'/logs'
tt = doaframe(flat,arc,objects,plugMapFile,inputDir=inputDir, $
 outputDir=outputDir, plugMapDir=plugMapDir)

