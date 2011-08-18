;+
; NAME:
;    combinebpm
;
; PURPOSE:
;     Takes a set of darks (preferably at least 5 to reduce cosmics)
;     and a set of biases; combines them each individually, gets a
;     bias subtracted dark, and then makes a bad pixel mask.  Also
;     checks the bias for any hot pixels and marks those
;
; CALLING SEQUENCE:
; combinebpm,docams=,nsig=,ncount=,osx=,osy=,darkstart=,darkend=,biasstart=,biasend=,maxsat=,output=;
;
; INPUTS:
;
; docams     - camera to look at  ; default to ['b1','b2','r1','r2']
; darkstart   - first dark exposure; default to 01786
; darkend     - last  dark exposure; default to 01790
; biasstart   - first bias exposure; default to 01768
; biasend     - last  bias exposure; default to 01783
; nsig        - Number of sigma away from median to determine bad
;               pixel ;default to 5
; ncount      - Number of counts away from median to determine bad
;               pixel; default to 10
; output      - Write the 3 intermediary fits files,
;               masterdark, masterbias, and darkminusbias
; maxsat      - Any pixel that has a higher count than this in the averaged
;               biases is counted as a bad pixel, default to 1500.
; OPTIONAL KEYWORDS
;  
; OUTPUTS:
; 
; 
;
;      'bpm-sigma.fit'
;
;      'bpm-count.fit'
;
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   
;  This will combine darks and biases to make a bad pixel mask
;
;
;
; EXAMPLES:
;   Create a bad pixel mask for b1 using darks and biases from 55094
;   and write out intermediary steps
;
;
;    (assuming the files exist in /data/spectro/:)
;     IDL> combinebpm,docams='b1',nsig=10,ncount=15,osx=200,osy=200,darkstart=01786,darkend=01790,biasstart=01768,biasend=01783,/output
;
; 
;
; PROCEDURES CALLED:
;   sxpar()
;   sdssproc()
;   djs_iterstat()
; REVISION HISTORY:
;
;   22-Jan-2010  M. Olmstead, Utah
;-


;bit mask
;good 0
;hot pixel  (gt nsig or ncount)     1
;saturated hot pixel                2
;weak pixel  (lt nsig or ncount)    4
;bad column                         8
;bad warm column                    16
;hand added                         32
; for r2 55300
; combinebpm,docams='r2',darkstart=113290,darkend=113297,biasstart=113263,biasend=113285,/output
; for year 2 combinebpm,darkstart=118462,darkend=118472,biasstart=118538,biasend=118543,/output
; for year 2 combinebpm,darkstart=118462,darkend=118472,biasstart=119380,biasend=119395,/output

; for year 1    combinebpm,/output
pro combinebpm,docams=docams,nsig=nsig,ncount=ncount,darkstart=darkstart,darkend=darkend,biasstart=biasstart,biasend=biasend,maxsat=maxsat,output=output,indir=indir ;

if (NOT keyword_set(docams))   then docams =['b1','b2','r1','r2']
if (NOT keyword_set(nsig))      then nsig=10
if (NOT keyword_set(ncount))    then ncount=25;15
if (NOT keyword_set(darkstart)) then darkstart=101786
if (NOT keyword_set(darkend))   then darkend=101790
if (NOT keyword_set(biasstart)) then biasstart=101768
if (NOT keyword_set(biasend))   then biasend=101783
if (NOT keyword_set(maxsat))    then maxsat=15000.
if (NOT keyword_set(indir))     then begin
    indir = getenv('BOSS_SPECTRO_DATA')
    if (NOT keyword_set(indir)) then $
      indir='/data/spectro'
    indir = indir + '/*'
endif

ncam = n_elements(docams)
print,biasstart
 quiet = !quiet
   !quiet = 1

  if (ncam GT 1) then begin
      for icam=0, ncam-1 do begin
    combinebpm,docams=docams[icam],nsig=nsig,ncount=ncount,darkstart=darkstart,darkend=darkend,output=output,biasstart=biasstart,biasend=biasend,maxsat=maxsat,indir=indir
     endfor
      !quiet = quiet
      return
  endif


print,'dark step'
n_elementsdark=darkend-darkstart ;determining number of dark images to look at
filename =strarr(n_elementsdark+1)

for i=0, n_elementsdark do begin     ;defining filename of each image as a vector
filename(i)=darkstart+i     
 filename(i) = 'sdR-' + docams[0] +'-'+ string(filename(i), format='(i8.8)') $
    + '.fit*'
print,filename(i)
filename(i)=(findfile(djs_filepath(filename(i), root_dir=indir), count=ct))[0]
 if (ct EQ 0) then begin
      print, 'File ',filename(i),' not found'
      return
  endif
endfor

sdssproc,filename[0],imdark,hdr=hdr
if (docams eq 'b1' or docams eq 'b2') then begin
nx=4096
ny=4112
endif else begin
nx=4114
ny=4128
endelse

nxm1=nx-1
nym1=ny-1
mjd=sxpar(hdr,"MJD")

pixarrdark=fltarr(nx,ny,n_elementsdark+1) ;making a data cube 
for i=0,n_elementsdark do begin						
    sdssproc,filename[i],imdark
    pixarrdark[*,*,i]=imdark
endfor

for i=0, nxm1 do begin
    for j=0, nym1 do begin
        djs_iterstat,pixarrdark[i,j,*],median=median,mean=mean,sigma=sigma
        imdark[i,j]=median      ;picking median as best
    endfor
endfor



filenamedark='masterdark-'+docams+'-'+string(darkstart,format='(i8.8)')+'.fit'
if keyword_set(output) then  writefits,filenamedark,imdark,hdr


;================== now same thing but for bias



print,'bias step'
n_elementsbias=biasend-biasstart
filenamebias =strarr(n_elementsbias+1)

for i=0, n_elementsbias do begin ;defining filename vector
    filenamebias(i)=biasstart+i
    filenamebias(i) = 'sdR-' + docams[0] +'-'+ string(filenamebias(i), format='(i8.8)') + '.fit*'
    print,filenamebias(i)
    filenamebias(i)=(findfile(djs_filepath(filenamebias(i), root_dir=indir), count=ct))[0]
    if (ct EQ 0) then begin
        print, 'File ',filenamebias(i),' not found'
        return
    endif
endfor

sdssproc,filenamebias[0],imbias,hdr=hdr
bpmhot=imbias*0.
mjd=sxpar(hdr,"MJD")

pixarrbias=fltarr(nx,ny,n_elementsbias+1)
for i=0,n_elementsbias do begin						
    sdssproc,filenamebias[i],imbias
    pixarrbias[*,*,i]=imbias
endfor

for i=0, nxm1 do begin
for j=0, nym1 do begin
    djs_iterstat,pixarrbias[i,j,*],median=median,mean=mean,sigma=sigma
    imbias[i,j]=median ;picking median as best
    if imbias[i,j] gt maxsat then bpmhot[i,j]=2.         ;bad pixel mask for sat pixels
endfor
endfor

filenamebias='masterbias-'+docams+'-'+string(biasstart,format='(i8.8)')+'.fit'
writefits,filenamebias,imbias,hdr
if keyword_set(output) then writefits,filenamebias,imbias,hdr

;-----------------------now do bias subtraction
print,'bias subtraction step'
im_subtract=intarr(nx,ny)
im_subtract[*,*]=0LL
im_subtract=imdark-imbias+1000

filenamesubtract='darkminusbias-'+docams+'-'+string(darkstart,format='(i8.8)')+'-'+string(biasstart,format='(i8.8)')+'.fit'
if keyword_set(output) then writefits,filenamesubtract,im_subtract,hdr

;=======================now do bpm
print,'bpm step'
mdark=im_subtract
;nxm=2176
;nym=2112
;nyt=4224
;nxr=4352
mdark2=mdark
counts=0LL               ;checking for pixel below 50% of median

osx=0.
osy=0.
bpms=intarr(nx,ny)
bpms=bpms*0
bpmc=intarr(nx,ny)
bpmc=bpmc*0
nxm=nx/2.
nym=ny/2.


mdarktl=fltarr(nx/2-osx,ny/2-osy);(4061,2038)
mdarkbl=fltarr(nx/2-osx,ny/2-osy);(4061,1938)
mdarktr=fltarr(nx/2-osx,ny/2-osy);(4061,2038)
mdarkbr=fltarr(nx/2-osx,ny/2-osy);(4061,1938)
mdarktl2=fltarr(nx/2-osx,ny/2-osy);(4061,2038)
mdarkbl2=fltarr(nx/2-osx,ny/2-osy);(4061,1938)
mdarktr2=fltarr(nx/2-osx,ny/2-osy);(4061,2038)
mdarkbr2=fltarr(nx/2-osx,ny/2-osy);(4061,1938)
;mdarkv=fltarr(nx*ny)

for i=0,nx/2-osx-1 do begin;separating into 4 quadrants
for j=0,ny/2-osy-1 do begin
mdarktl[i,j]=mdark[i+osx,j+nym]
mdarktr[i,j]=mdark[i+nxm,j+nym]
mdarkbl[i,j]=mdark[i+osx,j+osy]
mdarkbr[i,j]=mdark[i+nxm,j+osy]
endfor
endfor

djs_iterstat,mdarktl,sigma=dsigtl,median=dmeantl
djs_iterstat,mdarktr,sigma=dsigtr,median=dmeantr
djs_iterstat,mdarkbl,sigma=dsigbl,median=dmeanbl
djs_iterstat,mdarkbr,sigma=dsigbr,median=dmeanbr
;stop

djs_iterstat,im_subtract,sigma=sigma,median=median,mean=mean

for i=0,nx-1 do begin
    for j=0,ny-1 do begin
        if (i lt nxm) then begin
            if (j lt nym) then begin ;  setting up correct bsig
                dsig=dsigbl
                dmean=dmeanbl
            endif else begin
                dsig=dsigtl
                dmean=dmeantl
            endelse
            if mdark[i,j] lt dmean-nsig*dsig then bpms[i,j]=1 ;sigma
            if mdark[i,j] gt dmean+nsig*dsig then bpms[i,j]=4
            if mdark[i,j] lt dmean-ncount then bpmc[i,j]=1 ;counts
            if mdark[i,j] gt dmean+ncount then bpmc[i,j]=4
        endif else begin
            if (j lt nym) then begin
                dsig=dsigbr
                dmean=dmeanbr
            endif else begin
                dsig=dsigtr
                dmean=dmeantr   ;
            endelse
;                                  
            if mdark[i,j] lt dmean-nsig*dsig then bpms[i,j]=1
            if mdark[i,j] gt dmean+nsig*dsig then bpms[i,j]=4
            if mdark[i,j] lt dmean-ncount then bpmc[i,j]=1
            if mdark[i,j] gt dmean+ncount then bpmc[i,j]=4
        endelse
        if (j lt osy or j gt ny-osy) then mdark[i,j]=0
        if (i lt osx or i gt nx-osx) then mdark[i,j]=0 ;gets overscan region
    endfor
endfor

;for i=0,nxm1 do begin
;    for j=0,nym1 do begin
;        if im_subtract[i,j] lt mean - nsig*sigma then bpms[i,j]=0   ;sigma
;        if im_subtract[i,j] gt mean + nsig*sigma then bpms[i,j]=0 
;        if im_subtract[i,j] lt mean - ncount then bpmc[i,j]=0   ;count
;        if im_subtract[i,j] gt mean + ncount then bpmc[i,j]=0 ;

;    endfor
;endfor


;add ticket 1368,1369, 1363  not listed from here

if docams eq 'b1' then  bpmhot[3734,1087:1344]+=16 ;blocked column only appears at low flux levels.  It blocks light for a fixed period of time, and then releases it creating a warm column. The signal does not appear in pixel flats or calibration flats, but does appear in science frames and arcs. For an example, see 55209/sdR-b1-00107423.fit.gz.  ds9 coordinates

if docams eq 'b2' then bpmhot[1833,2048:2226]+=16 ;same defect as b1  55186/sdR-b2-00105398.fit.gz ds9 coordinates


if docams eq 'r1' then bpmhot[3635:3652,397:416]+=32
if docams eq 'r1' then bpmhot[3810:3812,2278:2282]+=32

if docams eq 'r1' and  mjd lt 55413 then bpmhot[481:488,2141:2150]+=32
if docams eq 'r1' and  mjd ge 55413 then  bpmhot[692:710,3194]+=32
if docams eq 'r1' and  mjd lt 55413 then bpmhot[780:798,3467:3486]+=32
if docams eq 'r1' and  mjd lt 55413 then bpmhot[783:793,3468:3485]+=32
if docams eq 'r1' and  mjd lt 55413 then bpmhot[1618:1630,2047:2080]+=32
if docams eq 'r1' and  mjd lt 55413 then bpmhot[2166:2184,1362:1383]+=32
;if docams eq 'r1' and  mjd lt 55413 then bpmhot[2171:2178,1367:1371]+=32 ;kyles region larger
if docams eq 'r1' and  mjd lt 55413 then bpmhot[3635:3651,397:416]+=32
if docams eq 'r1' and  mjd ge 55413 then  bpmhot[3648:3649,395:2063]+=8


if docams eq 'r2' and mjd lt 55300 then bpmhot[2171:2192,2115:2217]+=32

if docams eq 'r2' and  mjd ge 55300 and  mjd lt 55413 then begin
    bpmhot[508,2064:4127]+=8         ;from looking at darks to get bad columns
    bpmhot[509:523,3474:3548]+=8;from looking at darks to get bad columns
    bpmhot[2984,2064:4127]+=8         ;from looking at darks to get bad columns;doesn't get all
    bpmhot[3956,2064:4127]+=8         ;from looking at darks to get bad columns;doesn't get all
    bpmhot[3882,2064:4127]+=8         ;from looking at darks to get bad columns;doesn't get all
    bpmhot[3865:3884,2636:3533]+=2         ;from looking at darks to get bad columns;doesn't get all; like ticket 1363
    bpmhot[3936:3958,2846:3106]+=2         ;from looking at darks to get bad columns;doesn't get all; like ticket 1363
    bpmhot[2619:2628,2232:2304]+=2  ;;ticket 1363
    bpmhot[2672,2064:2295]+=16 ;;r2, fiber~840. 8300 \aa\, sdssproc images, flag x=2673, 2065<y<2296 for data after 55300 (r2 replacement). This is a warm column that appears to vary in magnitude. See 55539/sdR-r2-00123636.fit.gz for an example. ds9 coordinates
endif

if docams eq 'r2' and  mjd ge 55413 then  bpmhot[508,2064:3513]+=8
;if docams eq 'r2' and  mjd ge 55413 then  bpmhot[509:520,3510:3512]+=32
if docams eq 'r2' and  mjd ge 55413 then  bpmhot[508:523,3510:3513]+=32
if docams eq 'r2' and  mjd ge 55413 then  bpmhot[2616:2629:2283]+=32
if docams eq 'r2' and  mjd ge 55413 then  bpmhot[2661:2673,2292]+=32
if docams eq 'r2' and  mjd ge 55413 then bpmhot[2984,2064:4127]+=32 ;from looking at darks to get bad columns;doesn't get all
if docams eq 'r2' and  mjd ge 55413 then bpmhot[3882,2064:4127]+=32 ;from looking at darks to get bad columns;doesn't get all
if docams eq 'r2' and  mjd ge 55413 then  bpmhot[3940:3958,2850:3096]+=32


bpms=bpms+bpmhot
bpmc=bpmc+bpmhot

n_bpms=n_elements(bpms);percentage of 0's
n_0=n_elements(where(bpms eq 0))
n_0c=n_elements(where(bpmc eq 0))

print,n_0c
print,n_bpms
print,'sigma '+docams+'          ',n_0*1./n_bpms
print,'counts '+docams+'         ',n_0c*1./n_bpms

bpms=bpms
bpmc=bpmc

if mjd lt 55413 then use='55025'
if docams eq 'r2' and mjd ge 55300 and mjd lt 55413 then use='55300'
if mjd ge 55413  then use='55413'  
filenamesigma='bpm-'+use+'-'+docams+'-sigma.fit'
writefits,filenamesigma,bpms,hdr
filenamecount='bpm-'+use+'-'+docams+'-count.fit'
writefits,filenamecount,bpmc,hdr

    
filenamesigma='bpm-'+use+'-'+docams+'-'+string(darkstart,format='(i8.8)')+'-'+string(biasstart,format='(i8.8)')+'-sigma.fit'
writefits,filenamesigma,bpms,hdr
filenamecount='bpm-'+use+'-'+docams+'-'+string(darkstart,format='(i8.8)')+'-'+string(biasstart,format='(i8.8)')+'-count.fit'
writefits,filenamecount,bpmc,hdr


;stop

end
