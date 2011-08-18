;+
; NAME:
;    flatfield
;
; PURPOSE:
;    This takes an LED flatfield and compares pixel by pixel the value
;    of that pixel compared to the median in a subregion about it
;
; CALLING SEQUENCE:
; 
;flatfield, filename,filenamebias,docams=
; INPUTS:
;
; docams      - camera to look at  ; default to ['b1','b2','r1','r2']
; filename     - name of flat
; filenamebis  - name of bias
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
;  This will take an LED flat and a bias to make a flatfield
;
;
;
; EXAMPLES:
;   Create a flatfield using LED flats from 55141/55142
;
;
;    (assuming the files exist in /data/spectro/:)
;     IDL> flatfield,
;
; 
;
; PROCEDURES CALLED:
;   sxpar()

;
; REVISION HISTORY:
;
;   22-Jan-2010  M. Olmstead, Utah
;-




pro flatfield, expflat,expbias,docams=docams,indir=indir

if (NOT keyword_set(docams))   then docams =['b1','b2','r1','r2']
if (NOT keyword_set(indir))     then begin                         ;set BOSS_SPECTRO_DATA
    indir = getenv('BOSS_SPECTRO_DATA')
    if (NOT keyword_set(indir)) then $
      indir='/data/spectro'
    indir = indir + '/*'
endif

ncam = n_elements(docams)

quiet = !quiet
!quiet = 1

if (ncam GT 1) then begin
    for icam=0, ncam-1 do begin                                    ;cycle through cameras
        flatfield,expflat,expbias,docams=docams[icam],indir=indir
    endfor
    !quiet = quiet
    return
endif

filename = 'sdR-' + docams[0] + '-' + string(expflat, format='(i8.8)')  + '.fit*'
filenamebias = 'sdR-' + docams[0] + '-' + string(expbias, format='(i8.8)')  + '.fit*'            ;get files
filename = (findfile(djs_filepath(filename, root_dir=indir), count=ct1))[0]
filenamebias = (findfile(djs_filepath(filenamebias, root_dir=indir), count=ct2))[0]

if (ct1 EQ 0 ) then begin
    print, 'flat file not found'
stop
    return
endif
if (ct2 EQ 0 ) then begin
    print, 'bias file not found'
    return
endif


imflat2=readfits(filename,hdr)
imbias2=readfits(filenamebias)
imflat3=imflat2-imbias2                                   ;bias subtract
location=0.
nx=sxpar(hdr,'NAXIS1')
ny=sxpar(hdr,'NAXIS2')
outfitosbias=fltarr(nx,ny)
imflatosbias=fltarr(nx,ny)

if docams eq 'b1' then  gain = [1.048, 1.048, 1.018, 1.006] ; b1 gain
if docams eq 'b2' then  gain = [1.040, 0.994, 1.002, 1.010] ; b2 gain
if docams eq 'r1' then  gain = [1.966, 1.566, 1.542, 1.546] ;r1 gain
if docams eq 'r2' then  gain = [1.956, 1.618, 1.538, 1.538] ;r2 gain



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; pixel by pixel bias subtraction
;imflat[0:2175,0:2111]=imflat3[0:2175,0:2111]*gain[0]
;imflat[2176:4351,0:2111]=imflat3[2176:4351,0:2111]*gain[1]
;imflat[0:2175,2112:4223]=imflat3[0:2175,2112:4223]*gain[2]
;imflat[2176:4351,2112:4223]=imflat3[2176:4351,2112:4223]*gain[3]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;bias subtraction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;of overscan 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;after pixel by pixel

   os2 = [ median(imflat3[10:67,56:2111]), $
                  median(imflat3[4284:4340,56:2111]), $
                  median(imflat3[10:67,2112:4167]), $
                  median(imflat3[4284:4340,2112:4167]) ]


imflatosbias[0:2175,0:2111] = gain[0] * (imflat3[0:2175,0:2111] - os2[0])
imflatosbias[2176:4351,0:2111] = gain[1] * (imflat3[2176:4351,0:2111] - os2[1])
imflatosbias[0:2175,2112:4223] = gain[2] * (imflat3[0:2175,2112:4223] - os2[2])
imflatosbias[2176:4351,2112:4223] = gain[3] * (imflat3[2176:4351,2112:4223] - os2[3])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
imflat2=0
imbias2=0

;=====================================================================
mx=10.                             ;number of pixels or left and right for median count
my=10
counti=210                       
countj=210
;==================================
sampleosbias=fltarr(2*mx+1,2*my+1)

for i=200,nx-200 do begin        ;making a rectangle around the point and median comparing it
    for j=200,ny-200 do begin
        sampleosbias=imflatosbias(i-mx:i+mx,j-my:j+my)
        djs_iterstat,sampleosbias,median=medianosbias
        outfitosbias[i,j]=imflatosbias[i,j]/medianosbias
    endfor
    if (i gt 100+location) then begin    ;a status counter
        location=location+100
        print,i
    endif
endfor

if docams eq 'b2' then begin               ;setting 3 bad columns = 0 for b2
    for j=1693,2111 do begin
        outfitosbias[2094,j]=0
    endfor

    for j=2112,3457 do begin
        outfitosbias[2621,j]=0
    endfor

    for j=1392,2111 do begin
        outfitosbias[462,j]=0
    endfor  
endif

imout=outfitosbias[119:4232,48:4175]
filenameout='pixflat-'+docams+'.fit'
writefits, imout,outfitosbias

;

;============================================================================




end
