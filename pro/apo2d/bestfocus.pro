;+
; NAME:
;     bestfocus
;
; PURPOSE:
;   Compute the offset from Hartmann images of different focus ring positions
;
; CALLING SEQUENCE:
;   pro bestfocus  expstart=,expend=,red_focus1=, blue_focus1=, red_focus2=, blue_focus2=,camname=,indir=,skipcollimate=,viewplots=

;
; INPUTS:
;   expstart    - First exposure number of raw sdR file.
;   expend      - Last exposure number of raw sdR file
; OPTIONAL KEYWORDS
;   red_focus1  - The different red 1 focus ring positions; default to
;                 [-50.,-25.,0.,25.,50.,0.]
;   red_focus2  - The different red 2 focus ring positions; default to
;                red_focus1
;   blue_focus1  - The different blue 1 focus ring positions; default to
;                 red_focus1
;   blue_focus2  - The different blue 2 focus ring positions; default to
;                 red_focus1

;   camname    - Cameras to analyze; default to ['b1','b2','r1','r2'].
;   indir      - Input directory for files; default to searching for
;                files in $BOSS_SPECTRO_DATA/*.  If $BOSS_SPECTRO_DATA is not set,
;                then it is assumed to be /data/spectro.
;  
;   skipcollimate - If the Collimate*.log files already exist in the
;                  current directory, then you can skip collimate by
;                  setting skipcollimate=1.  By default, collimate
;                  will run
;
;   viewplots  - To look at each plot individually to determine what
;                data points to throw away, set viewplots=1.  By
;                default, you will not view each plot individually
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   
;   If required, Collimate is run on all files between expstart and expstop that go
;   in the sequence 'Left   ' 'Right   '.  The offsets in Collimate*.log files
;   are then used to figure out the offset in microns of the ccd as a
;   function of the focus ring position.  Since there is low number
;   statistics, the plots need to be looked at to see if the data is
;   good, many times the first and last focus ring position (very far
;   out of focus) is not a good fit.  The screen output is the out of
;   focus position from the mean for each subsection.  Going from left
;   to right is the focus offset in microns in the spatial direction
;   and top to bottom is the focus offset in microns in the spectral direction
;
;
;
; EXAMPLES:
;   Check for tilt on b2 with Hartmanns from 55218 
;    (assuming the files exist in /data/spectro/55218):
;     IDL> bestfocus, expstart=00107910,expend=00107921,camname='b2',red_focus1= [-50.,-25.,0.,25.,50.,0.]
;
; BUGS:
;
; PROCEDURES CALLED:
;   collimate
;   sdssproc
;   splog
;   sxpar()
;
; INTERNAL SUPPORT ROUTINES:
;   collimate_obscomm()
;   collimate_bad_hdr()
;
; REVISION HISTORY:
;
;   01-Sep-2009  K. Dawson, Utah
;   22-Jan-2010  M. Olmstead, Utah
;-



pro bestfocus  ,red_focus1=red_focus1, blue_focus1=blue_focus1, red_focus2=red_focus2, blue_focus2=blue_focus2,expstart=expstart,expend=expend,camname=camname,indir=indir, skipcollimate=skipcollimate,viewplots=viewplots

                                ;give focus ring positions and start and stop exposure numbers

if (NOT keyword_set(red_focus1)) then red_focus1= [-50.,-25.,0.,25.,50.,0.]          ;red 1 focus ring position
if (NOT keyword_set(blue_focus1)) then blue_focus1=red_focus1
if (NOT keyword_set(red_focus2)) then red_focus2=red_focus1
if (NOT keyword_set(blue_focus2)) then blue_focus2=red_focus1
if (NOT keyword_set(camname)) then camname =['b1','b2','r1','r2']             ;cameras to look at
if (NOT keyword_set(expstart)) then expstart=00107910                         ;First exposure while looking for Hartmanns
if (NOT keyword_set(expend)) then expend=00107921                             ;Last exposure while looking for Hartmanns
if (NOT keyword_set(skipcollimate)) then skipcollimate=0 ;If don't need to run collimate, set skipcollimate=1
if (NOT keyword_set(viewplots)) then viewplots=0 ;If you want to look at each plot used to find offset, set viewplots=1

ncam = n_elements(camname)

 quiet = !quiet                                     ;do for each camera
   !quiet = 1

  if (ncam GT 1) then begin
      for icam=0, ncam-1 do begin
    bestfocus,camname=camname[icam], red_focus1=red_focus1, blue_focus1=blue_focus1, red_focus2=red_focus2, blue_focus2=blue_focus2, expstart=expstart, expend=expend, indir=indir, skipcollimate=skipcollimate,viewplots=viewplots
     endfor
      !quiet = quiet
      return
  endif

 if (NOT keyword_set(indir)) then begin                        ;setting raw data directory
      indir = getenv('BOSS_SPECTRO_DATA')
      if (NOT keyword_set(indir)) then $
        indir='/data/spectro'
      indir = indir + '/*'
  endif



num=expend-expstart+1              ;total number of exposures
filename=strarr(num+1)

for i=0, num -1 do begin     ;defining filename of each image as a vector
filename(i)=expstart+i     
 filename(i) = 'sdR-' + camname[0] + '-' + string(filename(i), format='(i8.8)') $
    + '.fit*'
print,filename(i)
filename(i)=(findfile(djs_filepath(filename(i), root_dir=indir), count=ct))[0]
 if (ct EQ 0) then begin
      print, 'File ',filename(i),' not found'
      return
  endif
endfor

hexp=expstart                                ;want to keep good exposure numbers
goodones=fltarr(num)
hartmanns=0
hartmann1=strarr(1)
hartmann2=strarr(1)

for i = 0, num-2 do begin                        ;looking at headers, getting exposures that go LEFT, RIGHT
    hdr1=headfits(filename[i])
    hartmann1=sxpar(hdr1,'HARTMANN')
;print,hartmann1
    if (hartmann1 eq 'Left    ') then begin
        hdr2=headfits(filename[i+1])
        hartmann2=sxpar(hdr2,'HARTMANN')
        if hartmann2 eq 'Right   ' then begin
            goodones[hartmanns]=hexp ;ones that go left, right
            hartmanns=hartmanns+1
        endif
    endif
hexp=hexp+1
endfor

namegoodones=where(goodones ne 0)                             
numgoodones=n_elements(namegoodones)

if skipcollimate ne 1 then begin
for i=0, numgoodones-1 do begin
collimate,goodones[i],docams=camname,/debug
endfor
endif

offset=fltarr(numgoodones,8,8)          ;getting needed exposure numbers

mjd=sxpar(headfits(filename[0]),'MJD')

for i=0, numgoodones-1 do begin           ;reading in collimate*.log for LR pairs
readcol,'Collimate-'+string(mjd,format='(i5.5)')+'-'+camname+'-'+string(goodones[i],format='(i8.8)')+'.log',a,b1,b2,b3,b4,b5,b6,b7,b8,format='(a,f,f,f,f,f,f,f,f)'
offset[i,0,*]=b1
offset[i,1,*]=b2
offset[i,2,*]=b3
offset[i,3,*]=b4
offset[i,4,*]=b5
offset[i,5,*]=b6
offset[i,6,*]=b7
offset[i,7,*]=b8
endfor


if camname eq 'b1' then focus=blue_focus1
if camname eq 'b2' then focus=blue_focus2
if camname eq 'r1' then focus=red_focus1
if camname eq 'r2' then focus=red_focus2


; save offsets (in units of degrees on focus ring) as out(8,8)
offset2=fltarr(numgoodones)
out=fltarr(8,8)

for i=0,7 do begin
   for j=0,7 do begin
      plot,focus,offset[*,i,j],psym=2,ytitle='Distance between spots (pixels)',xtitle='degrees on focus ring for '+camname, $
                                       yran=[2,-2],symsize=2,charsize=1.5

;FLAG indices in offset vector, these fits are poor based on manual observations
      offset2[*]=offset[*,i,j]
                                ;offset[0]=0. ;if 0 bad by manual
                                ;observation, many times first and
                                ;last column are bad fits
     ;-------------------------------------------------------------------------------------------------------;
                                ;These are the ones that were not
                                ;counted for 55218 after going through
                                ;by hand

      if mjd eq 55218 then begin
          if camname eq 'b1' then begin
              offset2[3]=0
          endif
          if camname eq 'r1' then begin
              offset2[3]=0
          endif
          if camname eq 'b2' then begin                  ;b2 shows the most tilt and is clearly the worst
              if (i eq 0 and j eq 0) then offset2[3]=0 ;offset[3] = 25 degrees, offset[4]=50 degrees, offset [2] and [6] = 0 degrees
              if (i eq 1 and (j eq 0 or j eq 2 or j eq 4 or j eq 5 or j eq 6) ) then offset2[3]=0.
              if (i eq 2 and j ne 7) then offset2[3]=0.
              if (i eq 3 and (j eq 1 or j eq 2 or j eq 4 or j eq 5 or j eq 6 or j eq 7) ) then offset2[3]=0.
              if (i eq 3 and j eq 7) then offset[4]=0.
              if (i eq 4 and (j eq 1 or j eq 2 or j eq 3 or j eq 4 or j eq 5 or j eq 6)) then offset2[3]=0.0
              if (i eq 5 and j ne 7) then offset2[3]=0
              if (i eq 6 ) then offset2[3]=0
              if (i eq 7 and j eq 0) then offset2[3]=0. 
              if (i eq 4 and j eq 7) then offset2[0]=0.
              if (i eq 4 and j eq 7) then offset2[2]=0.
              if (i eq 4 and j eq 7) then offset2[4]=0.
              if (i eq 5 and j eq 7) then offset2[4]=0.
              if (i eq 5 and j eq 7) then offset2[1]=0.
              if (i eq 2 and j eq 7) then offset2[2]=0.
              if ((i eq 1 or i eq 2) and j eq 7) then offset2[5]=0.
              if (i eq 1 and j eq 7) then offset2[3]=0.
              if (i eq 1 and j eq 7) then offset2[4]=0.
              if ((i eq 4 and (j eq 6 )) or (i eq 5 and j eq 7) or (i eq 1 and j eq 7) or (i eq 3 and j eq 3)) then begin
                  plot,focus,offset[*,i,j],psym=2,ytitle='Distance between spots (pixels)',xtitle='degrees on focus ring', $
                    yran=[5,-5],symsize=2,charsize=1.5
              endif
          endif
      endif
;--------------------------------------------------------------------------------------------;

;FLAG all poor fits that have been assigned a value of 0
      ;a=where(offset2 ne 0.0 AND offset2 lt 2.5 AND offset2 gt -2.5)
      a=where(offset2 ne 0)
      result=linfit(focus[a],offset2[a],yfit=yfit)
      oplot,focus[a],yfit
      out[i,j]=-result[0]/result[1]
    ;Watch each fit, this is small number statistics, over-ride poor fits manually
      oplot,[out[i,j],out[i,j]],[0.,0.],psym=4,symsize=3
      if (viewplots eq 1) then stop

   endfor
endfor

outpos=(out-median(out))*25./11.25           ;offset in microns

print,camname
print,outpos

end

