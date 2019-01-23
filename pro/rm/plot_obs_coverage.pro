; plot the obs coverage of BOSS, bok, cfht and mayall

pro plot_obs_coverage

figfile=getenv('IDLRM_DIR')+'/etc/obs_cover.eps'
begplot, name=figfile,/color,/encap,/cmyk, xsize=14, ysize=6
charsize=2 & symthick=4

file=getenv('IDLRM_DIR')+'/etc/obs_cover.txt'

; imaging data
readcol, file,format='a,d,a,x,l,l,l,l,l,d',name, mjd, ut,obs_g,obs_i,obs_U,obs_r,obs_z,seeing
; boss data
readcol, file,format='x,d', boss_mjd, skipline=91L

; now plot
xrange=[56600, 56870]
plot, [0],[0], xrange=xrange, yrange=[0,1.1],xsty=5, /ysty,xtitle='Modified Julian Date [JD - 2400000.5 days]', $
 ytitle='Facility / Moon Illumination', charsize=charsize,pos = [0.08, 0.15, 0.98, 0.89], /norm, $
 yticklen=0.01, xticklen=0.03, $
 xtickname=['56600', '56650', '56700', '56750', '56800', '56850']
title='SDSS-RM Data Coverage'
xyouts, 0.5, 0.94, title,charsize=2, /norm,align=0.5

date=[[2013,12,01], $
      [2014,01,01], $
      [2014,02,01], $
      [2014,03,01], $
      [2014,04,01], $
      [2014,05,01], $
      [2014,06,01], $
      [2014,07,01], $
      [2014,08,01] ]
mjd0=dblarr(9)
mon=['2013 Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', '2014 Jul']
for i=0,8 do begin
  juldate, date[*,i], rjd  ; mjd=rjd-0.5
  mjd0[i]=rjd-0.5
  if i gt 0 then begin
    ; plot moon phase (illumination fraction)
    jdarr = mjd0[i-1] + findgen( (mjd0[i] - mjd0[i-1]) + 1 ) + 2400000.5
    mphase, jdarr, moon_illu
    ; make a shaded patch for the moon phase
    image=dblarr(n_elements(jdarr), 110)
    for j=0L, 109L do image[*,j]=moon_illu
    x_norm=convert_coord(jdarr - 2400000.5, replicate(0.,n_elements(jdarr)), /to_norm)
    x_norm=x_norm[0,*]
    y_norm=convert_coord([jdarr[0],jdarr[0]], [0,1.1], /to_norm)
    y_norm=y_norm[1,*]
    pos=[min(x_norm), min(y_norm), max(x_norm), max(y_norm)]
    ;message, 'stop'
    loadct, 0 ;, bottom=10  ;0
    ;cgimage, image, position=pos, /noerase
    ;tvscl, (image), pos[0],pos[1],xsize=pos[2]-pos[0],ysize=pos[3]-pos[1], /norm
    loadct, 0
    ;oplot, jdarr - 2400000.5, moon_illu, line=1
    ;xyouts, 0.5*(mjd0[i] + mjd0[i-1]), 1.02, mon[i-1],align=0.5,color=cgcolor('red')
  endif
  ;oplot, [mjd0[i], mjd0[i]], [0,1.1],line=2,color=cgcolor('cyan')
endfor
jdarr = mjd0[0] + findgen( (mjd0[8] - mjd0[0])/0.1 + 1 )*0.1 + 2400000.5
mphase, jdarr, moon_illu
image=dblarr(n_elements(jdarr), 110)
for j=0L, 109L do image[*,j]=moon_illu
x_norm=convert_coord(jdarr - 2400000.5, replicate(0.,n_elements(jdarr)), /to_norm)
x_norm=x_norm[0,*]
y_norm=convert_coord([jdarr[0],jdarr[0]], [0,1.1], /to_norm)
y_norm=y_norm[1,*]
pos=[min(x_norm), min(y_norm), max(x_norm), max(y_norm)]
loadct, 0 
tvscl, (image), pos[0],pos[1],xsize=pos[2]-pos[0],ysize=pos[3]-pos[1], /norm
loadct, 0
oplot, jdarr - 2400000.5, moon_illu, line=1
for i=0,8 do begin
  if i gt 0 then xyouts, 0.5*(mjd0[i] + mjd0[i-1]), 1.02, mon[i-1],align=0.5,color=cgcolor('red')
  oplot, [mjd0[i], mjd0[i]], [0,1.1],line=2,color=cgcolor('cyan')
endfor


axis, xaxis=0, charsize=charsize, xtickname=['56600', '56650', '56700', '56750', '56800', '56850'], /xsty, xrange=xrange, xtitle='Modified Julian Date [JD - 2400000.5 days]'
date2=[[2013, 12, 10], [2013, 12, 20], $
       [2014, 01, 10], [2014, 01, 20], $
       [2014, 02, 10], [2014, 02, 20], $
       [2014, 03, 10], [2014, 03, 20], $
       [2014, 04, 10], [2014, 04, 20], $
       [2014, 05, 10], [2014, 05, 20], $
       [2014, 06, 10], [2014, 06, 20], $
       [2014, 07, 10], [2014, 07, 20]]
mjd2=dblarr(16)
for i=0,15 do begin
  juldate, date2[*,i], rjd
  mjd2[i]=rjd-0.5
endfor
xtickname2=['10','20','10','20','10','20','10','20','10','20','10','20','10','20','10','20']
axis, xaxis=1, charsize=1, xtickv=mjd2, /xsty, xrange=xrange,xtickname=xtickname2,xticks=17
xyouts, 0.12,0.902, 'UT date',/norm,charsize=1, color=cgcolor('black')

; plot boss
nep=n_elements(boss_mjd)
oplot, boss_mjd, replicate(0.1, nep), psym=symcat(2), color=cgcolor('TG2')
xyouts, 56605., 0.1, 'BOSS',charsize=charsize

; plot bok
ind=where(name eq 'Bok' and obs_g eq 1,nnn)
oplot, mjd[ind], replicate(0.28, nnn), psym=symcat(16,thick=symthick),color=cgcolor('green')
ind=where(name eq 'Bok' and obs_g eq 2, nnn)
oplot, mjd[ind], replicate(0.28, nnn), psym=symcat(9,thick=symthick),color=cgcolor('green')
ind=where(name eq 'Bok' and obs_i eq 1, nnn)
oplot, mjd[ind], replicate(0.32, nnn), psym=symcat(16,thick=symthick),color=cgcolor('red')
ind=where(name eq 'Bok' and obs_i eq 2, nnn)
if nnn gt 0 then oplot, mjd[ind], replicate(0.32, nnn), psym=symcat(9,thick=symthick),color=cgcolor('red')
xyouts, 56605, 0.28, 'Bok',charsize=charsize

; plot cfht
ind=where(name eq 'CFHT' and obs_g eq 1,nnn)
oplot, mjd[ind], replicate(0.48, nnn), psym=symcat(14,thick=symthick),color=cgcolor('green')
ind=where(name eq 'CFHT' and obs_g eq 2, nnn)
oplot, mjd[ind], replicate(0.48, nnn), psym=symcat(4,thick=symthick),color=cgcolor('green')
ind=where(name eq 'CFHT' and obs_i eq 1, nnn)
oplot, mjd[ind], replicate(0.52, nnn), psym=symcat(14,thick=symthick),color=cgcolor('red')
ind=where(name eq 'CFHT' and obs_i eq 2, nnn)
oplot, mjd[ind], replicate(0.52, nnn), psym=symcat(4,thick=symthick),color=cgcolor('red')
xyouts, 56605, 0.48, 'CFHT',charsize=charsize

; plot Mayall
ind=where(name eq 'Mayall' and obs_u eq 1,nnn)
oplot, mjd[ind], replicate(0.88, nnn), psym=symcat(17,thick=symthick),color=cgcolor('BLU4')
ind=where(name eq 'Mayall' and obs_u eq 2, nnn)
if nnn gt 0 then oplot, [mjd[ind]], [replicate(0.88, nnn)], psym=symcat(5,thick=symthick),color=cgcolor('BLU4')
ind=where(name eq 'Mayall' and obs_g eq 1,nnn)
oplot, mjd[ind], replicate(0.84, nnn), psym=symcat(17,thick=symthick),color=cgcolor('green')
ind=where(name eq 'Mayall' and obs_g eq 2, nnn)
if nnn gt 0 then oplot, [mjd[ind]], [replicate(0.84, nnn)], psym=symcat(5,thick=symthick),color=cgcolor('green')
ind=where(name eq 'Mayall' and obs_r eq 1,nnn)
oplot, mjd[ind], replicate(0.80, nnn), psym=symcat(17,thick=symthick),color=cgcolor('cyan')
ind=where(name eq 'Mayall' and obs_r eq 2, nnn)
if nnn gt 0 then oplot, [mjd[ind]], [replicate(0.80, nnn)], psym=symcat(5,thick=symthick),color=cgcolor('cyan')
ind=where(name eq 'Mayall' and obs_i eq 1, nnn)
oplot, mjd[ind], replicate(0.76, nnn), psym=symcat(17,thick=symthick),color=cgcolor('red')
ind=where(name eq 'Mayall' and obs_i eq 2, nnn)
if nnn gt 0 then oplot, [mjd[ind]], [replicate(0.76, nnn)], psym=symcat(5,thick=symthick),color=cgcolor('red')
ind=where(name eq 'Mayall' and obs_z eq 1, nnn)
oplot, mjd[ind], replicate(0.72, nnn), psym=symcat(17,thick=symthick),color=cgcolor('magenta')
ind=where(name eq 'Mayall' and obs_z eq 2, nnn)
if nnn gt 0 then oplot, [mjd[ind]], [replicate(0.72, nnn)], psym=symcat(5,thick=symthick),color=cgcolor('magenta')
xyouts, 56605, 0.80, 'Mayall',charsize=charsize


xyouts, 56605, 0.95, 'U', charsize=charsize,color=cgcolor('BLU4')
xyouts, 56610, 0.95, 'g', charsize=charsize,color=cgcolor('green')
xyouts, 56614, 0.95, 'r', charsize=charsize,color=cgcolor('cyan')
xyouts, 56617, 0.95, 'i', charsize=charsize,color=cgcolor('red')
xyouts, 56620, 0.95, 'z', charsize=charsize,color=cgcolor('magenta')

;date=[[2013,12,01], $
;      [2014,01,01], $
;      [2014,02,01], $
;      [2014,03,01], $
;      [2014,04,01], $
;      [2014,05,01], $
;      [2014,06,01], $
;      [2014,07,01], $
;      [2014,08,01] ]
;mjd=dblarr(9)
;mon=['2013 Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul']
;for i=0,8 do begin
;  juldate, date[*,i], rjd  ; mjd=rjd-0.5
;  mjd[i]=rjd-0.5
;  oplot, [mjd[i], mjd[i]], [0,1.1],line=2
;  if i gt 0 then xyouts, 0.5*(mjd[i] + mjd[i-1]), 1.01, mon[i-1],align=0.5
;endfor




endplot

end
