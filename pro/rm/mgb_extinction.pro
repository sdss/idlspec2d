; Dec 13, 2004
; Find E(B-V) at the position given.
; Calculate reddening as a function of wavelength
; Output reddening spectrum for each object
; inputs ra, dec in degrees, wavel in angstroms
; wfind = 0 means return the entire spectrum
; rv=3.1 or 3.08 for average Milky Way extinction curve
; this uses the Cardelli, Clayton, Mathis (1989) extinction curve (CCM)
; for optical/NIR only: 3030-9091 A
; is the wavelength in vacuum or air?

; return, magnitude extinction A_lam at each wavelength
; dereddened flux = flux*10.0D^(0.4*A_lam)

function mgb_extinction,ra,dec,rv,wavel,wfind, eb_v=eb_v

 nnn = n_elements(wavel)

 if n_elements(eb_v) eq 0 then begin ; get E(B-V) from target ra/dec
   ;;get Galactic coordinates
   glactc, ra, dec, 2000, gl, gb, 1, /DEGREE
   ;;get reddening
   ; ipath="/u/jgreene/misc/extinction/"
   ; ipath="/u/yshen/Research/misc/extinction/"
   ; ipath = "/home/yshen/Research/misc/maps/"
   eb_v = dust_getval(gl,gb,ipath=ipath,/interp)
 endif
 ;;calculate A_lam
 wavel_micron=wavel*1.e-4
 xarr=1./wavel_micron
 yarr=xarr-1.82 ; used to be a typo: y=x-1.86 

 a_x = dblarr(nnn) & b_x = dblarr(nnn)

 ; Infrared
 ind = where(xarr ge 0.3 and xarr le 1.1)
 if ind[0] ne -1 then begin
   a_x[ind] = 0.574*xarr[ind]^1.61 & b_x[ind] = -0.527*xarr[ind]^1.61
 endif

; optical/NIR
ind = where(xarr gt 1.1 and xarr le 3.3)
if ind[0] ne -1 then begin
   y = yarr[ind]
   a_x[ind]=1.+0.17699*y-0.50447*y^2-0.02427*y^3+0.72085*y^4+0.01979*y^5-$
     0.77530*y^6+0.32999*y^7
   b_x[ind]=1.41338*y+2.28305*y^2+1.07233*y^3-5.38434*y^4-0.62251*y^5+$
     5.30260*y^6-2.09002*y^7
endif

 a_lam=eb_v*(rv*a_x+b_x)
 if wfind eq 0 then return,a_lam
 if wfind gt 0 then begin
     reg=where(wavel ge wfind-5. and wavel le wfind+5.)
     a_find=mean(a_lam[reg])
     return,a_find
 endif

end

