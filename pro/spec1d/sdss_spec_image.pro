;+
; NAME:
;  sdss_spec_image
; PURPOSE:
;  Create an image of a spectrum
; CALLING SEQUENCE:
;  sdss_spec_image, outbase, plate, fiber, mjd=, xra=
; INPUTS:
;  outbase - output file base
;  plate - plate #
;  fiber - fiber #
;  mjd - MJD of observation
; OPTIONAL INPUTS:
;  xra - [2] wavelength limits to plot (Ang)
; COMMENTS:
;  Creates the files:
;     [outbase].jpg
;     [outbase].thumb.jpg
;  which have a picture of the specified spectrum,
;  with some nice annotation.
; REVISION HISTORY:
;  Written by MRB, NYU 14-05-2010
;------------------------------------------------------------------------------
function sdss_spec_smooth, loglam,flux,vdisp

if vdisp GT 1.0 then begin
    nlambda= n_elements(loglam)
    pixsize= $
      abs(alog(10.)*2.99792e+5*(loglam[nlambda-1]-loglam[0])/double(nlambda))
    smoothing= vdisp/pixsize    ; pixels
    npix= long(4.0*ceil(smoothing))*2L+3
    klam= findgen(npix)-float(npix-1.)/2.
    kernel= exp(-0.5*(klam/smoothing)^2)/sqrt(2.*!DPI)/smoothing
    kernel= kernel/total(kernel)
    if(n_elements(kernel) ge n_elements(flux)) then $
      smoothed_spec= flux $
    else $
      smoothed_spec= convol(flux,kernel,/edge_truncate)
endif else begin
    smoothed_spec= flux
endelse

return, smoothed_spec

end
;;
pro sdss_spec_image, outbase, plate, fiber, mjd=mjd, run2d=run2d, run1d=run1d, $
                     topdir=topdir, xra=xra, silent=silent

common com_sdss_spec_image, plans

if(n_tags(plans) eq 0) then $
  plans= yanny_readone(getenv('PLATELIST_DIR')+'/platePlans.par')

sscale=1.5
if(NOT keyword_set(xsize)) then xsize= 10.5*sscale
if(NOT keyword_set(ysize)) then ysize= 7.5*sscale
if(NOT keyword_set(xra)) then xra= [3501., 10499.]

if(n_elements(outbase) eq 0 OR $
   n_elements(plate) eq 0 OR $
   n_elements(fiber) eq 0) then $
  message, 'Usage: sdss_spec_image, outbase, plate, fiber, mjd='

if(n_elements(outbase) gt 1 OR $
   n_elements(plate) gt 1 OR $
   n_elements(fiber) gt 1) then $
  message, 'OUTBASE, PLATE, FIBER must be scalar'

readspec, plate, fiber, mjd=mjd, zans=zans, flux=flux, wave=wave, $
  invvar=invvar, run2d=run2d, run1d=run1d, topdir=topdir, plug=plug, $
  silent=silent

ipl= where(plate eq plans.plateid, npl)

igd= where(invvar gt 0, ngd)
if(ngd gt 0) then begin
    ist= min(igd)
    ind= max(igd)
    flux= flux[ist:ind]
    invvar= invvar[ist:ind]
    wave= wave[ist:ind]
endif

if(NOT keyword_set(axis_char_scale)) then axis_char_scale= 1.75*sscale
if(NOT keyword_set(tiny)) then tiny=1.d-4

bangp=!P
bangx=!X
bangy=!Y
!P.FONT= -1
set_plot, "PS"
!P.BACKGROUND= djs_icolor('white')
!P.COLOR= djs_icolor('black')

device, file=outbase+'.ps',/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color
!P.THICK= 2.0*sscale
!P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
!P.CHARSIZE= 1.0
!P.PSYM= 0
!P.LINESTYLE= 0
!P.TITLE= ''
!X.STYLE= 1
!X.CHARSIZE= axis_char_scale
!X.MARGIN= [1,1]*0.5
!X.OMARGIN= [7,7]*axis_char_scale
!X.RANGE= 0
!X.TICKS= 0
!Y.STYLE= 1
!Y.CHARSIZE= !X.CHARSIZE
!Y.MARGIN= 0.6*!X.MARGIN
!Y.OMARGIN= 0.6*!X.OMARGIN
!Y.RANGE= 0
!Y.TICKS= !X.TICKS
!P.MULTI= [1,1,1]
xyouts, 0,0,'!6'
colorname= ['red','green','blue','magenta','cyan','dark yellow', $
            'purple','light green','orange','navy','light magenta', $
            'yellow green']
ncolor= n_elements(colorname)
loadct,0,silent=silent

sflux= sdss_spec_smooth(alog10(wave), flux, 100.)

igd= where(invvar gt 0 and abs(wave-5577.) gt 4. and wave lt 10000., ngd)

if(ngd gt 0) then $
  yra= minmax(sflux[igd]) $
else $
  yra= minmax(sflux)
size= 0.07*(yra[1]-yra[0])
yra=yra+[-1.2,1.7]*size*1.7
if(yra[0] lt -2.) then yra[0]=-1.999

djs_plot, [0], [0], xra=xra, yra=yra, /nodata, $
  xtitle='!6Wavelength (Angstroms)!6', $
  ytitle='!8f_\lambda !6(10^{-17} erg/s/cm^2/Ang)!6'

if(npl gt 0) then begin
    title0= 'Survey: !8'+strtrim(plans[ipl].survey,2)+ $
      '!6 Program: !8'+strtrim(plans[ipl].programname,2)
    targets=''
    if(plans[ipl].survey eq 'sdss') then $
      targets= strtrim(strjoin(sdss_flagname('TARGET', plug.primtarget),' '),2)
    if(plans[ipl].survey eq 'segue1') then $
      targets= strtrim(strjoin(sdss_flagname('SEGUE1_TARGET', plug.primtarget),' '),2)
    if(plans[ipl].survey eq 'boss') then $
      targets= strtrim(strjoin(sdss_flagname('BOSS_TARGET1', plug.boss_target1),' '),2)+ $
      ' '+strtrim(strjoin(sdss_flagname('ANCILLARY_TARGET1', plug.ancillary_target1),' '),2)
    title0= title0+' !6Target: !8'+targets+'!6'
endif
    
title1= 'RA='+strtrim(string(f='(f40.5)', zans.plug_ra),2)+', '+ $
  'Dec='+strtrim(string(f='(f40.5)', zans.plug_dec),2)+', '+ $
  'Plate='+strtrim(string(zans.plate),2)+', '+ $
  'Fiber='+strtrim(string(zans.fiberid),2)+', '+ $
  'MJD='+strtrim(string(zans.mjd),2)

if(zans.z lt 1000./299792.) then $
  zstr= '!8cz='+strtrim(string(long(zans.z*299792.)),2)+'+/-'+ $
  strtrim(string(long(zans.z_err*299792.)),2)+'!6 km/s' $
else $
  zstr= '!8z='+strtrim(string(f='(f40.5)',zans.z),2)+'\pm!8'+ $
  strtrim(string(f='(f40.5)',zans.z_err),2) +'!6'
title2= zstr+' Class='+strtrim(zans.class)+' '+strtrim(zans.subclass)

warnings= strtrim(strjoin(sdss_flagname('ZWARNING', zans.zwarning),' '),2)
;;; print, zans.zwarning
if(keyword_set(warnings) gt 0) then $
    title3= 'Warnings: '+warnings  $
else $
  title3='No warnings.'



xst= !X.CRANGE[0]+0.01*(!X.CRANGE[1]-!X.CRANGE[0])
if(keyword_set(title0)) then begin
    yst= !Y.CRANGE[0]+1.14*(!Y.CRANGE[1]-!Y.CRANGE[0])
    xyouts, xst, yst, title0, charsize=1.2*sscale
endif
yst= !Y.CRANGE[0]+1.10*(!Y.CRANGE[1]-!Y.CRANGE[0])
djs_xyouts, xst, yst, title1, charsize=1.2*sscale
yst= !Y.CRANGE[0]+1.06*(!Y.CRANGE[1]-!Y.CRANGE[0])
djs_xyouts, xst, yst, title2, charsize=1.2*sscale
yst= !Y.CRANGE[0]+1.02*(!Y.CRANGE[1]-!Y.CRANGE[0])
xyouts, xst, yst, title3, charsize=1.2*sscale

err= fltarr(n_elements(sflux))
igd= where(invvar gt 0, ngd)
if(ngd gt 0) then $
  err[igd]= sqrt(1./invvar[igd])
ibd= where(invvar le 0, nbd)
if(nbd gt 0) then $
  err[ibd]= (yra[1]-yra[0])
fillx= [wave, reverse(wave)]
filly= [sflux+err, reverse(sflux-err)]
fillx= [fillx, fillx[0]]
filly= [filly, filly[0]]

polyfill, fillx, filly, noclip=0, color=djs_icolor('light gray')

djs_oplot, wave, sflux, th=6

;;wave= [4160., 4227., 4300., 4383., 4455., 4531., 4668., 5015., 5100., 5175., 5270., 
;;5335., 5406., 5709., 5782., 5895., 5970., 6230. ], 
;;name= ["CN", "Ca4227", "G4300", "Fe4385", "Ca4455", "Fe4531", "Fe4668", "Fe5015", 
;;"Mg1", "Mg b", "Fe5270", "Fe5335", "Fe5406", "Fe5709", "Fe5782", "Na D", 
;;"TiO1", "TiO2"],

awave= [ 4300., 5895., 5175., 8498., 8542., 8662., 3968., 3938.]
aname= [ "G", "Na D", "Mg", "CaII", "", "", "H", "K"]

ewave= [3727., 3869.7867, 4105.8884, 4341.6803, 4364.3782, $
        4862.6778, 4960.2140, 5008.1666, 5876., 6301.9425, $
        6549.7689, 6564.6127, 6585.1583, 6718.1642, 6732.5382, $
        7137.6370, 2800., 1216., 1549., 1640., 1909., 2326., 1400. ]
        
ename= textoidl(['OII', 'NeIII', 'H\delta', 'H\gamma', 'OIII', $
                 'H\beta', '', 'OIII', 'HeI', 'OI', 'NII', $
                 'H\alpha', 'NII', '', 'SII', 'ArIII', 'Mg', 'Ly\alpha', $
                 'CIV', 'HeII', 'CIII', 'CII', 'SIV+OIV'])

xsize= 0.07*(xra[1]-xra[0])

yoff= replicate(0., n_elements(ewave))
xoff= replicate(0., n_elements(ewave))
ioff= where(abs(ewave-6585.) lt 3., noff)
xoff[ioff]=0.3*xsize
yoff[ioff]=-0.5*size
ioff= where(abs(ewave-6549.) lt 3., noff)
xoff[ioff]=-0.3*xsize
yoff[ioff]=-0.5*size
ioff= where(abs(ewave-4862.) lt 3., noff)
xoff[ioff]=-0.1*xsize
ioff= where(abs(ewave-5008.) lt 3., noff)
xoff[ioff]=0.05*xsize
ioff= where(abs(ewave-4364.) lt 3., noff)
xoff[ioff]=0.20*xsize
ioff= where(abs(ewave-4341.) lt 3., noff)
xoff[ioff]=-0.05*xsize
yoff[ioff]=size

ewave=ewave*(1.+zans.z)
iwave= where(ewave gt xra[0] and ewave lt xra[1], nwave)
for j=0L, nwave-1L do begin
    i= iwave[j]
    inear= where(wave gt ewave[i]-100. and wave lt ewave[i]+100 and $
                 invvar gt 0 and abs(wave-5577.) gt 4., nnear)
    if(nnear gt 0) then begin
        eval= max(sflux[inear])
        djs_oplot, [ewave[i], ewave[i]], $
          eval[0]*1.05+[0.,size]+yoff[i], color='blue', th=3
        djs_xyouts, ewave[i]+xoff[i], eval[0]*1.05+yoff[i]+size*1.2, align=0.5, $
          ename[i], noclip=0, charsize=sscale, charthick=2.
    endif
endfor

yoff= replicate(-size*0.2, n_elements(awave))
xoff= replicate(0., n_elements(awave))

ioff= where(abs(awave-3938.) lt 3., noff)
xoff[ioff]=-0.07*xsize
ioff= where(abs(awave-3968.) lt 3., noff)
xoff[ioff]=0.07*xsize

awave=awave*(1.+zans.z)
iwave= where(awave gt xra[0]+100. and awave lt xra[1]-100., nwave)
for j=0L, nwave-1L do begin
    i= iwave[j]
    inear= where(wave gt awave[i]-100. and wave lt awave[i]+100 and $
                 invvar gt 0 and abs(wave-5577.) gt 4., nnear)
    if(nnear gt 0) then begin
        aval= min(sflux[inear])
        djs_oplot, [awave[i], awave[i]], $
          aval[0]*0.95-[0.,size]+yoff[i], color='red', th=3
        djs_xyouts, awave[i]+xoff[i], aval[0]*0.95+yoff[i]-size*1.5, align=0.5, $
          aname[i], noclip=0, charsize=sscale, charthick=2.
    endif
endfor

device,/close
if(keyword_set(bangp)) then !P=bangp
if(keyword_set(bangx)) then !X=bangx
if(keyword_set(bangy)) then !Y=bangy
set_plot,'x'

if(not keyword_set(silent)) then splog, 'Converting to JPG image.'
spawn, ['convert',outbase+'.ps', '-quality', '100', outbase+'.jpg'], /nosh
spawn, ['convert',outbase+'.jpg', outbase+'.gif'], /nosh
spawn, ['convert',outbase+'.gif', outbase+'.png'], /nosh
if(not keyword_set(silent)) then splog, 'Making thumbnail image.'
spawn, ['convert',outbase+'.png', '-scale', '180x130', outbase+'.thumb.png'], $
       /nosh
rmfile, outbase+'.ps'
rmfile, outbase+'.jpg'
rmfile, outbase+'.gif'
if(not keyword_set(silent)) then splog, 'Done.'

end
