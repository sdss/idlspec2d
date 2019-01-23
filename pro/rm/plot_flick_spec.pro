; make a 2D image to show the "flicker" of all 849 quasars

pro plot_flick_spec, epoch = epoch, psplot=psplot


file = getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
target=mrdfits(file,1)
target=target[0:848]
ind=sort(target.zfinal)
rmid=ind
nobj = n_elements(rmid)
target = target[ind]
zz = target.zfinal

if n_elements(epoch) eq 0 then epoch = indgen(32) + 1L
nep = n_elements(epoch)

rm_readspec, 0, 1, mjd=56837,wave=wave
npix=n_elements(wave)

prepspec_dir = '/data3/yshen/ftp/sdssrm/collab/prepspec/2014a/'
if n_elements(maxvalue) eq 0 then maxvalue= 1.
if n_elements(minvalue) eq 0 then minvalue= -1.
range=[minvalue, maxvalue]


for ii=0, nep - 1 do begin

  img = dblarr(npix, nobj)
  i = epoch[ii] - 1L

  for j=0, nobj - 1 do begin
    rm_readspec, 0, rmid[j]+1, mjd=56837,wave=wave, flux=flux0
    plate = (target[j].plate)[i]
    mjd = (target[j].mjd)[i]
    fiber = (target[j].fiberid)[i]
    rm_readspec, plate, fiber, mjd=mjd, wave=wave1, flux=flux1, calibdir = 'wh_skysub/'
    ; pt correction
    pt = rm_get_pt(rmid[j], topdir=prepspec_dir, mjd_all=mjd)
    flux1 = flux1/pt[0]
    ; map onto the same wavelength grid
    flux = interpol(flux1, wave1, wave)
    amp = median(flux)/median(flux0)
    if amp gt 1 then scale = 2. else scale = 1./2
    
    img[*, j] = flux/median(flux0)*scale

  endfor

  ; now plot this 2D image for epoch epoch[ii]
  trim_npix=100L
  img=img[trim_npix:*, *] & wave=wave[trim_npix:*]
  img = alog10(img)
  tag = 'Epoch ' + string(epoch[ii],format='(i2.2)') + '(MJD ' + string(mjd,format='(i5.5)') + ')'

  nx=n_elements(wave) & ny=n_elements(zz)

  angstr=textoidl('\AA')
  if keyword_set(psplot) then begin
    if keyword_set(figfile) then fig=figfile else $
    fig='/data3/yshen/ftp/sdssrm/collab/prepspec/2014a/figs/epoch' + string(epoch[ii],format='(i2.2)') + '.ps'
    begplot, name=fig,/color,/landscape
    angstr=string(197B)
  endif

  lines=textoidl(['Ly\alpha', 'SiIV', 'CIV', 'CIII]', 'MgII', 'HeII', 'H\beta', 'H\alpha'])
  linewave=[1215.67, (1396.76 + 1402.06)*0.5, 1549.06,1908.73,2798.75,4687.02, 4862.68,6564.61, $
  4341.68, 4102.89, 3890., 3800., 3588.30]
  linepix=[800,700,600,500, 400, 150, 80, 30, $
    200,250,200,300,300]
  zline=interpol(zz, indgen(ny), linepix)
  linewave_pix=interpol(indgen(nx), wave, linewave*(1+zline))
  nline=n_elements(lines)
  xrange=[0,nx-1] & yrange=[0,ny-1]
  pos=[0.08,0.1,0.92,0.9]
  title='Epoch Spectra'
  charsize=1
  plot,[0],[0], xrange=xrange, yrange=yrange, xsty=5,ysty=5, pos=pos, /norm, /nodata
  cgLoadCT, 55, /reverse
  cgimage, img, maxvalue=maxvalue, minvalue=minvalue, pos=pos,/norm, /noerase
  cgColorbar, Position=[pos[0], pos[3], pos[2], pos[3]+0.05], range=range,/top
  cgLoadCT, 0

  for j=0, nline-1 do begin
    ;print, linewave_pix[j], linepix[j], lines[j]
    xyouts, linewave_pix[j], linepix[j], lines[j], color=cgcolor('blue'),charsize=charsize
  endfor

  ; map the linear wavelength and z grid onto the pixelized image
  zgrid=[0, 1, 2, 3, 4]
  zgrid_pix=interpol(indgen(ny), zz, zgrid)
  wgrid=[4000, 5000, 6000, 7000, 8000, 9000, 10000.]
  wgrid_pix=interpol(indgen(nx), wave, wgrid)

  axis, xaxis=0, xrange=xrange, /xsty, xtitle=textoidl('Observed Wavelength (')+angstr+')', xtickv=wgrid_pix, $
     xtickname=string(wgrid, format='(i0)'), xticks=n_elements(wgrid) - 1
  axis, yaxis=0, yrange=yrange, /ysty, ytitle='Redshift', ytickv=zgrid_pix, $
     ytickname=string(zgrid, format='(i0)'), yticks=n_elements(zgrid) - 1
  axis, yaxis=1, yrange=yrange, /ysty, ytitle='# of Quasars'
  xyouts, pos[0], pos[1]-0.08, tag, /norm,color=cgcolor('red')

  if keyword_set(psplot) then begin
    endplot
    cgfixps, fig
  endif


endfor

end
