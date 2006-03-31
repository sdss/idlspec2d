
   f = findfile('disp_bounds.fits')
   if f[0] EQ '' then begin

   public = ['EDR','DR1','DR2','DR3','DR4','DR5']
   platelist, plist=plist
   qkeep = bytarr(n_elements(plist))
   for i=0, n_elements(public)-1 do $
    qkeep = qkeep OR strmatch(plist.public,public[i])
   qkeep = qkeep AND strmatch(plist.statuscombine,'Done*')
   plist = plist[where(qkeep, nplate)]

   nplate = n_elements(plist)


   bluepix = 210
   blue_loglam = 3.5800 + dindgen(bluepix)*1d-3
   redpix = 190
   red_loglam = 3.7700 + dindgen(redpix)*1d-3

   plstr = string(plist.plate,format='(i4.4)')
   iplate =0
   spfiles = findfile('/raid/spectro/2d_v5_1/'+plstr[iplate]+'/spFr*gz')
   r1 = where(strpos(spfiles, 'r1') GT 0)
   r1w = mrdfits(spfiles[r1[0]],3)
   r1d = mrdfits(spfiles[r1[0]],4)

   stop
   wset = { func    :    'legendre' , xmin  : 0.0d, xmax :  2047.0d, $
       coeff   :    dblarr(5, 16) }
   dset = { func    :    'legendre' , xmin  : 0.0d, xmax :  2047.0d, $
       coeff   :    dblarr(4, 16) }

   fibers = lindgen(16)*20 + 10
   wset_all = replicate(wset, 4, nplate)
   dset_all = replicate(dset, 4, nplate)

   cameras = ['b1', 'b2', 'r1', 'r2']

   for iplate =0, nplate-1 do begin
      print, iplate
      spfiles = findfile('/raid/spectro/2d_v5_1/'+plstr[iplate]+'/spFr*gz')
      for ic=0,3 do begin
        first = where(strpos(spfiles, cameras[ic]) GT 0)
        ww = mrdfits(spfiles[first[0]],3, /silent)
        wset_all[ic,iplate].coeff = ww.coeff[*,fibers]
        tt = mrdfits(spfiles[first[0]],4, /silent)
        if ic LE 1 then dset_all[ic,iplate].coeff[0:2,*] = tt.coeff[*,fibers] $
        else dset_all[ic,iplate].coeff = tt.coeff[*,fibers]
      endfor
   endfor

   

   pix_arr  = fltarr(bluepix, 16, nplate)
   disp_arr = fltarr(bluepix, 16, nplate)
   del_arr  = fltarr(bluepix, 16, nplate)
  
   for iplate =0, nplate-1 do begin & $
     traceset2xy, wset_all[0,iplate], p, l & $
     traceset2xy, dset_all[0,iplate], p, d & $
     for i=0,15 do begin & $
        pix_arr[*,i,iplate] = interpol(p[*,i], l[*,i], blue_loglam) & $
        disp_arr[*,i,iplate] = interpolate(d[*,i], pix_arr[*,i,iplate]) & $
        del_arr[*,i,iplate] = abs(interpol(p[*,i], l[*,i], $
                            blue_loglam + 0.0001) - pix_arr[*,i,iplate]) & $
     endfor & $
   endfor

   b1_arr = disp_arr / del_arr *69.02 * 2.3548

   for iplate =0, nplate-1 do begin & $
     traceset2xy, wset_all[1,iplate], p, l & $
     traceset2xy, dset_all[1,iplate], p, d & $
     for i=0,15 do begin & $
        pix_arr[*,i,iplate] = interpol(p[*,i], $
                                       l[*,i], blue_loglam) & $
        disp_arr[*,i,iplate] = interpolate(d[*,i], pix_arr[*,i,iplate]) & $
        del_arr[*,i,iplate] = abs(interpol(p[*,i], l[*,i], $
                            blue_loglam + 0.0001) - pix_arr[*,i,iplate]) & $
     endfor & $
   endfor

   b2_arr = disp_arr / del_arr *69.02 * 2.3548

   pix_arr  = fltarr(redpix, 16, nplate)
   disp_arr = fltarr(redpix, 16, nplate)
   del_arr  = fltarr(redpix, 16, nplate)
  
   for iplate =0, nplate-1 do begin & $
     traceset2xy, wset_all[2,iplate], p, l & $
     traceset2xy, dset_all[2,iplate], p, d & $
     for i=0,15 do begin & $
        pix_arr[*,i,iplate] = interpol(p[*,i], $
                                       l[*,i], red_loglam) & $
        disp_arr[*,i,iplate] = interpolate(d[*,i], pix_arr[*,i,iplate]) & $
        del_arr[*,i,iplate] = abs(interpol(p[*,i], l[*,i], $
                            red_loglam + 0.0001) - pix_arr[*,i,iplate]) & $
     endfor & $
   endfor

   r1_arr = disp_arr / del_arr *69.02 * 2.3548

   for iplate =0, nplate-1 do begin & $
     traceset2xy, wset_all[3,iplate], p, l & $
     traceset2xy, dset_all[3,iplate], p, d & $
     for i=0,15 do begin & $
        pix_arr[*,i,iplate] = interpol(p[*,i], $
                                       l[*,i], red_loglam) & $
        disp_arr[*,i,iplate] = interpolate(d[*,i], pix_arr[*,i,iplate]) & $
        del_arr[*,i,iplate] = abs(interpol(p[*,i], l[*,i], $
                            red_loglam + 0.0001) - pix_arr[*,i,iplate]) & $
     endfor & $
   endfor

   r2_arr = disp_arr / del_arr *69.02 * 2.3548

   b1_bounds = fltarr(bluepix, 3)
   b2_bounds = fltarr(bluepix, 3)
   r1_bounds = fltarr(redpix, 3)
   r2_bounds = fltarr(redpix, 3)

   p10 = long(nplate * 16 * 0.10)
   p50 = long(nplate * 16 * 0.50)
   p90 = long(nplate * 16 * 0.90)
   for j=0,bluepix-1 do begin & $
      t = (b1_arr[j, *, *])[*] & s = sort(t) & $
      b1_bounds[j,*] = t[s[[p10,p50,p90]]] & $
      t = (b2_arr[j, *, *])[*] & s = sort(t) & $
      b2_bounds[j,*] = t[s[[p10,p50,p90]]] & $
   endfor
   for j=0,redpix-1 do begin & $
      t = (r1_arr[j, *, *])[*] & s = sort(t) & $
      r1_bounds[j,*] = t[s[[p10,p50,p90]]] & $
      t = (r2_arr[j, *, *])[*] & s = sort(t) & $
      r2_bounds[j,*] = t[s[[p10,p50,p90]]] & $
   endfor

     mwrfits, b1_bounds, 'disp_bounds.fits', /create
     mwrfits, b2_bounds, 'disp_bounds.fits'
     mwrfits, r1_bounds, 'disp_bounds.fits'
     mwrfits, r2_bounds, 'disp_bounds.fits'

   endif else begin
     b1_bounds = mrdfits('disp_bounds.fits',0)
     b2_bounds = mrdfits('disp_bounds.fits',1)
     r1_bounds = mrdfits('disp_bounds.fits',2)
     r2_bounds = mrdfits('disp_bounds.fits',3)

     bluepix = 210
     blue_loglam = 3.5800 + dindgen(bluepix)*1d-3
     redpix = 190
     red_loglam = 3.7700 + dindgen(redpix)*1d-3
   endelse

   
   dfpsplot, 'plot_disp.ps', /color, /square
   !p.thick = 6
   djs_plot, 10^blue_loglam, b1_bounds[*,1], /nodata, xr=[3800,9200], /xs, yr=[90,220], /ys, $
     charthick=4., chars=1.5, thick=6, xtitle='Wavelength (\AA)', $
     ytitle='FWHM (km/s)', xthick=6, ythick=6

   s = lindgen(11) * 20 + 3
   oplot, 10^blue_loglam, b1_bounds[*,1], color=djs_icolor('cyan'), thick=8
   errplot, 10^blue_loglam[s], b1_bounds[s,0], b1_bounds[s,2], color=djs_icolor('cyan')
   s = lindgen(11) * 20 + 7
   oplot, 10^blue_loglam, b2_bounds[*,1], color=djs_icolor('blue'), thick=8, lines=2
   errplot, 10^blue_loglam[s], b2_bounds[s,0], b2_bounds[s,2], color=djs_icolor('blue'), lines=2

   s = lindgen(10) * 20 + 3 
   oplot, 10^red_loglam, r1_bounds[*,1], color=djs_icolor('orange'), thick=8
   errplot, 10^red_loglam[s], r1_bounds[s,0], r1_bounds[s,2], color=djs_icolor('orange') 
   s = lindgen(10) * 20 + 7 
   oplot, 10^red_loglam, r2_bounds[*,1], color=djs_icolor('red'), thick=8, lines=2
   errplot, 10^red_loglam[s], r2_bounds[s,0], r2_bounds[s,2], color=djs_icolor('red') , lines=2
   dfpsclose
   
   
      
 end
