; since eBOSS the platemap has a mismatched fiber problem; this program aims to 
; find these mismatches
; epochs affected: since 57038 (ep33)

pro corr_mis_fiber, mjd_start, mjd_end

file = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
target = mrdfits(file, 1) ; the target fibermap
nobj = 1000L

nep = n_elements(target[0].mjd)
mjdall = target[0].mjd
ind1 = (where(mjdall eq mjd_start))[0] & ind2 = (where(mjdall eq mjd_end))[0]

ncheck = ind2 - ind1 + 1L ; # of eBOSS epochs to check

; note this starts from mjd_start
result = {rmid:0L, plate:lonarr(ncheck), mjd:lonarr(ncheck), old_fiber:lonarr(ncheck), $
   mis_flag:lonarr(ncheck), corr_fiber:lonarr(ncheck)}
result.mjd = (target[0].mjd)[ind1:ind2]
result.plate = (target[0].plate)[ind1:ind2]
result = replicate(result, nobj)

for iobj=0,nobj-1 do begin

  redchi2_arr = dblarr(nobj) + 999L
  mis_flag_arr = lonarr(ncheck)
  corr_fiber_arr = lonarr(ncheck)
  old_fiber_arr = lonarr(ncheck)

  ; readin the reference spectrum; do not specify calibdir
  rm_readspec, 0, iobj+1, $
        mjd = 56837, wave=wave_ref, flux=flux_ref, invvar=ivar_ref
  ; cut the reference wave window to ensure overlap with the comparison spectrum
  wave_ref = wave_ref[150:4000]
  flux_ref = flux_ref[150:4000]
  ivar_ref = ivar_ref[150:4000]
  npix = n_elements(wave_ref)

  print, 'redchi2(old,new), plate, mjd, fiber(old,new)'
  for iep=ind1, ind2 do begin

      ; read in all fibers in this epoch
      rm_readspec, (target[iobj].plate)[iep],  $
           mjd=(target[iobj].mjd)[iep], wave=wave, flux=flux, invvar=ivar

      ; now compare the flux with the reference flux - note the wavegrid is different
      ind = where(ivar_ref gt 1d-5)
      minwave = min(wave_ref[ind]) & maxwave = max(wave_ref[ind] )

      ; find the spectrum that has the minimum redchi2 compared with the reference spectrum
      for j=0, nobj-1 do begin 

         flux_use = interpol(flux[*,j], wave[*,j], wave_ref)
         ivar_use = interpol(ivar[*,j], wave[*,j], wave_ref)

         ; make sure this spectrum is not corrupted
         ind_good = where(ivar_use gt 1d-5, ngood)
         if ngood gt 0.9*npix then begin
            ind_hi = where(flux_ref gt min(flux_ref), n_hi)
            ;scale = median(flux_ref)/median(flux_use)
            ;flux_use = flux_use*scale
            ;ivar_use = ivar_use/scale^2
            chi2 = total((flux_use[ind_hi] - flux_ref[ind_hi])^2*ivar_use[ind_hi], /double)
            if chi2 gt 0 then redchi2_arr[j] = chi2/double(n_hi)
         endif
      endfor

      ; check if it is the same fiber id
      old_fiber = (target[iobj].fiberid)[iep]
      old_fiber_arr[iep - ind1] = old_fiber
      if redchi2_arr[old_fiber - 1] gt 15 then begin
         mis_flag_arr[iep - ind1] = 1L
         ttt = min(redchi2_arr, ind_min)
         new_fiber = ind_min+1
         print, redchi2_arr[old_fiber - 1], ttt, (target[iobj].plate)[iep], $
               (target[iobj].mjd)[iep], old_fiber, new_fiber
         corr_fiber_arr[iep - ind1] = new_fiber
      endif
  endfor 


  result[iobj].mis_flag = mis_flag_arr
  result[iobj].corr_fiber = corr_fiber_arr
  result[iobj].old_fiber = old_fiber_arr

  splog, 'Finished RMID: ', iobj
  ; print, result[iobj].mis_flag
  ; message, 'stop'
  ; pause

endfor

; NOTE that not all these findings are mismatches!! Visual inspection still needed!!
outfile = '/home/yshen/products/Linux/idlrm/etc/corr_fiber_eboss_' + string(mjd_start,format='(i5.5)') + $
  '-' + string(mjd_end, format='(i5.5)') + '.fits'
mwrfits, result, outfile, /create


end


pro output_corr_fiber_list, infile, outfile

if ~keyword_set(infile) then file = '/home/yshen/products/Linux/idlrm/etc/corr_fiber_eboss_57038.fits' $
  else file = infile
result = mrdfits(file,1)

if ~keyword_set(outfile) then outfile = '/home/yshen/products/Linux/idlrm/etc/corr_fiber_eboss_57038'
openw, lun, outfile, /get_lun
printf, lun, '# RMID  Plate MJD  old_fiber  new_fiber'
fmt = '(i3.3, " ", i4.4, " ", i5.5, " ", i4.4, " ", i4.4)'

for i=0L, 999 do begin
  plate = result[i].plate & mjd = result[i].mjd
  o_fiber = result[i].old_fiber & n_fiber = result[i].corr_fiber
  ind = where(n_fiber gt 0, nnn)
  if nnn gt 0 then begin
     plate = plate[ind] & mjd = mjd[ind] & o_fiber = o_fiber[ind]
     n_fiber = n_fiber[ind]
     for jj = 0, nnn-1 do begin
        printf, lun, format=fmt, i, plate[jj], mjd[jj], o_fiber[jj], n_fiber[jj]
     endfor
  endif

endfor

close, lun
free_lun, lun


end

; make diag plots for plate 7339-57518 to pin down the last few mismatches
pro plot_57518_check

file = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
target = mrdfits(file, 1) ; the target fibermap
rmid=lindgen(1000L)
nnn = n_elements(rmid)

plate=target.plate & fiber=target.fiberid & mjd=target.mjd
ind=where(mjd[*,0] eq 57518)
plate=plate[ind,*] & fiber=fiber[ind,*] & mjd=mjd[ind,*]
plate=reform(plate) & fiber=reform(fiber) & mjd=reform(mjd)
count=lonarr(1000)
for i=0,999L do begin & ind=where(fiber eq fiber[i], nn) & count[i]=nn & endfor

; find all the missing fibers
for i=0,999L do begin
  ind=where(fiber eq i+1, nn)
  if nn eq 0 then begin
     if n_elements(missing_fiber) eq 0 then missing_fiber=i+1 else $
       missing_fiber=[missing_fiber,i+1]
  endif
endfor
figfile='/data3/yshen/ftp/sdssrm/collab/bossredux/fiber_mismatch_eboss/57518_tmp/dup_fiber_spec.ps'
begplot, name=figfile,/color
ind = where(count gt 1)
rmid_use = rmid[ind] & fiber_use = fiber[ind]
ind = sort(fiber_use)
rmid_use = rmid_use[ind] & fiber_use = fiber_use[ind]
for i=0, n_elements(rmid_use) - 1 do begin

  if (i mod 2) eq 0 then begin 
    noerase = 0 & pos = [0.12,0.55, 0.96, 0.92]
  endif else begin
    noerase = 1 & pos = [0.12, 0.08, 0.96, 0.45]
  endelse
  rm_readspec,0,rmid_use[i]+1,mjd=56837,wave=wave0,flux=flux0
  rm_readspec,7339,fiber_use[i],mjd=57518,wave=wave,flux=flux
  title='RMID'+string(rmid_use[i],format='(i3.3)')+', ' + '7339-57518-'+string(fiber_use[i],format='(i4.4)')
  plot, wave0,flux0/median(flux0), xrange=[3500,10500.],/xsty, thick=3,title=title, $
    pos=pos,noerase=noerase,/nodata, yrange=[-2,10]
  oplot, wave, flux/median(flux), color=cgcolor('red'),thick=3
  oplot, wave0,flux0/median(flux0)
  flux_sm = flux/median(flux)
  oplot, wave, median(flux_sm,20), color=cgcolor('cyan'),thick=3
endfor

endplot

; plot the spectra of all missing fibers
figfile='/data3/yshen/ftp/sdssrm/collab/bossredux/fiber_mismatch_eboss/57518_tmp/missing_fiber_spec.ps'
begplot, name=figfile,/color, /landscape

for i=0,n_elements(missing_fiber)-1 do begin
  title = '7339-57518-' + string(missing_fiber[i],format='(i4.4)')
  rm_readspec, 7339,missing_fiber[i],mjd=57518, wave=wave,flux=flux
  plot, wave, flux/median(flux), yrange=[-2,10], title=title, xrange=[3500,10500.],/xsty, thick=3
  flux_sm = flux/median(flux)
  oplot, wave, median(flux_sm, 15), thick=5,color=cgcolor('red')
endfor

endplot
cgfixps, figfile

end


; take a close look at the mismatched fibers
pro check_corr_fiber_close, rmid, plate, mjd, o_fiber,start_obj=start_obj, plotnextfiber=plotnextfiber,run2d=run2d

file = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
target = mrdfits(file, 1) ; the target fibermap
nnn = n_elements(rmid)

if ~keyword_set(start_obj) then start_obj=1L
i=start_obj - 1L

if ~keyword_set(outfile) then outfile = 'output.marked'
openw, lun, outfile, /get_lun
printf, lun, '# RMID  Plate MJD  old_fiber new_fiber'
fmt = '(i3.3, " ", i4.4, " ", i5.5, " ", i4.4, " ", a0)'


while i le nnn - 1L do begin

  rm_readspec, 0, rmid[i]+1, mjd=56837, wave=wave0, flux=flux0
  flux0 = flux0/median(flux0)
  title = 'RIMD'+string(rmid[i],format='(i3.3)')+',' + string(plate[i],format='(i4.4)')+ $
    '-' + string(mjd[i], format='(i5.5)')+'-'+string(o_fiber[i],format='(i3.3)')
  plot, wave0, flux0, yrange=[-1,8], title=title,/nodata

  rm_readspec, plate[i], o_fiber[i],mjd=mjd[i],calibdir='/wh_skysub/',wave=wave,flux=flux,run2d=run2d
  flux=flux/median(flux)
  oplot, wave, median(flux, 5), color=cgcolor('blue')

  rm_readspec, plate[i], o_fiber[i]+1,mjd=mjd[i],calibdir='/wh_skysub/',wave=wave,flux=flux,run2d=run2d
  flux=flux/median(flux)
  if keyword_set(plotnextfiber) then begin
    oplot, wave, median(flux, 5), color=cgcolor('red')
    xyouts, 0.5, 0.85, string(o_fiber[i]+1,format='(i3.3)'),color=cgcolor('red'), /norm
  endif
  oplot, wave0, flux0

  i=i+1  ; move to the next object
  mark = 'NO'
  mark_input = ''
  read, mark_input, prompt='QUIT? [NO] '

  if strmatch(strmid(mark_input,0,1),'#') then begin  ; contral commands rather than comments
    if strmatch(strmid(mark_input,1,1), 'z') then begin
       new_z = double(strmid(mark_input,2))
    endif else new_z = -1
    i = i-2  ; move back one more object
  endif else begin
    if strlen(mark_input) gt 1L then mark = mark_input
    if strmatch(mark_input, 'y') then mark = string(o_fiber[i-1] +1, format='(i4.4)')
    if strmatch(mark_input, 'e') or strmatch(mark_input, 'q') then begin
       close, lun
       free_lun, lun
       return
    endif

    ; only write out if the mark is not "NO"
    if not strmatch(mark,'NO') then printf, lun, format = fmt, rmid[i-1], plate[i-1], mjd[i-1], $
      o_fiber[i-1], mark
  endelse
endwhile

close, lun
free_lun, lun

end


; look through all potential mismatched fibers and flag the ones that are indeed mismatches
pro check_corr_fiber, infile=infile, outfile=outfile,start_obj=start_obj,plate=plate,mjd=mjd,o_fiber=o_fiber, $
   rmid_all=rmid_all, output=output

if ~keyword_set(infile) then file = '/home/yshen/products/Linux/idlrm/etc/corr_fiber_eboss_57038-57576' $
  else file = infile
if n_elements(rmid_all) eq 0 then begin
  readcol, file, format='l, l,l,l', rmid_all, plate, mjd, o_fiber
endif
nnn=n_elements(rmid_all)
if ~keyword_set(start_obj) then start_obj=1L


if keyword_set(output) then begin
  if ~keyword_set(outfile) then $
    outfile = file + '.marked'
  openw, lun, outfile, /get_lun
  printf, lun, '# RMID  Plate MJD  old_fiber new_fiber'
endif
fmt = '(i3.3, " ", i4.4, " ", i5.5, " ", i4.4, " ", a0)'
charsize=1.5

i=start_obj - 1L
while i le nnn - 1L do begin
  
  rmid = rmid_all[i]

  ; get the ref spec
  rm_readspec, 0, rmid+1, mjd=56837, wave=wave0, flux=flux0
  ; get the new spec
  rm_readspec, plate[i],o_fiber[i], mjd=mjd[i], wave=wave1, flux=flux1
  if o_fiber[i]+1 le 1000 then begin
     rm_readspec, plate[i],o_fiber[i]+1, mjd=mjd[i], wave=wave2, flux=flux2
  endif else begin ; no fiber spectrum is read
     wave2 = wave1 & flux2 = wave2*0.D
  endelse

  title = string(i,format='(i0)') + ' RMID '+string(rmid,format='(i3.3)') + ' ' + string(plate[i], format='(i4.4)') + $
     ' ' + string(mjd[i], format='(i5.5)')

  sflux = median(flux0[50:4000], 15)
  sflux1 = median(flux1[50:4000], 15)
  sflux2 = median(flux2[50:4000], 15)
  ;yrange = [min([sflux,sflux1,sflux2]), max([sflux,sflux1,sflux2])]
  yrange = [-1, max([sflux,sflux1,sflux2])]

  plot, wave0, flux0, xrange=[3600, 1d4],/xsty, title=title,charsize=charsize, yrange=yrange,/ysty
  oplot, wave1, smooth(flux1,5), color=cgcolor('red')
  oplot, wave2, smooth(flux2,5), color=cgcolor('green')
  xyouts, 0.2, 0.9, /norm, 'old_fiber='+string(o_fiber[i],format='(i4.4)'),color=cgcolor('red'),charsize=charsize
  xyouts, 0.6, 0.9, /norm, 'fiber+1='+string(o_fiber[i]+1,format='(i4.4)'),color=cgcolor('green'), charsize=charsize
 
  i=i+1  ; move to the next object

  mark = 'NO'
  mark_input = ''
  read, mark_input, prompt='FLAG? [NO] '

  if strmatch(strmid(mark_input,0,1),'#') then begin  ; contral commands rather than comments
    if strmatch(strmid(mark_input,1,1), 'z') then begin
       new_z = double(strmid(mark_input,2))
    endif else new_z = -1
    i = i-2  ; move back one more object
  endif else begin
    if strlen(mark_input) gt 1L then mark = mark_input
    if strmatch(mark_input, 'y') then mark = string(o_fiber[i-1] +1, format='(i4.4)')
    if strmatch(mark_input, 'e') or strmatch(mark_input, 'q') then begin
       close, lun
       free_lun, lun
       return
    endif

    ; only write out if the mark is not "NO"
    if keyword_set(output) then begin
      if not strmatch(mark,'NO') then printf, lun, format = fmt, rmid, plate[i-1], mjd[i-1], $
        o_fiber[i-1], mark
    endif
  endelse

endwhile

if keyword_set(output) then begin
  close, lun
  free_lun, lun
endif

end


pro plot_misfiber
; plot the possible mismatched fibers

file = '/home/yshen/products/Linux/idlrm/etc/target_fibermap.fits'
target = mrdfits(file, 1) ; the target fibermap
nobj = 1000L

file = '/home/yshen/products/Linux/idlrm/etc/corr_fiber_eboss_57038.txt'
readcol, file, format='l,l,l,l', rmid, plate, mjd, fiber


file = '/home/yshen/products/Linux/idlrm/etc/corr_fiber_eboss_57038.txt.flagged'
readcol,file, format='x,x,a', flag

; ind_mis = [125, 197, 209, 214, 230, 247, 1596, 1599, 1600, 1618, 1625] - 1L
; ind_mis = where(plate eq 7339 and mjd eq 57518)
ind_mis = where(plate eq 7338 and mjd eq 57038) ;; looks all OK
ind_mis = where(plate eq 7338 and mjd eq 57135);; a few affected
;ind_mis = where(plate eq 7340 and mjd eq 57196) ; nothing
;ind_mis = where(plate eq 7340 and mjd eq 57576) ; nothing
;ind_mis = where(plate eq 7338 and mjd eq 57166) ; nothing
ind_mis = where(plate eq 7339 and mjd eq 57481) ; nothing
;ind_mis = where(plate eq 7339 and mjd eq 57510) ; nothing
;ind_mis = where(plate eq 7339 and mjd eq 57082) ; nothing
;ind_mis = where(plate eq 7339 and mjd eq 57544) ; nothing
;ind_mis = where(plate eq 7339 and mjd eq 57492) ; nothing
 

rmid = rmid[ind_mis] & plate = plate[ind_mis] & mjd = mjd[ind_mis] & fiber = fiber[ind_mis]
nnn = n_elements(rmid)

figfile='/home/yshen/products/Linux/idlrm/etc/plot_misfiber.ps'
begplot, name=figfile, /color, /landscape

colors = cgcolor(['black', 'red', 'green', 'cyan'])
for i=0, nnn - 1 do begin

  str = string(plate[i],format='(i4.4)')+'-'+string(mjd[i],format='(i5.5)')+'-'+string(fiber[i],format='(i4.4)')
  items = ['RMID '+string(rmid[i],format='(i3.3)'), str, 'Fiber+1','Fiber-1']
  rm_readspec, 0, rmid[i] + 1, mjd=56837, wave=wave,flux=flux, calibdir='wh_skysub/'
  plot, wave, median(flux, 15), xtitle='Wavelength', ytitle = 'Flux', color=colors[0]
  rm_readspec, plate[i], fiber[i], mjd = mjd[i], wave=wave,flux=flux, calibdir='wh_skysub/'
  oplot, wave, median(flux, 15), color=colors[1]
  if fiber[i]+1 le 1000 then $
    rm_readspec, plate[i], fiber[i]+1, mjd = mjd[i], wave=wave,flux=flux, calibdir='wh_skysub/'
  oplot, wave, median(flux, 15), color=colors[2]
  rm_readspec, plate[i], fiber[i]-1, mjd = mjd[i], wave=wave,flux=flux, calibdir='wh_skysub/'
  oplot, wave, median(flux, 15), color=colors[3]

  legend, items, color=colors, pos=[0.7, 0.92], box=0, /norm, textcolor=colors

endfor

endplot
cgfixps, figfile


end

; check the 2018 spectroscopy
pro check_2018


end
