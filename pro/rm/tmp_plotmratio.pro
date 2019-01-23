;+
; NAME: 
;   tmp_plotmratio
; PURPOSE:
;   plot the mratio vectors for all the good std in an exposure 
; --------------------------------------------------------------------
pro tmp_plotmratio, objname, calibdir = calibdir, linear=linear, _extra=extra

    if not keyword_set(calibdir) then $
        calibdir = '/data3/quasar/yshen/spectro/bossredux/v5_6_0/7338/recalib/'

    file = calibdir + objname
    kindx = mrdfits(file, 2, /silent)

    ; good stars used in fluxing
    ind_good = where(kindx.qgood eq 1, ngood)
    colors = fsc_color(['opposite', 'red', 'green', 'magenta', 'cyan'])
    ncolor = n_elements(colors)

    if strmatch(objname, '*-b*') then begin
      xrange=[3.55, 3.8] & yrange = [0,20]
      if keyword_set(linear) then xrange = [3500, 6500]
    endif
    if strmatch(objname, '*-r*') then begin
      xrange=[3.7, 4.05] & yrange = [0,40]
      if keyword_set(linear) then xrange = [5000, 11000]
    endif

    xpos = xrange[0] + 0.05*(xrange[1]-xrange[0])
    ypos = yrange[1] - 0.1*(yrange[1]-yrange[0])
    dypos = 0.05*(yrange[1]-yrange[0])
    dxpos = 0.05*(xrange[1]-xrange[0])

    plot, [0],[0], xrange=xrange, yrange=yrange, xtitle=textoidl('log \lambda') $
      , ytitle = 'mratio/mratio0', /nodata, /xsty
   

    for i=0L, ngood - 1L do begin
        if not keyword_set(linear) then $
          oplot, kindx[ind_good[i]].loglam, kindx[ind_good[i]].mratio $
             , color=colors[i mod ncolor], _extra=extra else $
          oplot, 10.D^kindx[ind_good[i]].loglam, kindx[ind_good[i]].mratio $
             , color=colors[i mod ncolor], _extra=extra
        xyouts, xpos, ypos, string(kindx[ind_good[i]].fiberid,format='(i0)') $
           , color=colors[i mod ncolor], charsize=1.5
        xpos = xpos + dxpos

    endfor

end
