; This reads the save save written by RLINE_SAMPLE.
; Set PLOTNUM to 0, 1, or 2 for plot number
pro rline_plotsome, plotnum

   if (NOT keyword_set(plotnum)) then plotnum = 0
;   nsmooth = 3

   restore, 'goodpks.ss'

   junk = {pinfo, iselect: 0L, string1: '', string2: ''}
   pinfo = replicate( {pinfo}, 4 )

   if (plotnum EQ 0) then begin
      pinfo.iselect= [ 15        , 213       , 226       , 140        ]
      pinfo.string1= [ 'H\alpha' , 'H\alpha' , 'H\alpha' , 'H\alpha'  ]
      pinfo.string2= [ 'z=0.072' , 'z=0.136' , 'z=0.039' , 'z=0.018'  ]
   endif else if (plotnum EQ 1) then begin
      pinfo.iselect= [ 106       , 37        , 35        , 276        ]
      pinfo.string1= [ 'O II'    , 'O II'    , 'O II'    , 'O II'     ]
      pinfo.string2= [ 'z=0.272' , 'z=0.441' , 'z=0.806' , 'z=0.156'  ]
   endif else if (plotnum EQ 2) then begin
      pinfo.iselect= [ 179          , 93           , 143          , 66            ]
      pinfo.string1= [ 'Ly\alpha ??', 'Ly\alpha ??', 'Ly\alpha ??', 'Ly\alpha ??' ]
      pinfo.string2= [ 'z=4.22'     , 'z=2.17'     , 'z=2.24'     , 'z=5.13'      ]
   endif

   dfpsplot, 'rline-panel-'+strtrim(string(plotnum),2)+'.ps', /color

   iselect = pinfo.iselect
   psize = 1.5
   pthick = 2
pthick=1

   ; Interpolate over bad sky regions
   invvar = skymask(invvar, andmask)
   fluxinterp = djs_maskinterp(flux, invvar EQ 0, /const, iaxis=0)

   nplot = n_elements(iselect)
   yoffset = 0.90 / nplot
;   ypos = [0.10, 0.10+yoffset] ; Do this to not separate each spectrum in Y
   ypos = [0.10, 0.06+yoffset] ; Do this to separate each spectrum in Y
   for iplot=0, nplot-1 do begin
      if (iplot EQ 0) then begin
         xtitle='Wavelength [Ang]'
         noerase = 0
         xtickname = ''
      endif else begin
         xtitle=''
         noerase = 1
         xtickname = replicate(' ',30)
      endelse
      if (iplot EQ nplot-1) then begin
         title1 = 'FULL SDSS SPECTRUM'
         title2 = 'ROGUE LINES'
      endif
      ytitle = 'Flux'

      ; Info on the rogue line
      thispeak = goodpks[iselect[iplot]]
      thiswave = 10.^thispeak.xloglam

      ; Info on the foreground/background galaxy
      thisz = thispeak.zans.z

      iw = where(wave[*,iselect[iplot]] NE 0)
      xplot = wave[iw,iselect[iplot]]
      yplot1 = fluxinterp[iw,iselect[iplot]]
      yplot2 = flux[iw,iselect[iplot]]
      synplot = synflux[iw,iselect[iplot]]
      errplot = flerr[iw,iselect[iplot]]

      ; Smooth...
      if (keyword_set(nsmooth)) then begin
         yplot1 = smooth(yplot1,nsmooth)
         synplot = smooth(synplot,nsmooth)
         errplot = smooth(errplot,nsmooth)
      endif

      ; Plot the whole spectrum

      xrange = [3700,9200]
      yrange = 1.3 * minmax(smooth(yplot1,5))
      yrange[0] = yrange[0] < 0
      yrange[0] = yrange[0] > (-0.1 * yrange[1])

      plot, xplot, yplot1, /nodata, noerase=noerase, $
       position=[0.10,ypos[0],0.65,ypos[1]], $
       xtitle=xtitle, ytitle=ytitle, $
       xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
       charsize=psize, charthick=pthick, $
       xtickname=xtickname, title=title1
      djs_oplot, xplot, yplot1, thick=pthick
      djs_oplot, xplot, synplot, color='blue', thick=pthick
      djs_oplot, xplot, errplot, color='red', thick=pthick
      djs_xyouts, 0.9*xrange[0]+0.1*xrange[1], 0.2*yrange[0]+0.8*yrange[1], $
       'Galaxy @ z='+string(thisz,format='(f5.2)'), $
       charsize=psize, charthick=pthick

      ; Indicate the rogue em. line

      yarr = [yrange[0], 0.75*yrange[0]+0.25*yrange[1]]
      djs_arrow, thiswave, yarr[0], thiswave, yarr[1], $
       color='green', thick=pthick, /data

      ; Plot the spectrum around the rogue em. line

      xrange = thiswave+[-100,100]
      plot, xplot, yplot2, /nodata, /noerase, $
       position=[0.65,ypos[0],0.97,ypos[1]], $
       xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
       charsize=psize, charthick=pthick, $
       ytickname=replicate(' ',30), title=title2
      djs_oplot, xplot, yplot2, thick=pthick, ps=10
      djs_oplot, xplot, synplot, color='blue', thick=pthick

      djs_xyouts, 0.35*xrange[0]+0.65*xrange[1], 0.2*yrange[0]+0.8*yrange[1], $
       pinfo[iplot].string1, charsize=psize, charthick=pthick, color='green'
      djs_xyouts, 0.35*xrange[0]+0.65*xrange[1], 0.35*yrange[0]+0.65*yrange[1], $
       pinfo[iplot].string2, charsize=psize, charthick=pthick, color='green'

      ; Indicate the rogue em. line

      djs_arrow, thiswave, yarr[0], thiswave, yarr[1], $
       color='green', thick=pthick, /data

      ypos = ypos + yoffset
   endfor

   dfpsclose

   return
end
