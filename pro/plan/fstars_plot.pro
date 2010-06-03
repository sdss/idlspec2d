;+
; NAME:
;   fstars_plot
;
; PURPOSE:
;   Make color-color plots of spectro-photometric stars as compared
;   to the stellar locus and BD+17.
;
; CALLING SEQUENCE:
;   fstars_plot, [ /deredden ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   deredden  - If set, the de-redden the spectroscopic F stars,
;               and the stellar locus stars from run 94.  Do not
;               de-redden the colors of BD+17.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   For finding colors of stars on the stellar locus, this routine
;   reads tsObj files for run 94 rerun 7, which must be on disk.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $BOSS_SPECTRO_REDUX/spAll.fits
;   $PHOTO_DATA/94/7/calibChunks/$CAMCOL/tsObj-000094-$CAMCOL-7-$FIELD.fit
;
; PROCEDURES CALLED:
;   dfpsclose
;   dfpsplot
;   djs_filepath()
;   djs_oplot
;   djs_plot
;   djs_xyouts
;   mrdfits()
;   objc_select()
;   rdss_obj()
;
; REVISION HISTORY:
;   12-Sep-2002  Written by D. Schlegel, Princeton
;------------------------------------------------------------------------------
pro fstars_plot, deredden=deredden

   cname = ['u','g','r','i','z']

   ;----------
   ; Read spAll file, and trim to plates >= 432.

   spfile = djs_filepath('spAll-public.fits', root_dir=getenv('BOSS_SPECTRO_REDUX'))
   spall = mrdfits(spfile, 1)
   indx = where(spall.plate GE 432 AND spall.specprimary EQ 1)
   spall = spall[indx]

   ;----------
   ; Get colors of the F spectro-photo stars

   iredden = where((spall.sectarget AND 2L^1) NE 0)
   isphoto = where((spall.sectarget AND 2L^5) NE 0)

   spcolor = spall.psfcounts[0:3] - spall.psfcounts[1:4]

   ; De-redden these stars, as is done in the target-selection code
   if (keyword_set(deredden)) then begin
      redden = spall.reddening[0:3] - spall.reddening[1:4]
      spcolor = spcolor - redden
      for i=0, 4 do $
       splog, 'Median reddening ' + cname[i] + '-band = ', $
        median(spall[isphoto].reddening[i])
   endif

   ;----------
   ; Get colors of BD+17

   bd17mag = [10.560, 9.626, 9.346, 9.245, 9.232]
   bd17color = bd17mag[0:3] - bd17mag[1:4]

   ;----------
   ; Colors in target selection code (not defined for i-z)

   targetcolor = [0.80, 0.30, 0.10, 99]

   ;----------
   ; Get colors for stars on the stellar locus for some run

   objs = rdss_obj(94, 7, [2,3,4,5], 200+lindgen(100), ftype='tsObj')
   indx = objc_select(objs, ancestry='unique', objtype='star', /trim)
   objs = objs[indx]
   indx = where(objs.psfcounts[2] LT 21)
   objs = objs[indx]

   locuscolor = objs.psfcounts[0:3,*] - objs.psfcounts[1:4,*]

   ; De-redden these stars
   if (keyword_set(deredden)) then begin
      redden = objs.reddening[0:3] - objs.reddening[1:4]
      locuscolor = locuscolor - redden
   endif

   ;----------
   ; Make color-color plots

   dfpsplot, 'fstars.ps', /color, /square
   csize = 1.2

   for iplot=0, 1 do begin
      !p.multi = [0,2,2]

      if (iplot EQ 0) then begin
         indx = isphoto
         starname = 'Spectro-photo stars'
      endif else begin
         indx = isphoto
         starname = 'Reddening stars'
      endelse

      ; Loop over different color-color plots
      for ic1=0, 2 do begin
         ic2 = ic1+1
         xrange = median(spcolor[ic1,iredden]) + [-1,1] * 0.4
         yrange = median(spcolor[ic2,iredden]) + [-1,1] * 0.4

         xtitle='('+cname[ic1]+'-'+cname[ic1+1]+')'
         ytitle='('+cname[ic2]+'-'+cname[ic2+1]+')'
         if (keyword_set(deredden)) then begin
            xtitle = 'De-reddened ' + xtitle
            ytitle = 'De-reddened ' + ytitle
         endif

         djs_plot, [0], [0], /nodata, charsize=csize, $
          xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
          xtitle=xtitle, ytitle=ytitle, title=starname
         djs_contourpts, locuscolor[ic1,*], locuscolor[ic2,*], $
          bin1=0.05, bin2=0.05, xrange=xrange, yrange=yrange, /overplot
;         djs_oplot, locuscolor[ic1,*], locuscolor[ic2,*], $
;          psym=3
         djs_oplot, spcolor[ic1,indx], spcolor[ic2,indx], $
          psym=3, color='red'
         djs_oplot, [targetcolor[ic1]], [targetcolor[ic2]], $
          psym=4, color='blue', symsize=2
         djs_oplot, [bd17color[ic1]], [bd17color[ic2]], $
          psym=2, color='green', symsize=2

         xpos = total(!x.crange*[0.95,0.05])
         ypos = total(!y.crange*[0.05,0.95])
         dy = 0.06 * (!y.crange[1] - !y.crange[0])
         djs_xyouts, xpos, ypos, 'Stellar locus', charsize=csize, color='default'
         djs_xyouts, xpos, ypos-dy, starname, charsize=csize, color='red'
         djs_xyouts, xpos, ypos-2*dy, 'Target priority', charsize=csize, color='blue'
         djs_xyouts, xpos, ypos-3*dy, 'BD+17', charsize=csize, color='green'
      endfor
   endfor

   dfpsclose

   return
end
;------------------------------------------------------------------------------
