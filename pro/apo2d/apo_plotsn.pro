pro apo_plotsn, plate, plugmapfile, outdir=outdir, plotfile=plotfile

   if (NOT keyword_set(plate)) then return

   platestr = string(plate, format='(i4.4)')

   bands = [1,3]
   snarray = fltarr(2, 640)
   camnames = ['b1','b2','r1','r2']

   ;----------
   ; Read the plug map file for all 640 fibers

   yanny_read, plugmapfile, pdata
   plugmap = *pdata[0]
   yanny_free, pdata
   plugsort = sortplugmap(plugmap, fibermask=fibermask)

   ;----------
   ; Loop through reductions for all science frames, and add S/N
   ; in quadrature, e.g. sum (S/N)^2

   for icam=0, n_elements(camnames)-1 do begin
      outsci = filepath('sci-'+platestr+'-'+camnames[icam]+'-*.fits', $
       root_dir=outdir)
      filenames = findfile(outsci, count=ct)
      meansn = 0
      for ifile=0, ct-1 do begin
         meansn = sqrt( meansn^2 + (mrdfits(filenames[ifile], 2))^2 )
      endfor
      case icam of
         0: snarray[0,0:319] = meansn
         1: snarray[0,320:639] = meansn
         2: snarray[1,0:319] = meansn
         3: snarray[1,320:639] = meansn
      endcase
   endfor

   ;----------
   ; Make the plot

   ; Lock the file to do this.
   if (keyword_set(plotfile)) then $
    while(djs_lockfile(plotfile, lun=plot_lun) EQ 0) do wait, 5

   title = 'APO SPECTRO PLATE=' + strtrim(string(plate),2)
   plotsn, snarray, plugsort, bands=bands, title=title, plotfile=plotfile

   if (keyword_set(plotfile)) then $
    djs_unlockfile, plotfile, lun=plot_lun

   return
end
