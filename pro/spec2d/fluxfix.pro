pro fluxfix

   spec = '1'

   sciname = findfile('spFrame-?'+spec+'*.fits')
   spflux, sciname

   for iname=0, n_elements(sciname)-1 do begin
      spawn, 'cp '+sciname[iname]+' '+sciname[iname]+'.BAK'

      hdr = headfits(sciname[iname])
      expstr = string(sxpar(hdr,'EXPOSURE'), format='(i8.8)')
      platestr = string(sxpar(hdr,'PLATEID'), format='(i4.4)')
      mjdstr = string(sxpar(hdr,'MJD'), format='(i5.5)')
      camname = strtrim(sxpar(hdr,'CAMERAS'),2)

      objflux = mrdfits(sciname[iname],0) * mrdfits(sciname[iname],6)
      objivar = mrdfits(sciname[iname],1) * mrdfits(sciname[iname],6)
      wset = mrdfits(sciname[iname],3)
      traceset2xy, wset, xx, loglam

      calibset = $
       mrdfits('spFluxcalib-'+platestr+'-'+mjdstr+'-'+camname+'.fits',1)
      calibfac = bspline_valu(loglam,calibset)

      corrset = mrdfits('spFluxcorr-'+expstr+'-'+spec+'.fits',1)
      traceset2xy, corrset, loglam, corrimg

      divideflat, objflux, objivar, calibfac, minval=0.05*mean(calibfac)
      divideflat, objflux, objivar, 1.0/corrimg, minval=0.05*mean(1.0/corrimg)

      modfits, sciname[iname], objflux, exten_no=0
      modfits, sciname[iname], objivar, exten_no=1
   endfor

end

