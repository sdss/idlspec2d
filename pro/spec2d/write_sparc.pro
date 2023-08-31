pro write_sparc, arcinfoname, iarc, arcstruct, archdr, $
     flatname, arcname, fbadpix, bestcorr, tai, $
     lambda, xpeak, ntrace, width_final, ilamp, $
     arcimg, arcivar, ymodel, nowrite=nowrite, $
     writearcmodel=writearcmodel, outdir=outdir

    sxaddpar, archdr, 'FBADPIX', fbadpix, 'Fraction of bad pixels in raw image'
    sxaddpar, archdr, 'BESTCORR', bestcorr, 'Best Correlation coefficient'
    sxaddpar, archdr, 'TAI', tai, 'Tai of arc'
    sxaddpar, archdr, 'TSEP', arcstruct[iarc].tsep, 'Tai seperation from associated flat'
    sxaddpar, archdr, 'FLATNAME', flatname, 'Name of associated flat
    sxaddpar, archdr, 'ARCNAME', arcname[iarc], 'Name of arc
    sxaddpar, archdr, 'NMATCH',arcstruct[iarc].nmatch, 'Number of lamp lines traced'
    
    quad = ['LL','LR','UL','UR']
    foreach mw,arcstruct[iarc].medwidth, idx do $
        sxaddpar, archdr, 'MEDWIDT'+strtrim(idx,2), mw, $
                ' Median spectral dispersion widths in '+quad[idx]+' quadrant'

    foreach mw,arcstruct[iarc].MEDRESOL, idx do $
        sxaddpar, archdr, 'MEDREX'+strtrim(idx,2), mw, $
                ' Median resolution widths in '+quad[idx]+' quadrant'

    arcstruct[iarc].hdr = ptr_new(archdr)

    if not keyword_set(nowrite) then begin
        if keyword_set(outdir) then begin
            file_mkdir, outdir
            arcinfofile = filepath(string(format='(a,i8.8,a)',arcinfoname, $
                                               sxpar(archdr, 'EXPOSURE'),'.fits'),$
                                root_dir=outdir)
        endif else begin
            arcinfofile = string(format='(a,i8.8,a)',arcinfoname, $
                                    sxpar(archdr, 'EXPOSURE'), '.fits')
        endelse
        
        sxaddpar, archdr, 'EXTNAME', 'flux'
        mwrfits, flux, arcinfofile, archdr, /create

        sxaddpar, hdr, 'EXTNAME', 'lambda'
        mwrfits, [transpose(lambda), xpeak], arcinfofile

        sxaddpar, hdr, 'EXTNAME', 'wset'
        mwrfits, *arcstruct[iarc].wset, arcinfofile

        sxaddpar, hdr, 'EXTNAME', 'fibermask'
        mwrfits, *arcstruct[iarc].fibermask, arcinfofile

        sxaddpar, hdr, 'EXTNAME', 'dispset'
        mwrfits, *arcstruct[iarc].dispset, arcinfofile

        sxaddpar, hdr, 'EXTNAME', 'reslset'
        mwrfits, *arcstruct[iarc].reslset, arcinfofile

        sxaddpar, hdr, 'EXTNAME', 'rejline'
        mwrfits, *arcstruct[iarc].rejline, arcinfofile

        sxaddpar, hdr, 'EXTNAME', 'XDIF_TSET'
        mwrfits, *arcstruct[iarc].XDIF_TSET, arcinfofile

        sxaddpar, hdr, 'EXTNAME', 'FIBERMASK'
        mwrfits, *arcstruct[iarc].FIBERMASK, arcinfofile

        ;width = fltarr(n_elements(lambda), ntrace)
        ;width[ilamp, *] = width_final
        ;mwrfits, width, arcinfofile ;--- !!!!!!!!!!!!! for debug purposes only
  
        spawn, ['gzip', '-f', arcinfofile], /noshell

        ; ASB: write arc image model info if requested:
        if keyword_set(writearcmodel) then begin
            arcmodelfile = string(format='(a,i8.8,a)',arcinfoname + $
                        'MODELIMG-', sxpar(archdr, 'EXPOSURE'), '.fits')
            mwrfits, arcimg, arcmodelfile, /create
            mwrfits, arcivar, arcmodelfile
            mwrfits, ymodel, arcmodelfile
            spawn, ['gzip', '-f', arcmodelfile], /noshell
        endif
    endif
    ymodel = 0
    return
end
