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
                ' Median spectral dispersion width in '+quad[idx]+' quadrant'

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
        
        ; Reformating rejline into a structure for writting to fits file
        ; previous it was convering the strings to numeric ASCII arrays 
        rejline_temp = *arcstruct[iarc].rejline
        n = N_ELEMENTS(rejline_temp)
        REJLINE = REPLICATE({LAMBDA: 0.0, REJLINE:' '}, n)
        FOR i = 0, n-1 DO begin
            val = rejline_temp[i]
            if strlen(val) eq 0 then val = ' '
            REJLINE[i].REJLINE = val
            REJLINE[i].LAMBDA = lambda[i]
        ENDFOR
        
        mwrfits_named, flux, arcinfofile, hdr = archdr, name = 'FLUX', /create
        mwrfits_named, [transpose(lambda), xpeak], arcinfofile, name = 'LAMBDA'
        mwrfits_named, *arcstruct[iarc].wset, arcinfofile, name = 'WSET'
        mwrfits_named, *arcstruct[iarc].fibermask, arcinfofile, name = 'FIBERMASK'
        mwrfits_named, *arcstruct[iarc].dispset, arcinfofile, name = 'DISPSET'
        mwrfits_named, *arcstruct[iarc].reslset, arcinfofile, name = 'RESLSET'
        mwrfits_named, REJLINE, arcinfofile, name = 'REJLINE'
        mwrfits_named, *arcstruct[iarc].XDIF_TSET, arcinfofile, name = 'XDIF_TSET'

        ;width = fltarr(n_elements(lambda), ntrace)
        ;width[ilamp, *] = width_final
        ;mwrfits_named, width, arcinfofile, name='WIDTH' ;--- !!!!!!!!!!!!! for debug purposes only
  
        spawn, ['gzip', '-f', arcinfofile], /noshell

        ; ASB: write arc image model info if requested:
        if keyword_set(writearcmodel) then begin
            arcmodelfile = string(format='(a,i8.8,a)',arcinfoname + $
                        'MODELIMG-', sxpar(archdr, 'EXPOSURE'), '.fits')
            mwrfits_named, arcimg, arcmodelfile, name= 'ARCIMG', /create
            mwrfits_named, arcivar, arcmodelfile, name= 'ARCIVAR'
            mwrfits_named, ymodel, arcmodelfile, name= 'YMODEL'
            spawn, ['gzip', '-f', arcmodelfile], /noshell
        endif
    endif
    ymodel = 0
    return
end
