
pro write_spflat, flatinfoname, iflat, flatstruct, flathdr, $
                  arcname, nbright, ymodel, scatter, $
                  nowrite=nowrite, writeflatmodel=writeflatmodel, outdir=outdir, $
                  master_flat = master_flat, pad_blue=pad_blue, pad_red=pad_red, $
                  buildTraceFlat=buildTraceFlat

                  

    sxaddpar, flathdr, 'NBRIGHT', nbright, $
          'Number of bright pixels (>10^5) in extracted flat-field'

    sxaddpar, flathdr, 'FLATNAME', flatstruct[iflat].name, $
            'Name of flat'
    sxaddpar, flathdr, 'ARCNAME', arcname[flatstruct[iflat].iarc], $
            'Name of associated arc'
    sxaddpar, flathdr, 'PROFTYPE', flatstruct[iflat].PROFTYPE, $
            'Extract Profile (1:guass 2:ExpCubic 3:DoubleGauss 4:2+3)'
    if keyword_set(master_flat) then begin
        sxaddpar, flathdr, 'MASTERFLAT', FILE_BASENAME(master_flat),'Flat used to pad FFLAT'
    endif else sxaddpar, flathdr, 'MASTERFLAT', '','Flat used to pad FFLAT'
    sxaddpar, flathdr, 'FFBLUEMASTER', pad_blue, 'Raw FFLAT Blue padded'
    sxaddpar, flathdr, 'FFBLUEMASTER', pad_red, 'Raw FFLAT Red padded'
    sxaddpar, flathdr, 'SFLATMIN', flatstruct[iflat].superflat_minval, ' Superflat Minimum'
    quad = ['LL','LR','UL','UR']
    foreach mw,flatstruct[iflat].medwidth, idx do $
            sxaddpar, flathdr, 'MEDWIDT'+strtrim(idx,2), mw, $
                    ' Median spatial dispersion widths in '+quad[idx]+' quadrant'



    if not keyword_set(nowrite) then begin
        if keyword_set(outdir) then begin
            file_mkdir, outdir
            flatinfofile = filepath(string(format='(a,i8.8,a)',flatinfoname, $
                                               sxpar(flathdr, 'EXPOSURE'),'.fits'),$
                                root_dir=outdir)
        endif else begin
            flatinfofile = string(format='(a,i8.8,a)',flatinfoname, $
                sxpar(flathdr, 'EXPOSURE'), '.fits')
        endelse
        mwrfits_named, *flatstruct[iflat].fflat, flatinfofile, hdr = flathdr, name='FFLAT', /create
        mwrfits_named, *flatstruct[iflat].tset, flatinfofile, name='TSET'
        mwrfits_named, *flatstruct[iflat].fibermask, flatinfofile, name='FIBERMASK'
        mwrfits_named, *flatstruct[iflat].widthset, flatinfofile, name='WIDTHSET'
        mwrfits_named, *flatstruct[iflat].superflatset, flatinfofile, name='SUPERFLATSET'
        mwrfits_named, *flatstruct[iflat].xsol, flatinfofile, name='XSOL'
        if keyword_set(buildTraceFlat) then begin
            mwrfits_named, *flatstruct[iflat].flux, flatinfofile, name='FLUX'
            mwrfits_named, *flatstruct[iflat].ivar, flatinfofile, name='IVAR'
            mwrfits_named, *flatstruct[iflat].WSET, flatinfofile, name='WSET'
        endif

        spawn, ['gzip', '-f', flatinfofile], /noshell
        
        ; ASB: write flat image model info if requested:
        if keyword_set(writeflatmodel) then begin
            flatmodelfile = string(format='(a,i8.8,a)',flatinfoname + $
                'MODELIMG-', sxpar(flathdr, 'EXPOSURE'), '.fits')
            mwrfits_named, flatimg, flatmodelfile, name = 'FLATIMG'/create
            mwrfits_named, flativar, flatmodelfile, name = 'FLATIVAR'
            mwrfits_named, ymodel + scatter,flatmodelfile, name = 'ymodel + scatter'
            spawn, ['gzip', '-f', flatmodelfile], /noshell
        endif
    endif
    
    flatstruct[iflat].hdr = ptr_new(flathdr)

    return
end
