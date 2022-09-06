;------------------------------------------------------------------------------

pro update_fitsTable, fitsfile, name, extdescription, hstruct
    fits_info, fitsfile, EXTNAME=fitsexts, textout=fitsfile+'.prt'
    file_delete, fitsfile+'.tmp', fitsfile+'.prt',  /ALLOW_NONEXISTENT, /QUIET

    foreach ext, fitsexts, i do begin
        if strmatch(strtrim(ext,2), name, /fold_case) eq 0 then begin
            fibermap=MRDFITS(fitsfile, i, fits_hdr, /silent)
            MWRFITS, fibermap, fitsfile+'.tmp', fits_hdr, Status=Status, /silent
        endif else begin
            undefine, sumhdr
            sxaddpar, sumhdr, 'EXTNAME', name, extdescription
            MWRFITS, hstruct, fitsfile+'.tmp', sumhdr, Status=Status, /silent
        endelse
    endforeach
    FILE_MOVE, fitsfile+'.tmp', fitsfile, /overwrite
end
