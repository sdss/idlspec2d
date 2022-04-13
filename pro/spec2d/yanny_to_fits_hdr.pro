function yanny_to_fits_hdr, hdr
    @plugmapkeys.idl

    nhead=n_elements(hdr)
    for i=0, nhead-1 do begin
        if strmid(hdr[i],0,1) EQ "#" then continue
        if strmid(hdr[i],0,8) EQ "EVILSCAN" then continue
        key = (str_sep( strtrim(hdr[i],2), ' '))[0]
        if strlen(key) eq 0 then continue
        if not key_match_dict.haskey(key) then begin
            splog, key+' not saved to fibermap fits header'
            continue
        endif
        matched_key =key_match_dict[key]
        val=yanny_par(hdr, key)
        if n_elements(val) gt 1 then begin
            outval=''
            for j =0, n_elements(val)-1 do outval=outval + (yanny_par(hdr, key))[j]+' '
            val=strmid(outval,0,strlen(outval)-1)
        endif
        sxaddpar, fits_hdr, matched_key, val, ' '+key
    endfor
    if strtrim(yanny_par(hdr,'observatory'),2) eq 'APO' then begin
        sxaddpar, fits_hdr, 'CARTID', 'FPS-N', ' cartridgeId'
    endif else sxaddpar, fits_hdr, 'CARTID', 'FPS-S', ' cartridgeId'
    return, fits_hdr
end
