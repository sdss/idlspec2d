function strct_to_struct, inStructure,inIndx, inTag,$
                          outStructure,outIndx, outTag=outTag, $
                          nomatch=nomatch, altTag=altTag, altOutTag=altOutTag
    nomatch = 0
    if not keyword_set(outTag) then outTag=inTag

    if (tag_exist(outStructure,outTag) AND tag_exist(inStructure,inTag)) then begin
        inTaginx = where(strmatch(tag_names(inStructure), inTag,/FOLD_CASE) eq 1)
        outTaginx = where(strmatch(tag_names(outStructure), outTag,/FOLD_CASE) eq 1)
        if isa(outIndx, /string) then begin
            if isa(inIndx,/string) then outStructure.(outTaginx) = inStructure.(inTaginx) $
            else outStructure.(outTaginx) = inStructure[inIndx].(inTaginx)
        endif else begin
            if isa(inIndx,/string) then outStructure[outIndx].(outTaginx) = inStructure.(inTaginx) $
            else outStructure[outIndx].(outTaginx) = inStructure[inIndx].(inTaginx)
        endelse
    endif else begin
        if keyword_set(altTag) then begin
            if not keyword_set(altOutTag) then altOutTag=altTag
            if (tag_exist(outStructure,altOutTag) AND tag_exist(inStructure,altTag)) then begin
                inTaginx = where(strmatch(tag_names(inStructure), altTag,/FOLD_CASE) eq 1)
                outTaginx = where(strmatch(tag_names(outStructure), altOutTag,/FOLD_CASE) eq 1)

                if isa(inIndx,/string) then outStructure[outIndx].(outTaginx) = inStructure.(inTaginx) $
                else outStructure[outIndx].(outTaginx) = inStructure[inIndx].(inTaginx)
            endif else nomatch=1
        endif else nomatch=1
    endelse
    return, outStructure
end
