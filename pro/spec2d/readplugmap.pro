;+
; NAME:
;   readplugmap
;
; PURPOSE:
;   readplugmap file and use append tags as needed
;
; CALLING SEQUENCE:
;   plugmap = readplugmap ( 'plPlugMapM.par' )
;
; INPUTS:
;   plugmap - anonymous plugmap structure which includes header cards
;                and other goodies
;
; OPTIONAL KEYWORDS:
;   deredden - apply dereddening to fiber magnitudes 
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   struct_addtags()
;   yanny_free
;   yanny_read
;   yanny_par()
;
; REVISION HISTORY:
;   29-Jan-2001  Written by S. Burles, FNAL
;-
;------------------------------------------------------------------------------
function readplugmap, filename, deredden=deredden

    yanny_read, filename, pdata, hdr=hdr

    if size(pdata,/tname) EQ 'INT' then return, 0

    rawplug = *pdata[0]
    yanny_free, pdata

    nplug = n_elements(rawplug)

    tagstemplate = { $
        cartid   : 0L, $
        plateid  : 0L, $ 
        tileid   : 0L, $
        raplate  : 0.0, $
        decplate : 0.0, $
        redden_med : fltarr(5), $
        fibersn    : fltarr(3), $
        synthmag   : fltarr(3) }

    tagstemplate.cartid   = (yanny_par(hdr, 'cartridgeId'))[0]
    tagstemplate.plateid  = (yanny_par(hdr, 'plateId'))[0]
    tagstemplate.tileid   = (yanny_par(hdr, 'tileId'))[0]
    tagstemplate.redden_med = (yanny_par(hdr, 'reddeningMed'))
    tagstemplate.raplate  = (yanny_par(hdr, 'raCen'))[0]
    tagstemplate.decplate = (yanny_par(hdr, 'decCen'))[0]

    plugmap = struct_addtags(rawplug, replicate(tagstemplate, nplug))

    if keyword_set(deredden) then begin
       splog, 'Applying reddening vector ', plugmap.redden_med 
       plugmap.mag = plugmap.mag - plugmap.redden_med 
    endif
     
    return, plugmap
end
