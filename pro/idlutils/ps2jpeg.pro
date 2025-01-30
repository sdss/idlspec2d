;+
; NAME:
;   ps2jepg
;
; PURPOSE:
;   Convert PS plots to jpeg plots
;
; CALLING SEQUENCE:
;   ps2jpeg, plotfile
;
; INPUTS:
;   plotfile      - PS plot filename
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   spawn
;   splog
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   22-Jan-2025 Witten by S. Morrison, UIUC
;-
;------------------------------------------------------------------------------

pro ps2jpeg, plotfile, jpegfiletmp, jpegfile, flag = flag, extra = extra
    if not file_test(plotfile) then return
    if not keyword_set(extra) then extra = ''
    if not keyword_set(flag) then flag = ''
    cmd = '/usr/bin/convert '+flag+' '+plotfile+' '+jpegfiletmp+' ; \mv '+jpegfiletmp+' '+jpegfile+extra+' &'
    splog, 'SPAWN '+cmd
    spawn, cmd, sh_out, sh_err
    splog, 'SPAWN out=', sh_out
    splog, 'SPAWN err=', sh_err
    splog, 'Done convert '+plotfile+' to pdf'
      
    return
end
