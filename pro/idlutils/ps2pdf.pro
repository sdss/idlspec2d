;+
; NAME:
;   ps2pdf
;
; PURPOSE:
;   Convert PS plots to PDF plots
;
; CALLING SEQUENCE:
;   ps2pdf, plotfile
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
;   19-Jan-2024 Witten by S. Morrison, UIUC
;-
;------------------------------------------------------------------------------

pro ps2pdf, plotfile
      
    cmd = 'ps2pdf '+plotfile+ ' '+repstr(plotfile,'.ps','.pdf')
    splog, 'SPAWN '+cmd
    spawn, cmd, sh_out, sh_err
    splog, 'SPAWN out=', sh_out
    splog, 'SPAWN err=', sh_err
    splog, 'Done convert '+plotfile+' to pdf'
      
    return
end
