PRO setupplot, type, help=help, test=test, true=true, invbw=invbw, $
               greyscale=greyscale

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    SETUPPLOT
;       
; PURPOSE:
;    Set up default plotting parameters. Parameters are taken from
;    the !pslayout system variable, created by the pslayout procedure.
;    The user should put a copy of pslayout in their path and change
;    the defaults to suit them. Note pslayout contains tags for setting
;    up X and Z buffer as well as PS
;
;    Also runs simpctable, to create a set of colors, and defsymbols
;    to define the system variables !tsym for true-type font symbols
;    and !vsym for vector drawn font symbols, and !csym for use with
;    either true-type or vector drawn fonts.
;
; CALLING SEQUENCE:
;    setupplot [type, /help, /test, true=true, /invbw]
;
; INPUTS: 
;    NONE
;
; KEYWORD PARAMETERS:
;    /true: Can set /true to use true-type fonts. Can be used to override
;           the x_true and ps_true tags in !pslayout
;    /test: run a test showing all the symbols created by defsymbols.pro
;    /help: print simple syntax/help.
;    /invbw: flip colors
;       
; OPTIONAL INPUTS:
;
;    type: if given, the type is set using set_plot, type
;
; OUTPUTS: 
;    None unless /test, in which case some plots are made.
;
; OPTIONAL OUTPUTS:
;    None
;
; CALLED ROUTINES:
;    PSLAYOUT
;    SIMPCTABLE
;    DEFSYMBOLS
; 
; PROCEDURE: 
;    
;       
;
; REVISION HISTORY:
;    19-Mar-2001 Erin Scott Sheldon UofMich
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF keyword_set(help) THEN BEGIN 
     print,'-Syntax: setupplot, display_type, help=help, test=test, true=true, $'
     print,'   invbw=invbw, greyscale=greyscale'
     print
     print,' Will setup plotting parameters; sets system variables'
     print,' for various plotting symbols. '
     print,' use type="ps" for postscript "x" for x-window'
     print,' If type is not given, then it is determined from the !d.flags'
     print,'Use doc_library,"msetupplot"  for more help.'  
     return
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Check input device type (optional)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(type) NE 0 THEN BEGIN 
      IF datatype(type) NE 'STR' THEN BEGIN
          print,'type must be a string'
          return
      ENDIF 
      set_plot,type
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; run pslayout. Will define !pslayout (if not already defined)
  ;; which contains defaults for postscript output and other stuff,
  ;; including whether or not we should use true-type fonts (or 
  ;; postscript fonts if device is 'ps')
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; the .true tag is obsolete, x_true and ps_true are better, 
  ;; can control each 
  pslayout
  IF n_elements(true) EQ 0 THEN BEGIN
      IF !pslayout.true EQ 1 THEN true=1 ELSE true=0
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; set up a simple color table
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(greyscale) THEN BEGIN 
      loadct,0
  ENDIF ELSE BEGIN 
      simpctable
  ENDELSE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; set default background to white if requested
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(invbw) THEN BEGIN 
      !p.background=!white
      !p.color = !black
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; set system variables in device dependent way
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  CASE !d.name OF
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Defaults for the postscript output
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      'PS': BEGIN 
          !p.thick = !pslayout.ps_thick
          !x.thick = !pslayout.ps_xthick
          !y.thick = !pslayout.ps_ythick
          !x.ticklen = !pslayout.ps_xticklen
          !y.ticklen = !pslayout.ps_yticklen
          !p.charsize = !pslayout.ps_charsize
          !p.charthick = !pslayout.ps_charthick

          ;; symbols/font defined same way for true and postscript fonts
          IF keyword_set(true) OR !pslayout.ps_true THEN !p.font=0

      END
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Defaults for the X-windows display
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      'X': BEGIN
          !p.thick = !pslayout.x_thick
          !x.thick = !pslayout.x_xthick
          !y.thick = !pslayout.x_ythick
          !x.ticklen = !pslayout.x_xticklen
          !y.ticklen = !pslayout.x_yticklen
          !p.charsize = !pslayout.x_charsize
          !p.charthick= !pslayout.x_charthick

          ;; use true-type fonts in X?
          IF keyword_set(true) OR !pslayout.x_true THEN BEGIN
              !p.font = 1 

              ;; set the default font
              ;; make a dummy window
;              IF display_exists() THEN Begin 
;                  window,/free,/pixmap,xsize=1,ysize=1
;                  fset=!pslayout.font
;                  IF !pslayout.bold THEN fset=fset+' bold'
;                  IF !pslayout.italic THEN fset=fset+' italic'
;                  device,set_font=fset,/tt_font
;                  wdelete,!d.window
;              Endif

          ENDIF ELSE BEGIN
              !p.font=-1
          ENDELSE 

      END 
      'Z': BEGIN
          !p.thick = !pslayout.z_thick
          !x.thick = !pslayout.z_xthick
          !y.thick = !pslayout.z_ythick
          !x.ticklen = !pslayout.z_xticklen
          !y.ticklen = !pslayout.z_yticklen
          !p.charsize = !pslayout.z_charsize
          !p.charthick= !pslayout.z_charthick

          ;; use true-type fonts in X?
          IF keyword_set(true) OR !pslayout.x_true THEN BEGIN
              !p.font = 1 
              fset=!pslayout.font
              IF !pslayout.bold THEN fset=fset+' bold'
              IF !pslayout.italic THEN fset=fset+' italic'
              device,set_font=fset,/tt_font
          ENDIF ELSE BEGIN
              !p.font=-1
          ENDELSE 
          
          ;; set the resolution of the z-buffer
          device,set_resolution = !pslayout.z_resolution
      END 
      ELSE: message,'type '+type+' unknown'
  Endcase
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Define new system variables to aid plotting 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; golden ratio: send to aplot
  defsysv,'!gratio', exists=exists
  IF NOT exists THEN defsysv,'!gratio',1.36603

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; define the plotting symbols: !tsym for true and !vsym for vector
  ;; and the common symbols in !csym using the set !p.font to determine
  ;; which to use
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  defsymbols

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; demonstrate the symbols if requested 
; go 3 columns per page
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dopage = 1
  column = 1
  IF keyword_set(test) THEN BEGIN 
      
      IF type EQ 'PS' THEN size=0.85 ELSE size=1.5

      IF keyword_set(true) THEN BEGIN
          print,'Testing true type fonts: switching to times'
          sym=!TSYM 
          symtags = tag_names(sym)
          symtags = '!!TSYM.' + symtags
          firstmess='Symbols defined in !!TSYM system variable'
      ENDIF ELSE BEGIN
          print,'Testing vector drawn fonts: simplex roman'
          sym=!VSYM
          symtags = tag_names(sym)
          symtags = '!!VSYM.' + symtags
          firstmess='Symbols defined in !!VSYM system variable'
      ENDELSE 
      
      ntags = n_elements(symtags)

      plot,[0],/nodata,ystyle=4,xstyle=4
      xyouts,0,1,firstmess

      xstep = .35
      ystep = .1

      ystart = 0.9
      y=ystart
      x = 0.0
      FOR i=0, ntags-1 DO BEGIN 

          IF ( (i+1) MOD 10) EQ 0 THEN BEGIN

              column=column+1
              IF (column EQ 4) OR (i EQ 0)  THEN BEGIN
                  column=1
                  x = 0
                  IF strupcase(type) EQ 'X' THEN key=get_kbrd(1)
                  plot,[0],/nodata,ystyle=4,xstyle=4
              ENDIF ELSE  x = x + xstep
              y = ystart
             
              xyouts, x, y, symtags[i]+'  '+SYM.(i),charsize=size
          ENDIF ELSE BEGIN
              IF i NE 0 THEN y = y - ystep
              xyouts, x, y, symtags[i]+'  '+SYM.(i),charsize=size
          ENDELSE 
      ENDFOR 

  ENDIF 

  return 
END 
