
;+
; NAME:
;       BEGPLOT 
;
; PURPOSE:
;       This procedure redirects plots from the present display/
;       graphics device to another display/graphics device. It
;       is a 'macro' to the SET_PLOT and DEVICE graphics routines.  At
;       this time only the 'PS' device is fully supported.
;
;       See also KILLPLOT and ENDPLOT.
;
; CATEGORY:
;
;       Misc.
;
; CALLING SEQUENCE:
;
;       BEGPLOT
;
; INPUTS:
;
;       None.
;
; KEYWORD PARAMETERS:
;
;       These keywords are applicable to all graphics/display devices.
;
;       LAYOUT:  Name of file containing keywords.  See example below.
;       DEVICE:  Name of graphics device.
;       NAME:  Name to give output file.
;       INVBW: Invert b & w color table and make axes labels black for
;              printing b & w spectrograms.
;
;       These keywords are applicable to the 'PS', 'PCL', and 'HP'
;       graphics device.  See also "IDL Graphics Devices" in IDL
;       Reference Manual.
;
;       INCHES:  Set to express sizes and offsets in inches.
;       LANDSCAPE:  Set for 'landscape' orientation.  
;       XOFFSET: Offset relative to lower left hand corner in
;                portrait(default) orientation. Default is 0.75 inches.
;       YOFFSET: Offset relative to lower left hand corner in
;                portrait(default) orientation. Default is 1.0 inches.
;       XSIZE:  Width of graphics output. Default = 5.0 inches.
;       YSIZE:  Height of graphics output. Default = 9.0 inches.
;
;       These keywords are specific to the 'PS' graphics device.
;
;       BITS_PER_PIXEL: Integer (1,2,4,8) specifing resolution.
;                       Default = 4 unless COLOR keyword set then
;                       default = 8.
;       BOLD: Set to use bold version of current PostScript font.
;       COLOR: Set to enable color output.
;       DEMI: Set to use demi version of current PostScript font.
;       ENCAPSULATED: Set to create encapsulated PostScript font.
;                     Suitable for importing into other documents.
;       CMYK:         Set to create colors in cmyk format - for ApJ purposes
;       FONT: Name of PostScript font. See IDL Reference Manual for
;             complete list.
;       ITALIC: Set to use italic version of current PostScript font.
;       LIGHT: Set to use light version of current PostScript font.
;       MEDIUM: Set to use meduim version of current PostScript font.
;       NARROW: Set to use narrow version of current PostScript font.
;       OBLIQUE: Set to use oblique versin of current PostScript font.
;       PREVIEW: Set to add 'device independent screen preview' to
;                output file when using ENCAPSULATED.
;       SCALE_FACTOR: Specifies scale factor applied to entire plot.
;                     Default = 1.
;       CT: Number of color table to apply when using COLOR. Default =
;           6 (prism).
;
; COMMON BLOCKS:
;
;       PLOTFILE_COM:
;               plotfile - name of output file.
;               save_device - name of previous graphics device.
;               spool - flag to prevent printing if ENCAPSULATED.
;
; RESTRICTIONS:
;       
;       !P.FONT must be set to 0 (zero) for the hardware-drawn (PostScript)
;       fonts. If !P.FONT = -1, Hershey vector drawn fonts are used.  Vector 
;       fonts should always be used for 3-dimensional plots.
;
; EXAMPLES:
;
;       To redirect a plot to a PostScript printer named 'printer':
;
;       BEGPLOT
;       PLOT,x,y
;       ENDPLOT,P='printer'
;
;       To redirect a plot to an Encapsulated PostScript (EPS) file
;       and inhibit printing:
;
;       BEGPLOT,NAME='test.eps',/ENCAPSUATED
;       PLOT,x,y
;       ENDPLOT,/NOPRINT
;
;       To use a 'layout file':
;
;       BEGPLOT,LAYOUT='layoutfile'
;       PLOT,x,y
;       ENDPLOT
;
;       Layout files can to used to set up the graphics device
;       keywords. A syntax similar to Xdefaults files is used.  Each
;       line looks like "device.keyword:value". 'device' is the
;       name of graphcis device used by SET_PLOT. ('ps' is currently
;       the only device supported).  'keyword' is the name of the
;       keyword to send to DEVICE, and 'value' is the value to send
;       with the keyword (1 is used to set and 0 to unset). A layout
;       file to produce a half-page, portrait, EPS file  using the
;       Times-Roman 12 point font and PREVIEW would look like:
;
;       ps.font:times
;       ps.font_size:12
;       ps.portrait:1
;       ps.encapsulated:1
;       ps.preview:1
;       ps.inches:1
;       ps.xsize:5.0
;       ps.ysize:5.0
;
;       All of the keywords defined above can be set in this manner.
;
;       To abort after running BEGPLOT but before running ENDPLOT, run
;       KILLPLOT. This will close the plotfile and reset the graphics
;       device to the one being used before BEGPLOT was called.
;
;       See also the documentation for KILLPLOT and ENDPLOT.
;
; COPYRIGHT NOTICE:
;
;       Copyright 1993, The Regents of the University of California. This
;       software was produced under U.S. Government contract (W-7405-ENG-36)
;       by Los Alamos National Laboratory, which is operated by the
;       University of California for the U.S. Department of Energy.
;       The U.S. Government is licensed to use, reproduce, and distribute
;       this software. Neither the Government nor the University makes
;       any warranty, express or implied, or assumes any liability or
;       responsibility for the use of this software.
;
; MODIFICATION HISTORY:
;
;       Based on BEGPLT by Phil Klingner, NIS-3, LANL.
;
;       Modified and renamed BEGPLOT. Michael J. Carter, NIS-1, LANL,
;       March, 1994.
;
;       Fixed landscape mode bug.  mjc, 1/24/95.
;
;       Added INVBW keyword to load inverse black  &  white color map - 
;       useful for printing spectrograms.
;  
;       Big rewrite. Use !pslayout system veriable.  Added parameters
;       to layout. Run setupplot if necessary.  
;              3-Jul-2002  Erin S. Sheldon UofMich

;       Added CMYK color option: Desika Narayanan, Steward
;       Observatory, March 9th, 2006
;-

PRO Begplot, layout=layout, $ ; Layout file containing device and keywords.
             device=device, $
             name=name, $  ; Common printer/plotter keywords.
             inches=inches, $
             landscape=landscape, $
             xoffset=xoffset, $
             xsize=xsize, $
             yoffset=yoffset, $
             ysize=ysize, $
             bits_per_pixel=bits_per_pixel, $ ; PostScript specific keywords.
             bold=bold, $
             color=color, $
             demi=demi, $
             encapsulated=encapsulated, $
             cmyk=cmyk,$
             font=font, $
             italic=italic, $
             light=light, $
             medium=medium, $
             narrow=narrow, $
             oblique=oblique, $
             preview=preview, $
             scale_factor=scale_factor, $
             ct=ct, $      ; Color table for color plots.
             invbw=invbw, $
             setup=setup,$
             help=help

   COMMON plotfile_com, plotfile, save_device, spool
   
   spool = 1

   IF keyword_set(help) THEN BEGIN
      doc_library, 'begplot'
      return
   ENDIF

   ;; default plotting 
   if n_elements(CMYK) eq 0 then CMYK=1
   if n_elements(color) eq 0 then color=1

   ;;;;;;;;;;;;;;;;;;;;;;;
   ;; run pslayout
   ;;;;;;;;;;;;;;;;;;;;;;;

   pslayout

   ;; Save current graphics device.
   sv = size(save_device)
   IF sv(sv(0)+1) EQ 0 THEN save_device = !D.name
   
   ;; Clean up leftover plotfile.   
   sv = size(plotfile)
   IF sv(sv(0)+1) EQ 0 THEN plotfile = 'null' ; Null 'plotfile' if no leftovers
   
   IF plotfile NE 'null' THEN BEGIN ; Close device if plotfile defined.
      message, 'Closing current plot file '+plotfile, /continue 
      IF !D.name NE 'X' THEN device, /close
   ENDIF

   ;; Set the new graphics device (always set to 'ps' for now).
   IF keyword_set(device) THEN BEGIN
       device=strlowcase(device)
       CASE device OF
           'ps': graphics_device = !pslayout ;ps is structure for 'PS' device keywords
           
           ELSE: print, 'Invalid graphics device: '+device ; return
       ENDCASE
   ENDIF ELSE BEGIN
       graphics_device = !pslayout
   ENDELSE
   set_plot, graphics_device.name


   ;; Get keywords from 'layoutfile' if keyword 'layout' set.
   defsysv, '!pslayout',EXIST=EXIST
   IF exist THEN BEGIN 
       graphics_device = !pslayout
       name_set=0
   ENDIF ELSE IF keyword_set(layout) THEN BEGIN
      name_set = 0
      layoutfile = layout
      setlayout, layoutfile, graphics_device, name_set
      print,'Using layout file '+layoutfile
   ENDIF ELSE BEGIN
      name_set = 0              ; name_set=1 if 'filename' specified
                                ; in layoutfile.
   ENDELSE

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;; now run setupplot if requested
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   IF !pslayout.runsetup OR keyword_set(setup) THEN setupplot

   ;; Check for keywords which are common to PCL, PostScript, and HP-GL,
   ;; devices and update graphics_device structure.
   IF ((graphics_device.name EQ 'pcl') OR (graphics_device.name EQ 'ps') OR $
       (graphics_device.name EQ 'hp')) THEN BEGIN
      IF keyword_set(name) THEN BEGIN ; Filename given as 'begplot' keyword.
         plotfile = name 
      ENDIF ELSE IF (name_set EQ 1) THEN BEGIN ; Filename set in 'setlayout'.
         plotfile = graphics_device.filename
      END ELSE BEGIN ; Use 'tmpfile' to create filename (default).
         prefix =  '/tmp/'
         suffix = '.' + strlowcase(graphics_device.name)
         plotfile = tmpfile(prefix = prefix, suffix = suffix)
      ENDELSE
      IF keyword_set(inches) THEN  graphics_device.inches = 1
      IF keyword_set(landscape) THEN BEGIN
          IF NOT keyword_set(layout) THEN BEGIN
              graphics_device.yoffset = 10.25
              graphics_device.xsize = 9.5
              graphics_device.ysize = 7
          ENDIF  
          graphics_device.landscape = 1
          graphics_device.portrait = 0
      ENDIF
      IF keyword_set(xoffset) THEN graphics_device.xoffset = xoffset
      IF keyword_set(xsize) THEN graphics_device.xsize = xsize
      IF keyword_set(yoffset) THEN graphics_device.yoffset = yoffset
      IF keyword_set(ysize) THEN graphics_device.ysize = ysize
   ENDIF
   
   ;; Check for keywords specific to the 'ps' graphics device and update
   ;; graphics_device structure.
   IF (graphics_device.name EQ 'ps') THEN BEGIN
      IF keyword_set(bits_per_pixel) THEN $
       graphics_device.bits_per_pixel = bits_per_pixel
      IF keyword_set(color) THEN graphics_device.color = 1
      IF keyword_set(encapsulated) THEN graphics_device.encapsulated = 1
      ;DTN EDITED TO ADD CMYK 03/09/2006
      IF keyword_set(cmyk) THEN graphics_device.cmyk=1
      IF keyword_set(font) THEN graphics_device.font = font
      IF keyword_set(font_index) THEN graphics_device.font_index = font_index
      IF keyword_set(font_size) THEN graphics_device.font_size = $
       font_size
      ;; cannot use preview if not in encap mode
      graphics_device.preview = 0
      IF graphics_device.encapsulated EQ 1 THEN BEGIN 
          IF keyword_set(preview) THEN graphics_device.preview = 1 
      ENDIF
      IF keyword_set(scale_factor) THEN graphics_device.scale_factor = $
       scale_factor
      IF keyword_set(user_font) THEN graphics_device.user_font = user_font
   ENDIF

; Set keywords for PS, PCL, and HP devices based on graphics_device structure.
   device, filename = plotfile
   IF (graphics_device.landscape EQ 1) THEN device, /landscape $
    ELSE device, /portrait
   device, xoffset = graphics_device.xoffset, inches = graphics_device.inches
   device, xsize = graphics_device.xsize, inches = graphics_device.inches
   device, yoffset = graphics_device.yoffset, inches = graphics_device.inches
   device, ysize = graphics_device.ysize, inches = graphics_device.inches
   
; Set keywords specific to PS device based on graphics_device structure.
   device, bits_per_pixel = graphics_device.bits_per_pixel
   device, color = graphics_device.color, bits_per_pixel = 8
   device, encapsulated = graphics_device.encapsulated
   device, cmyk = graphics_device.cmyk
   device, preview = graphics_device.preview
   device, scale_factor = graphics_device.scale_factor
   CASE graphics_device.font OF 
       'avantgarde': device, /avantgarde, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique, $
         isolatin1 = graphics_device.isolatin1
       
       'bkman': device, /bkman, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique, $
         isolatin1 = graphics_device.isolatin1
       
       'courier': device, /courier, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique, $
         isolatin1 = graphics_device.isolatin1
       
       'helvetica': device, /helvetica, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique, $
         isolatin1 = graphics_device.isolatin1
       
       'palatino': device, /palatino, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique, $
         isolatin1 = graphics_device.isolatin1
       
       'schoolbook': device, /schoolbook, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique, $
         isolatin1 = graphics_device.isolatin1
       
       'symbol': device, /symbol, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique
       
       'times': device, /times, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique, $
         isolatin1 = graphics_device.isolatin1
       
       'zapfchancery': device, /zpfchancery, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique, $
         isolatin1 = graphics_device.isolatin1
       
       'zapfdingbats': device, /zapfdingbats, bold = $
         graphics_device.bold, book = graphics_device.book, demi = $
         graphics_device.demi, italic = graphics_device.italic, light = $
         graphics_device.light, medium = graphics_device.medium, narrow = $
         graphics_device.narrow, oblique = graphics_device.oblique
   ENDCASE
   IF graphics_device.encapsulated EQ 1 THEN spool = 0

   if keyword_set(invbw) then begin
       device, color = 1
       loadct,41,file="$IDL_DIR/resource/colors/colors_nis.tbl"
       !P.COLOR=255
       
   endif

   message, 'Plotting directed to file '+plotfile, /continue
      
end
