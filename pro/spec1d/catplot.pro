;+
; NAME:
;   catplot
;
; PURPOSE:
;   Modified version of SPLOT by Tremonti for inspecting spectra.
;
; CALLING SEQUENCE:
;   catplot, [x], y, $
;    [color=, psym=, symsize=, thick= ]
;
;   soplot, [x], y, [/autoscale], $
;    [color=, psym=, symsize=, thick= ]
;
;   sxyouts, x, y, string, [alignment=, charsize=, charthick=, color=, $
;    font=, orientation= ]
;
;   serase, [nerase, /norefresh]
;
; INPUTS:
;
; OUTPUT:
;
; COMMENTS:
;   This code is based upon Aaron Barth's ATV procedure.
;
;   SpInspect added to menu bar.  "Open SpInspect File" brings up a
;   dialogue box.  The spInspect File must have the format 
;   spInspect-pppp-mmmmm-inspector.par.  Upon loading this file, plotspec
;   is called to loop through the fibers on the specified plate.  No
;   changes are written to the spInspect File until "Update SpInspect File"
;   is selected from the menu bar.
;
; EXAMPLES:
;
; BUGS:
;   Doesn't use the passed XRANGE, YRANGE properly yet...
;   Move around widgets to be more compact above plotting window.
;   Write splot_readfits.
;   Make POSITION= changeable based upon CHARSIZE.
;   Gaussian fitting or integrated gaussian fitting.
;   Allow one to step through an image row at a time? Or link to ATV?
;   Use the WCS in splot_gettrack.
;   Add widget button option to fix Y range or let it float, or fix YMIN=0.
;   Include options for plotting contours, etc?
;   Options for XLOG, YLOG
;   For FITS files, take XTITLE, YTITLE from header
;   Option to pass header as param in SPLOT
;   SpInspect only reads the 1st of semi-colon separated variables in
;     the manual_comments field.  (Writes them just fine!) 
;
; PROCEDURES CALLED:
;   fits_read
;
; INTERNAL SUPPORT ROUTINES:
;   splot_gausspix
;   splot_startup
;   splot_clearkeylist
;   splot_displayall
;   splot_readfits
;   splot_writeeps
;   splot_cleartext
;   splot_zoom
;   splot_gettrack
;   splot_event
;   splot_shutdown
;   splot_resize
;   splot_icolor()
;   splot_setheader
;   splot_headinfo
;   splot_headinfo_event
;   splot_plot1plot
;   splot_plot1text
;   splot_plotwindow
;   splot_plotall
;   sxyouts
;   splot_move_cursor
;   splot_set_minmax
;   splot_get_minmax
;   splot_refresh
;   splot_help
;   splot_help_event
;   splot_plotparam_refresh
;   splot_plotparam
;   splot_plotparam_event
;   serase
;   splot_autoscale
;   soplot
;   splot
;   splot_spinspect_new
;   splot_spinspect_update
;   splot_spinspect_event
;   splot_spinspect
;
; REVISION HISTORY:
;   28-Sep-1999  Written by David Schlegel, Princeton.
;   19-Jun-2001  Gaussfit amplitude was wrong - fixed - D. Finkbeiner
;   17-Jun-2002  Added spInspect procedures - C. Tremonti
;                Modified splot_startup, splot_events
;-
;------------------------------------------------------------------------------
; Routine to evaluate a gaussian function integrated over a pixel of width 1,
; plus any number of polynomial terms for a background level
;    A[0] = center of the Gaussian
;    A[1] = sigma width of the Ith Gaussian
;    A[2] = normalization of the Gaussian
;    A[3...] = polynomial coefficients for background terms
; 19 June 2001 factor of sqrt(2*!pi) added to make amplitude correct - DPF

pro splot_gausspix, x, a, f, pder

   ncoeff = N_elements(a)
   fac = size(a, /tname) EQ 'DOUBLE' ? sqrt(2. * !dpi) : sqrt(2. * !pi)
   f = (fac * a[2] * a[1]) * $
     (gaussint((x+0.5-a[0])/a[1]) -gaussint((x-0.5-a[0])/a[1]))

   if (ncoeff GT 3) then begin
      f = f + poly(x, a[3:ncoeff-1])
   endif

   return
end

;------------------------------------------------------------------------------

pro splot_startup 

   common splot_state, state, graphkeys 
   common splot_images, main_image
   common splot_pdata, nplot, maxplot, plot_ptr
   common splot_wcs, astr, equinox

   keylist = { key:' ', x:0.0, y:0.0 }

   state = {                         $
    version: '0.9'                 , $ ; Version number of this release
    head_ptr: ptr_new()            , $ ; Pointer to FITS image header
    imagename: ''                  , $ ; FITS image file name
    base_id: 0L                    , $ ; ID of top-level base
    base_min_size: [600L, 400L]    , $ ; Min size for top-level base
    draw_base_id: 0L               , $ ; ID of base holding draw window
    draw_window_id: 0L             , $ ; Window ID of draw window
    draw_widget_id: 0L             , $ ; Widget ID of draw widget
    location_bar_id: 0L            , $ ; ID of (x,y) label
    wcs_bar_id: 0L                 , $ ; ID of WCS label
    xmin_text_id: 0L               , $ ; ID of XMIN= widget
    xmax_text_id: 0L               , $ ; ID of XMAX= widget
    ymin_text_id: 0L               , $ ; ID of YMIN= widget
    ymax_text_id: 0L               , $ ; ID of YMAX= widget
    comments_text_id: 0L           , $ ; ID of comments output widget
    keyboard_text_id: 0L           , $ ; ID of keyboard input widget
    nkey: 0                        , $ ; Number of elements in keylist
    keylist: replicate(keylist,5)  , $ ; Record of keystrokes + cursor in plot
    xrange: [0.0,1.0]              , $ ; X range of plot window
    yrange: [0.0,1.0]              , $ ; Y range of plot window
    xfix: 'float'                  , $ ; 0=fix, 1=float
    yfix: 'fix'                    , $ ; 0=fix, 1=float
    position: [0.15,0.15,0.95,0.95], $ ; POSITION for PLOT procedure
    draw_window_size: [650L, 512L] , $ ; Size of main draw window
    menu_ids: lonarr(25)           , $ ; List of top menu items
    mouse: [0L, 0L]                , $ ; Cursor position in device coords
    mphys: [0.0, 0.0]              , $ ; Cursor position in data coordinates
    base_pad: [0L, 0L]             , $ ; Padding around draw base
    pad: [0L, 0L]                  , $ ; Padding around draw widget
    headinfo_base_id: 0L           , $ ; headinfo base widget id
    plotparam_base_id: 0L          , $ ; plotting parameter base widget id
    spinspect_base_id: 0L            $ ; plotting parameter base widget id
   }

; Trim this list for the time being... ???
; How to deal with elements that are arrays?
;   graphlist = ['BACKGROUND', 'CHARSIZE', 'XCHARSIZE', 'YCHARSIZE', $
;    'CHARTHICK', 'CLIP', 'COLOR', 'FONT', 'XGRIDSTYLE', 'YGRIDSTYLE', $
;    'LINESTYLE', 'XMARGIN', 'YMARGIN', 'XMINOR', 'YMINOR', 'NOCLIP', $
;    'NORMAL', 'ORIENTATION', 'POSITION', 'PSYM', 'XSTYLE', 'YSTYLE', $
;    'SUBTITLE', 'SYMSIZE', 'THICK', 'XTHICK', 'YTHICK', $
;    'XTICKFORMAT', 'YTICKFORMAT', 'TICKLEN', 'XTICKLEN', 'YTICKLEN', $
;    'XTICKNAME', 'YTICKNAME', 'XTICKS', 'YTICKS', 'XTICKV', 'YTICKV', $
;    'TITLE', 'XTITLE', 'YTITLE']
   graphlist = ['CHARSIZE', 'XCHARSIZE', 'YCHARSIZE', $
    'CHARTHICK', 'XGRIDSTYLE', 'YGRIDSTYLE', $
    'LINESTYLE', 'XMINOR', 'YMINOR', $
    'PSYM', 'XSTYLE', 'YSTYLE', $
    'SYMSIZE', 'THICK', 'XTHICK', 'YTHICK', $
    'TITLE', 'XTITLE', 'YTITLE']
   graphkeys = $
    replicate( { box_id: 0L, $
    keyword: '', value: ptr_new() }, N_elements(graphlist) )
   graphkeys.keyword = graphlist


   nplot = 0
   maxplot = 5000
   plot_ptr = ptrarr(maxplot)

   ; Load a simple color table with the basic 8 colors
   red   = [0, 1, 0, 0, 0, 1, 1, 1]
   green = [0, 0, 1, 0, 1, 0, 1, 1]
   blue  = [0, 0, 0, 1, 1, 1, 0, 1]
   tvlct, 255*red, 255*green, 255*blue

   ; Define the widgets.  For the widgets that need to be modified later
   ; on, save their widget ID's in state variables

   base = widget_base(title = 'splot', $
    /column, /base_align_right, app_mbar = top_menu, $
    uvalue = 'splot_base', /tlb_size_events)
   state.base_id = base

   tmp_struct = {cw_pdmenu_s, flags:0, name:''}
   top_menu_desc = [ $
                  {cw_pdmenu_s, 1, 'File'}, $         ; file menu
                  {cw_pdmenu_s, 0, 'ReadFits'}, $
                  {cw_pdmenu_s, 0, 'WriteEPS'},  $
                  {cw_pdmenu_s, 2, 'Quit'}, $
                  {cw_pdmenu_s, 1, 'Erase'}, $        ; erase menu
                  {cw_pdmenu_s, 0, 'EraseLast'}, $
                  {cw_pdmenu_s, 0, 'EraseAllButFirst'}, $
                  {cw_pdmenu_s, 2, 'EraseAll'}, $
                  {cw_pdmenu_s, 1, 'ImageInfo'}, $    ; info menu
                  {cw_pdmenu_s, 2, 'ImageHeader'}, $
                  {cw_pdmenu_s, 1, 'Plot'}, $         ; plot menu
                  {cw_pdmenu_s, 2, 'PlotParams'}, $
                  {cw_pdmenu_s, 1, 'SpInspect'}, $    ; spInspect menu
                  {cw_pdmenu_s, 0, 'Open SpInspect File'}, $
                  {cw_pdmenu_s, 2, 'Update SpInspect File'}, $
                  {cw_pdmenu_s, 1, 'Help'}, $         ; help menu
                  {cw_pdmenu_s, 2, 'SPLOT Help'} $
                ]

   top_menu = cw_pdmenu(top_menu, top_menu_desc, $
                     ids = state.menu_ids, $
                     /mbar, $
                     /help, $
                     /return_id, $
                     uvalue = 'top_menu')


   button_base1 =  widget_base(base, /row, /base_align_right)
   button_base2 =  widget_base(base, /row, /base_align_right)
   button_base3 =  widget_base(base, /row, /base_align_right)


   state.xmin_text_id = cw_field(button_base1, $
           uvalue = 'xmin_text', /float,  $
           title = 'XMIN=', $
           value = string(state.xrange[0]),  $
           /return_events, $
           xsize = 12)

   state.xmax_text_id = cw_field(button_base1, $
           uvalue = 'xmax_text', /float,  $
           title = 'XMAX=', $
           value = string(state.xrange[1]),  $
           /return_events, $
           xsize = 12)

   fixlist = ['Fix', 'Float']
   mode_droplist_id = widget_droplist(button_base1, $
                                   title = '', $
                                   uvalue = 'xfix', $
                                   value = fixlist)

   state.ymin_text_id = cw_field(button_base2, $
           uvalue = 'ymin_text', /float,  $
           title = 'YMIN=', $
           value = string(state.yrange[0]),  $
           /return_events, $
           xsize = 12)

   state.ymax_text_id = cw_field(button_base2, $
           uvalue = 'ymax_text', /float,  $
           title = 'YMAX=', $
           value = string(state.xrange[1]),  $
           /return_events, $
           xsize = 12)

   fixlist = ['Fix', 'Float']
   mode_droplist_id = widget_droplist(button_base2, $
                                   title = '', $
                                   uvalue = 'yfix', $
                                   value = fixlist)

   tmp_string = string('',format='(a30)')
   state.location_bar_id = widget_label(button_base1, $
                value = tmp_string,  $
                uvalue = 'location_bar', frame=1)

   tmp_string = string('',format='(a30)')
   state.wcs_bar_id = widget_label(button_base2, $
                value = tmp_string,  $
                uvalue = 'wcs_bar', frame=1)

   state.comments_text_id = widget_label(button_base3, $
;           value = '', $
           value = tmp_string, $
           uvalue='comments_text', $
           xsize=state.draw_window_size[0]-90, frame=1)

   zoomone_button = widget_button(button_base3, $
              value = 'Zoom1', $
              uvalue = 'zoom_one')

   done_button = widget_button(button_base3, $
              value = 'Done', $
              uvalue = 'done')

   state.keyboard_text_id =  widget_text(button_base3, $
              /all_events, $
              scr_xsize = 1, $
              scr_ysize = 1, $
              units = 0, $
              uvalue = 'keyboard_text', $
              value = '')

   state.draw_base_id = widget_base(base, $
              /column, /base_align_left, $
              /tracking_events, $
              uvalue = 'draw_base', $
              frame = 2)

   state.draw_widget_id = widget_draw(state.draw_base_id, $
              uvalue = 'draw_window', $
              /motion_events,  /button_events, $
              scr_xsize = state.draw_window_size[0], $
              scr_ysize = state.draw_window_size[1])

   ; Create the widgets on screen
   widget_control, base, /realize

   widget_control, state.draw_widget_id, get_value = tmp_value
   state.draw_window_id = tmp_value

   ; Find window padding sizes needed for resizing routines.
   ; Add extra padding for menu bar, since this isn't included in
   ; the geometry returned by widget_info.
   ; Also add extra padding for margin (frame) in draw base.

   basegeom = widget_info(state.base_id, /geometry)
   drawbasegeom = widget_info(state.draw_base_id, /geometry)
   state.pad[0] = basegeom.xsize - state.draw_window_size[0]
   state.pad[1] = basegeom.ysize - state.draw_window_size[1] + 30
   state.base_pad[0] = basegeom.xsize - drawbasegeom.xsize $
    + (2 * basegeom.margin)
   state.base_pad[1] = basegeom.ysize - drawbasegeom.ysize + 30 $
    + (2 * basegeom.margin)

   xmanager, 'splot', state.base_id, /no_block

   return
end

;------------------------------------------------------------------------------

pro splot_clearkeylist

   common splot_state

   state.nkey = 0
   state.keylist.key = ' '
   state.keylist.x = 0.0
   state.keylist.y = 0.0

   return
end

;------------------------------------------------------------------------------

pro splot_displayall
   splot_refresh
   return
end

;------------------------------------------------------------------------------

pro splot_readfits
; ???

   return
end

;------------------------------------------------------------------------------

pro splot_writeeps

   common splot_state

   widget_control, /hourglass

   filename = dialog_pickfile(filter = '*.eps', $ 
               file = 'splot.eps', $
               group =  state.base_id, $
               /write)
   tmp_result = findfile(filename, count = nfiles)

   result = ''
   if (nfiles GT 0 and filename NE '') then begin
      mesg = strarr(2)
      mesg[0] = 'Overwrite existing file:'
      tmp_string = strmid(filename, strpos(filename, '/', /reverse_search) + 1)
      mesg[1] = strcompress(tmp_string + '?', /remove_all)
      result =  dialog_message(mesg, $
                 /default_no, $
                 dialog_parent = state.base_id, $
                 /question)                 
   endif

   if ((nfiles EQ 0 OR result EQ 'Yes') AND filename NE '') then begin

      screen_device = !d.name
      wset, state.draw_window_id

      aspect_ratio = $
       state.draw_window_size[1] / float(state.draw_window_size[0])
    
      set_plot, 'ps'
      device, $
       filename = filename, $
       /color, $
       bits_per_pixel = 8, $
;       /encapsul, $
       encapsul=0, $
       /inches, $
       xsize = 6.0, $
       ysize = 6.0 * aspect_ratio
    
      splot_plotall

      device, /close
      set_plot, screen_device
   endif

   splot_cleartext

   return
end

;------------------------------------------------------------------------------

pro splot_cleartext

   ; Routine to clear the widget for keyboard input when the mouse is in
   ; the text window.  This de-allocates the input focus from the text
   ; input widget.

   common splot_state

   widget_control, state.draw_base_id, /clear_events
   widget_control, state.keyboard_text_id, set_value = ''

   return
end

;------------------------------------------------------------------------------

pro splot_zoom, zchange, recenter = recenter

   common splot_state

   ; Routine to do zoom in/out and recentering of image

   case zchange of
      'in': begin
         state.xrange = state.mphys[0] $
          + [-0.25, 0.25] * (state.xrange[1] - state.xrange[0])
         if (state.yfix EQ 'float') then splot_autoscale_y
      end
      'out': begin
         state.xrange = state.mphys[0] $
          + [-1.0, 1.0] * (state.xrange[1] - state.xrange[0])
         if (state.yfix EQ 'float') then splot_autoscale_y
      end
      'one': begin
         splot_autoscale_x
         splot_autoscale_y
      end
      'none': begin ; no change to zoom level: recenter on current mouse pos'n
         state.xrange = state.mphys[0] $
          + [-0.5, 0.5] * (state.xrange[1] - state.xrange[0])
      end
      else: print, 'Problem in splot_zoom!'
   endcase

   splot_set_minmax
   splot_refresh

   return
end

;------------------------------------------------------------------------------

pro splot_gettrack

   common splot_state
   common splot_wcs

   ; Update location bar with x, y

   xphysize = state.xrange[1] - state.xrange[0]
   xdevsize = state.draw_window_size[0] $
    * (state.position[2] - state.position[0])
   xdev0 = state.draw_window_size[0] * state.position[0]
   state.mphys[0] = $
    (state.mouse[0] - xdev0) * xphysize / xdevsize + state.xrange[0]

   yphysize = state.yrange[1] - state.yrange[0]
   ydevsize = state.draw_window_size[1] $
    * (state.position[3] - state.position[1])
   ydev0 = state.draw_window_size[1] * state.position[1]
   state.mphys[1] = $
    (state.mouse[1] - ydev0) * yphysize / ydevsize + state.yrange[0]

   loc_string = strcompress( string(state.mphys[0], state.mphys[1]) )

   widget_control, state.location_bar_id, set_value=loc_string

   return
end

;------------------------------------------------------------------------------

pro splot_event, event

  ; Main event loop for SPLOT widgets

   common splot_state
   common splot_pdata
   common splot_spinspect_state, spinspect_state, pdata, yanny_structs, $
                                 yanny_hdr

   widget_control, event.id, get_uvalue=uvalue

   case uvalue of
   'splot_base': begin  ; main window resize: preserve display range
      splot_resize, event
      splot_refresh
      splot_cleartext
   end

   'xfix': case event.index of
      0: state.xfix = 'fix'
      1: state.xfix = 'float'
      else: print, 'Unknown selection!'
   endcase

   'yfix': case event.index of
      0: state.yfix = 'fix'
      1: state.yfix = 'float'
      else: print, 'Unknown selection!'
   endcase

   'top_menu': begin       ; selection from menu bar
      widget_control, event.value, get_value = event_name

      case event_name of

         ; File menu options:
         'ReadFits'  : splot_readfits
         'WriteEPS'  : splot_writeeps
         'Quit'      : splot_shutdown

         ; Erase options:
         'EraseLast' : if (nplot GE 1) then serase, 1
         'EraseAllButFirst' : if (nplot GE 1) then serase, nplot-1
         'EraseAll'  : serase

         ; Info options:
         'ImageHeader': splot_headinfo

         ; Plot options:
         'PlotParams': splot_plotparam

         ; SpInspect options:
	  'Open SpInspect File': splot_spinspect
	  'Update SpInspect File': begin 
               widget_control, spinspect_state.spinspect_file_text_id, $
	                       get_value = filename
               print, 'Writing file ' + filename
	       yanny_write, filename, pdata, structs = yanny_structs, $
	                    hdr = yanny_hdr
             end

         ; Help options:
         'SPLOT Help': splot_help

         else: print, 'Unknown event in file menu!'
      endcase

   end   ; end of file menu options

   ; If the mouse enters the main draw base, set the input focus to
   ; the invisible text widget, for keyboard input.
   ; When the mouse leaves the main draw base, de-allocate the input
   ; focus by setting the text widget value.

   'draw_base': begin
      case event.enter of
      0: begin
         widget_control, state.keyboard_text_id, set_value = ''
         end
      1: begin
         widget_control, state.keyboard_text_id, /input_focus
         end
      endcase
   end

   'draw_window': begin  ; mouse movement or button press

      if (event.type EQ 2) then begin   ; motion event
         tmp_event = [event.x, event.y]
         state.mouse = tmp_event
         splot_gettrack
      endif

      if (event.type EQ 0) then begin
         case event.press of
            1: splot_zoom, 'in', /recenter
            2: splot_zoom, 'none', /recenter
            4: splot_zoom, 'out', /recenter
            else: print,  'trouble in splot_event, mouse zoom'
         endcase
      endif

      widget_control, state.keyboard_text_id, /input_focus

   end

   'xmin_text': begin     ; text entry in 'XMIN= ' box
      splot_get_minmax, uvalue, event.value
      splot_displayall
   end

   'xmax_text': begin     ; text entry in 'XMAX= ' box
      splot_get_minmax, uvalue, event.value
      splot_displayall
   end

   'ymin_text': begin     ; text entry in 'YMIN= ' box
      splot_get_minmax, uvalue, event.value
      splot_displayall
   end

   'ymax_text': begin     ; text entry in 'YMAX= ' box
      splot_get_minmax, uvalue, event.value
      splot_displayall
   end

   'keyboard_text': begin  ; keyboard input with mouse in display window
      eventchar = string(event.ch)

      if (state.nkey LT N_elements(state.keylist)) then begin
         state.keylist[state.nkey].key = eventchar
         state.keylist[state.nkey].x = state.mphys[0]
         state.keylist[state.nkey].y = state.mphys[1]
         state.nkey = state.nkey + 1
      endif else begin
         splot_clearkeylist
      endelse

      case eventchar of
         '1': splot_move_cursor, eventchar
         '2': splot_move_cursor, eventchar
         '3': splot_move_cursor, eventchar
         '4': splot_move_cursor, eventchar
         '6': splot_move_cursor, eventchar
         '7': splot_move_cursor, eventchar
         '8': splot_move_cursor, eventchar
         '9': splot_move_cursor, eventchar
         'g': splot_gaussfit
         else:  ;any other key press does nothing
      endcase
      widget_control, state.keyboard_text_id, /clear_events
   end

   'zoom_one': splot_zoom, 'one'

   'done':  splot_shutdown

   else:  print, 'No match for uvalue....'  ; bad news if this happens

   endcase

   return
end

;------------------------------------------------------------------------------

pro splot_shutdown

   ; Routine to kill the splot window(s) and clear variables to conserve
   ; memory when quitting splot.  Since we can't delvar the splot internal
   ; variables, just set them equal to zero so they don't take up a lot
   ; of space.  Also clear the state and the color map vectors.

   common splot_state
   common splot_images
   common splot_pdata

   if (xregistered ('splot')) then begin
      widget_control, state.base_id, /destroy
   endif

   if (nplot GT 0) then begin
      serase, /norefresh
      plot_ptr = 0
   endif

   if (size(state, /tname) EQ 'STRUCT') then begin
      if (size(state.head_ptr, /tname) EQ 'POINTER') then $
       ptr_free, state.head_ptr
      main_image = 0
      state = 0
   endif

   return
end

;------------------------------------------------------------------------------

pro splot_resize, event

   ; Routine to resize the draw window when a top-level resize event occurs.

   common splot_state

   tmp_event = [event.x, event.y]

   window = (state.base_min_size > tmp_event)

   newbase = window - state.base_pad

   newsize = window - state.pad

  widget_control, state.draw_base_id, $
   xsize = newbase[0], ysize = newbase[1]
  widget_control, state.draw_widget_id, $
   xsize = newsize[0], ysize = newsize[1]

   state.draw_window_size = newsize

   return
end

;------------------------------------------------------------------------------

function splot_icolor, color

   if (n_elements(color) EQ 0) then color='default'

   ncolor = N_elements(color)

   ; If COLOR is a string or array of strings, then convert color names
   ; to integer values
   if (size(color,/tname) EQ 'STRING') then begin ; Test if COLOR is a string

      ; Detemine the default color for the current device
      if (!d.name EQ 'X') then defcolor = 7 $ ; white for X-windows
       else defcolor = 0 ; black otherwise

      icolor = 0 * (color EQ 'black') $
             + 1 * (color EQ 'red') $
             + 2 * (color EQ 'green') $
             + 3 * (color EQ 'blue') $
             + 4 * (color EQ 'cyan') $
             + 5 * (color EQ 'magenta') $
             + 6 * (color EQ 'yellow') $
             + 7 * (color EQ 'white') $
             + defcolor * (color EQ 'default')

     if (!d.N_colors EQ 16777216) then begin
      red   = [0, 1, 0, 0, 0, 1, 1, 1]
      green = [0, 0, 1, 0, 1, 0, 1, 1]
      blue  = [0, 0, 0, 1, 1, 1, 0, 1]
        colors = 255L*red + ishft(255L*green,8) + ishft(255L*blue,16)
        icolor = colors[icolor]
      endif

   endif else begin
      icolor = long(color)
   endelse

   return, icolor
end

;------------------------------------------------------------------------------

pro splot_setheader, head

   ; Routine to set the image header using a pointer to a 
   ; heap variable.  If there is no header (i.e. if SPLOT has just been
   ; passed a data array rather than a filename), then make the
   ; header pointer a null pointer.
   ; The reason for doing it this way is that we don't know in advance
   ; how many lines the header will have, so we can't just save the
   ; header itself as an array in the state structure.  So instead,
   ; save a pointer to the header in the state structure.

   common splot_state
   common splot_wcs

   ; Kill the header info window when a new image is read in

   if (xregistered('splot_headinfo')) then begin
      widget_control, state.headinfo_base_id, /destroy
   endif

   if (n_elements(head) GT 1) then begin
      ptr_free, state.head_ptr
      state.head_ptr = ptr_new(head)

      ; Get astrometry information from header, if it exists
      extast, head, astr, noparams
      if (noparams EQ -1) then astr = 0
    
      equ = get_equinox(head, code)
      if (code EQ -1) then begin
         astr = 0
         equinox = 'J2000'
      endif else begin
         if (equ EQ 2000.0) then equinox = 'J2000'
         if (equ EQ 1950.0) then equinox = 'B1950'
         if (equ NE 2000.0 and equ NE 1950.0) then $
          equinox = string(equ, format = '(f6.4)')
      endelse

   endif else begin
      ; If there's no image header...
      ptr_free, state.head_ptr
      state.head_ptr = ptr_new()
      astr = 0
    
      widget_control, state.wcs_bar_id, set_value = $
       '---No WCS Info---'
   endelse

   return
end

;------------------------------------------------------------------------------

pro splot_headinfo

   common splot_state

   ; If there's no header, kill the headinfo window and exit this routine.
   if (NOT ptr_valid(state.head_ptr)) then begin
      if (xregistered('splot_headinfo')) then begin
         widget_control, state.headinfo_base_id, /destroy
      endif

      mesg = 'No header information available for this image!'
      junk = dialog_message(mesg, /error, $
                          dialog_parent = state.base_id)
      return
   endif

   ; If there is header information but not headinfo window,
   ; create the headinfo window.
   if (NOT xregistered('splot_headinfo')) then begin

      headinfo_base = $
       widget_base(/floating, $
                  /base_align_right, $
                  group_leader = state.base_id, $
                  /column, $
                  title = 'SPLOT image header information', $
                  uvalue = 'headinfo_base')
      state.headinfo_base_id = headinfo_base

      h = *(state.head_ptr)

      headinfo_text = widget_text(headinfo_base, $
                            /scroll, $
                            value = h, $
                            xsize = 85, $
                            ysize = 24)
    
      headinfo_done = widget_button(headinfo_base, $
                              value = 'Done', $
                              uvalue = 'headinfo_done')

      widget_control, headinfo_base, /realize
      xmanager, 'splot_headinfo', headinfo_base, /no_block

   endif

   return
end

;------------------------------------------------------------------------------

pro splot_headinfo_event, event

   common splot_state

   widget_control, event.id, get_uvalue=uvalue

   case uvalue of
      'headinfo_done': widget_control, event.top, /destroy
      else:
   endcase

   return
end

;------------------------------------------------------------------------------

pro splot_plot1plot, iplot

   common splot_state
   common splot_pdata

   widget_control, /hourglass

   ; Convert color names to index numbers
   options = (*(plot_ptr[iplot])).options
   c = where(tag_names(options) EQ 'COLOR', ct)
   if (ct EQ 1) then options.color = splot_icolor(options.color)

;   plot, [(*(plot_ptr[iplot])).x], [(*(plot_ptr[iplot])).y], $
;    /noerase, xstyle=5, ystyle=5, xrange=!x.crange, yrange=!y.crange, $
;    _EXTRA=options
;   oplot, [(*(plot_ptr[iplot])).x], [(*(plot_ptr[iplot])).y], $
;    _EXTRA=options
   plot, [(*(plot_ptr[iplot])).x], [(*(plot_ptr[iplot])).y], $
    /noerase, xstyle=5, ystyle=5, xrange=state.xrange, yrange=state.yrange, $
    _EXTRA=options

   return
end

;----------------------------------------------------------------------

pro splot_plot1text, iplot

   common splot_pdata

   widget_control, /hourglass

   ; Convert color names to index numbers
   options = (*(plot_ptr[iplot])).options
   c = where(tag_names(options) EQ 'COLOR', ct)
   if (ct EQ 1) then options.color = splot_icolor(options.color)

   string_tex = (*(plot_ptr[iplot])).text
;   if (keyword_set(string_tex)) then string_tex = TeXtoIDL(string_tex)
   xyouts, (*(plot_ptr[iplot])).x, (*(plot_ptr[iplot])).y, $
    string_tex, _EXTRA=options

   return
end

;----------------------------------------------------------------------

pro splot_plotwindow

   common splot_state
   common splot_pdata

   ; Set plot window - draw box

   !p.position = state.position ; ???
   if (nplot GT 0) then begin

      ; Always plot the box with color='default'
      iplot = 0
      options = (*(plot_ptr[iplot])).options
      c = where(tag_names(options) EQ 'COLOR', ct)
      if (ct EQ 1) then options.color = splot_icolor('default')

      plot, [0], [0], /nodata, $
       xrange=state.xrange, yrange=state.yrange, xstyle=1, ystyle=1, $
       _EXTRA=options
   endif else begin
      plot, [0], [0], /nodata, $
       xrange=state.xrange, yrange=state.yrange, xstyle=1, ystyle=1
   endelse

   return
end

;----------------------------------------------------------------------

pro splot_plotall
   common splot_state
   common splot_pdata

   ; Routine to overplot line plots from SPLOT and text from SXYOUTS
   splot_plotwindow

   for iplot=0, nplot-1 do begin
      case (*(plot_ptr[iplot])).type of
        'points'  : splot_plot1plot, iplot
        'text'    : splot_plot1text, iplot
        else      : print, 'Problem in splot_plotall!'
      endcase
   endfor

   return
end

;------------------------------------------------------------------------------

pro sxyouts, x, y, text, _EXTRA=options

   common splot_pdata
   common splot_state

   ; Routine to overplot text

   if (NOT xregistered('splot')) then begin
      print, 'You need to start SPLOT first!'
      return
   endif

   if (N_params() LT 3) then begin
      print, 'Too few parameters for SXYOUTS'
      return
   endif

   if (nplot LT maxplot) then begin

      ; Set default font to 1
      if (N_elements(options) EQ 0) then begin
         options = {font: 1}
      endif else begin
         c = where(tag_names(options) EQ 'FONT', ct)
         if (ct EQ 0) then options = create_struct(options, 'font', 1)
      endelse

      pstruct = { $
       type: 'text',   $       ; type of plot
       x: x,             $     ; x coordinate
       y: y,             $     ; y coordinate
       text: text,       $     ; text to plot
       options: options  $     ; plot keyword options
      }

      plot_ptr[nplot] = ptr_new(pstruct)
      nplot = nplot + 1

      wset, state.draw_window_id
      splot_plot1text, nplot-1
   endif else begin
      print, 'Too many calls to SXYOUTS'
   endelse

   return
end

;------------------------------------------------------------------------------

pro splot_move_cursor, direction

   ; Use keypad arrow keys to step cursor one pixel at a time.
   ; Get the new track image, and update the cursor position.

   common splot_state

   i = 1L

   case direction of
      '2': state.mouse[1] = (state.mouse[1] - i) > 0
      '4': state.mouse[0] = (state.mouse[0] - i) > 0
      '8': state.mouse[1] = (state.mouse[1] + i) < (state.draw_window_size[1]-1)
      '6': state.mouse[0] = (state.mouse[0] + i) < (state.draw_window_size[0]-1)
      '7': begin
         state.mouse[1] = (state.mouse[1] + i) < (state.draw_window_size[1]-1)
         state.mouse[0] = (state.mouse[0] - i) > 0
      end
      '9': begin
         state.mouse[1] = (state.mouse[1] + i) < (state.draw_window_size[1]-1)
         state.mouse[0] = (state.mouse[0] + i) < (state.draw_window_size[0]-1)
      end
      '3': begin
         state.mouse[1] = (state.mouse[1] - i) > 0
         state.mouse[0] = (state.mouse[0] + i) < (state.draw_window_size[0]-1)
      end
      '1': begin
         state.mouse[1] = (state.mouse[1] - i) > 0
         state.mouse[0] = (state.mouse[0] - i) > 0
      end

   endcase

   wset,  state.draw_window_id
   tvcrs, state.mouse[0], state.mouse[1], /device

   splot_gettrack

   ; Prevent the cursor move from causing a mouse event in the draw window
   widget_control, state.draw_widget_id, /clear_events

   return
end

;----------------------------------------------------------------------

pro splot_set_minmax

   ; Updates the MIN and MAX text boxes with new values.

   common splot_state

   widget_control, state.xmin_text_id, set_value=state.xrange[0]
   widget_control, state.xmax_text_id, set_value=state.xrange[1]
   widget_control, state.ymin_text_id, set_value=state.yrange[0]
   widget_control, state.ymax_text_id, set_value=state.yrange[1]

end

;----------------------------------------------------------------------

pro splot_get_minmax, uvalue, newvalue

   ; Change the min and max state variables when user inputs new numbers
   ; in the text boxes.

   common splot_state

   case uvalue of

      'xmin_text': begin
         reads, newvalue, tmp1
         state.xrange[0] = tmp1
      end

      'xmax_text': begin
         reads, newvalue, tmp1
         state.xrange[1] = tmp1
      end

      'ymin_text': begin
         reads, newvalue, tmp1
         state.yrange[0] = tmp1
      end

      'ymax_text': begin
         reads, newvalue, tmp1
         state.yrange[1] = tmp1
      end

   endcase

   splot_set_minmax

   return
end

;------------------------------------------------------------------------------

pro splot_refresh

   common splot_state

   widget_control, /hourglass

   ; Display all plots
   wset, state.draw_window_id
   splot_plotall
   splot_plotparam_refresh ; ???

   splot_gettrack

   ; prevent unwanted mouse clicks
   widget_control, state.draw_base_id, /clear_events


   return
end
;------------------------------------------------------------------------------

pro splot_help

   common splot_state

   h =  'SPLOT HELP'
   h = [h, '']
   h = [h, 'MENU BAR:']
   h = [h, 'File->ReadFits:         Read in a new fits image from disk']
   h = [h, 'File->WritePS:          Write a PostScript file of the current display']
   h = [h, 'File->Quit:             Quits SPLOT']
   h = [h, 'Erase->EraseLast:       Erases the most recent plot label']
   h = [h, 'Erase->EraseAll:        Erases all plot labels']
   h = [h, 'ImageInfo->ImageHeader: Display the FITS header, if there is one.']
   h = [h, 'Plot->PlotParams:       Edit plot options']
   h = [h, '']
   h = [h, 'CONTROL PANEL ITEMS:']
   h = [h,'XMIN:            Shows X minimum for display; click to modify']
   h = [h,'XMAX:            Shows X maximum for display; click to modify']
   h = [h,'YMIN:            Shows Y minimum for display; click to modify']
   h = [h,'YMAX:            Shows Y maximum for display; click to modify']
   h = [h, '']
   h = [h,'MOUSE:']
   h = [h,'                 Button1 = Zoom in & center']
   h = [h,'                 Button2 = Center on current position']
   h = [h,'                 Button3 = Zoom out & center']
   h = [h,'BUTTONS:']
   h = [h,'Zoom1:           Rescale XRANGE and YRANGE to show all data']
   h = [h,'Done:            Quits SPLOT']
   h = [h, '']
   h = [h,'Keyboard commands in display window:']
   h = [h,'    Numeric keypad (with NUM LOCK on) moves cursor']
   h = [h,'    g: gaussian fit (gaussian integrated over each pixel)']
   h = [h, '']
   h = [h,'IDL COMMAND LINE HELP:']
   h = [h,'To plot a 1-D array:']
   h = [h,'    splot, xvector, yvector [, options]']
   h = [h,'To plot a 1-D array from a FITS file:']
   h = [h,'    splot, filename [, options] (enclose file name in single quotes)']
   h = [h,'To overplot text on the draw window: ']
   h = [h,'    sxyouts, x, y, text_string [, options]  (enclose string in single quotes)']
   h = [h, '']
   h = [h,'The options for SPLOT and SXYOUTS are essentially']
   h = [h, 'the same as those for the IDL PLOT and XYOUTS commands.']
   h = [h,'The default color is red for overplots done from the IDL command line.']
   h = [h, '']
   h = [h,'Other commands:']
   h = [h,'serase [, N]:       erases all (or last N) plots and text']
   h = [h,'splot_shutdown:   quits SPLOT']
   h = [h,'NOTE: If SPLOT should crash, type splot_shutdown at the IDL prompt.']
   h = [h, '']
   h = [h, '']
   h = [h,strcompress('SPLOT.PRO version '+state.version+' by D. Schlegel')]
   h = [h,'For full instructions, or to download the most recent version, go to:']
   h = [h,'http://www.astro.princeton.edu/~schlegel/index.html']
   h = [h, '']
   h = [h,'For the companion program ATV, go to:']
   h = [h,'http://cfa-www.harvard.edu/~abarth/atv/atv.html']

   if (NOT xregistered('splot_help')) then begin

      helptitle = strcompress('splot v' + state.version + ' help')

      help_base =  widget_base(/floating, $
       group_leader = state.base_id, $
       /column, $
       /base_align_right, $
       title = helptitle, $
       uvalue = 'help_base')

      help_text = widget_text(help_base, $
                            /scroll, $
                            value = h, $
                            xsize = 85, $
                            ysize = 24)
    
      help_done = widget_button(help_base, $
                              value = 'Done', $
                              uvalue = 'help_done')

      widget_control, help_base, /realize
      xmanager, 'splot_help', help_base, /no_block

   endif

   return 
end

;------------------------------------------------------------------------------

pro splot_help_event, event

   widget_control, event.id, get_uvalue = uvalue

   case uvalue of
      'help_done': widget_control, event.top, /destroy
      else:
   endcase

   return
end

;------------------------------------------------------------------------------

pro splot_plotparam_refresh

   common splot_state
   common splot_pdata

   if (nplot GT 0 AND xregistered('splot_plotparam')) then begin
      options = (*(plot_ptr[0])).options ; Test if exists first ???
      opnames = tag_names(options)

      for i=0, N_elements(graphkeys)-1 do begin
         c = where(opnames EQ graphkeys[i].keyword, ct)
         if (ct EQ 1) then $
          widget_control, graphkeys[i].box_id, set_value=options.(c[0])
      endfor
   endif

   return
end

;------------------------------------------------------------------------------

pro splot_plotparam_event, event

   common splot_state
   common splot_pdata

   options = (*plot_ptr[0]).options ; Test if exists first ???

   widget_control, event.id, get_uvalue=uvalue

   case uvalue of

      'plotparam_done': widget_control, event.top, /destroy

      else: begin
         uvalue = strupcase(uvalue)
         if (N_elements(options) EQ 0) then begin
            options = create_struct(uvalue, event.value)
         endif else begin
            c = where(tag_names(options) EQ uvalue, ct)
            if (ct EQ 0) then $
             options = create_struct(options, uvalue, event.value) $
             else options.(c[0]) = event.value
         endelse

         (*plot_ptr[0]) = { $
          type: (*plot_ptr[0]).type, $
          x:    (*plot_ptr[0]).x,    $
          y:    (*plot_ptr[0]).y,    $
          options: options           $
         }

         splot_refresh
      end

   endcase

   return
end

;------------------------------------------------------------------------------

pro splot_plotparam

   common splot_state

   if (NOT xregistered('splot_plotparam')) then begin

      plotparam_base = $
      widget_base(/floating, $
                  /base_align_left, $
                  group_leader = state.base_id, $
                  /column, $
                  title = 'SPLOT plot params', $
                  uvalue = 'plotparam_base')

      tmp_string = string('',format='(a30)')
      for i=0, N_elements(graphkeys)-1 do begin
         graphkeys[i].box_id = $ ; ???
          cw_field(plotparam_base, $
                      /string, $
                      /return_events, $
                      title = graphkeys[i].keyword+':', $
                      uvalue = graphkeys[i].keyword, $
                      value = tmp_string)
      endfor

      plotparam_done = $
       widget_button(plotparam_base, $
                    value = 'Done', $
                    uvalue = 'plotparam_done')

      widget_control, plotparam_base, /realize
      xmanager, 'splot_plotparam', plotparam_base, /no_block
   endif

   splot_plotparam_refresh

   return
end

;------------------------------------------------------------------------------

pro serase, nerase, norefresh=norefresh

   common splot_pdata

   ; Routine to erase line plots from SPLOT and text from SXYOUTS.

   ; The norefresh keyword is used  when a new image
   ; has just been read in, and by splot_shutdown.

   if (N_params() LT 1) then nerase = nplot $
   else if (nerase GT nplot) then nerase = nplot

   for iplot=nplot-nerase, nplot-1 do begin
      ptr_free, plot_ptr[iplot]
      plot_ptr[iplot] = ptr_new()
   endfor

   nplot = nplot - nerase

   if (NOT keyword_set(norefresh) ) then splot_refresh

   return
end

;------------------------------------------------------------------------------

pro splot_gaussfit

   common splot_state
   common splot_pdata

   i = where(state.keylist.key EQ 'g', ct)
   if (ct EQ 1) then begin
      widget_control, state.comments_text_id, $
       set_value='GAUSSFIT: Press g at other side of feature to fit'
   endif else if (ct EQ 2) then begin
      ; Select all data points in the first PDATA array within the
      ; selected X boundaries.
      xmin = min([state.keylist[i].x, state.keylist[i].x])
      xmax = max([state.keylist[i].x, state.keylist[i].x])
      j = where((*(plot_ptr[0])).x GE xmin $
       AND (*(plot_ptr[0])).x LE xmax)
      if (N_elements(j) GT 3) then begin
         xtemp = (*(plot_ptr[0])).x[j]
         xmid = median(xtemp)
         xtemp = xtemp - xmid ; Do this for numerical stability in the fit
         ytemp = (*(plot_ptr[0])).y[j]
         ymin = min(ytemp)
         ymax = max(ytemp, imax)

         ; Set initial guess for fitting coefficients
         ; Fit a gaussian + a constant sky term
         a = [xtemp[imax], 0.2*(max(xtemp)-min(xtemp)), ymax-ymin, ymin]
         yfit = curvefit(xtemp, ytemp, xtemp*0+1.0, a, $
          /noderivative, function_name='splot_gausspix', chisq=chisq)

         ; If an absorption line is a better fit, use that
         a_abs = [xtemp[imax], 0.2*(max(xtemp)-min(xtemp)), ymin-ymax, ymax]
         yfit_abs = curvefit(xtemp, ytemp, xtemp*0+1.0, a_abs, $
          /noderivative, function_name='splot_gausspix', chisq=chisq_abs)
         if (chisq_abs LT chisq AND chisq_abs NE 0) then begin
            yfit = yfit_abs
            a = a_abs
         endif

         xtemp = xtemp + xmid
         a[0] = a[0] + xmid
         area = a[2] * sqrt(2.*!pi)
         out_string = 'GAUSSFIT: ' $
          + ' x0= ' + strtrim(string(a[0]),2) $
          + ' sig= ' + strtrim(string(a[1]),2) $
          + ' Area= ' + strtrim(string(area),2) $
          + ' sky= ' + strtrim(string(a[3]),2)
         widget_control, state.comments_text_id, set_value=out_string
         soplot, xtemp, yfit, color='red', psym=10
      endif else begin
         widget_control, state.comments_text_id, $
          set_value='GAUSSFIT: Too few points to fit'
      endelse
      splot_clearkeylist
   endif else begin
      splot_clearkeylist
   endelse

   return
end

;------------------------------------------------------------------------------

pro splot_autoscale_x

   common splot_state
   common splot_pdata

   ; First set default values if no data
   state.xrange[0] = 0.0
   state.xrange[1] = 1.0

   if (nplot GT 0) then begin

      ; Set plotting limits for first SPLOT
      state.xrange[0] = min( (*plot_ptr[0]).x )
      state.xrange[1] = max( (*plot_ptr[0]).x )

      ; Enlarge plotting limits if necessary for other calls to SOPLOT
      for i=1, nplot-1 do begin
         if ((*plot_ptr[i]).type EQ 'points') then begin
            state.xrange[0] = min( [state.xrange[0], (*plot_ptr[i]).x] )
            state.xrange[1] = max( [state.xrange[1], (*plot_ptr[i]).x] )
         endif
      endfor
   endif

   splot_set_minmax

   return
end

;------------------------------------------------------------------------------

pro splot_autoscale_y

   common splot_state
   common splot_pdata

   ; When determining Y plotting limits, only look at those data points
   ; within the X plotting limits.
   if (state.xrange[1] GE state.xrange[0]) then xrange = state.xrange $
    else xrange = state.xrange[[1,0]]

   ; First set default values if no data
   state.yrange[0] = 0.0
   state.yrange[1] = 1.0

   if (nplot GT 0) then begin

      ; Set plotting limits for first SPLOT
      igood = where( (*plot_ptr[0]).x GE xrange[0] $
       AND (*plot_ptr[0]).x LE xrange[1])
      if (igood[0] NE -1) then begin
         state.yrange[0] = min( (*plot_ptr[0]).y[igood] )
         state.yrange[1] = max( (*plot_ptr[0]).y[igood] )
      endif

      ; Enlarge plotting limits if necessary for other calls to SOPLOT
      for i=1, nplot-1 do begin
         if ((*plot_ptr[i]).type EQ 'points') then begin
            igood = where( (*plot_ptr[i]).x GE xrange[0] $
             AND (*plot_ptr[i]).x LE xrange[1])
            if (igood[0] NE -1) then begin
               state.yrange[0] = $
                min( [state.yrange[0], (*plot_ptr[i]).y[igood]] )
               state.yrange[1] = $
                max( [state.yrange[1], (*plot_ptr[i]).y[igood]] )
            endif
         endif
      endfor
   endif

   splot_set_minmax

   return
end

;------------------------------------------------------------------------------

pro soplot, x, y, autoscale=autoscale, replot=replot, $
 xrange=xrange, yrange=yrange, $
 position=position, _EXTRA=options

   common splot_state
   common splot_pdata

   if (N_params() LT 1) then begin
      print, 'Too few parameters for SOPLOT'
      return
   endif

   if (NOT xregistered('splot')) then begin
      print, 'Must use SPLOT before SOPLOT'
   endif

   if (keyword_set(xrange)) then state.xrange = xrange
   if (keyword_set(yrange)) then state.yrange = yrange
   if (keyword_set(position)) then state.position = position

   if (nplot LT maxplot) then begin

      if (N_elements(options) EQ 0) then $
       options = create_struct('color', 'default')

      if (N_params() EQ 1) then begin
         pstruct = { $
          type: 'points',   $     ; points
          x: lindgen(N_elements(x)),      $     ; x coordinate
          y: x,             $     ; y coordinate
          options: options  $     ; plot keyword options
         }
      endif else begin
         pstruct = { $
          type: 'points',   $     ; points
          x: x,             $     ; x coordinate
          y: y,             $     ; y coordinate
          options: options  $     ; plot keyword options
         }
      endelse
      plot_ptr[nplot] = ptr_new(pstruct)
      nplot = nplot + 1

      wset, state.draw_window_id

      if (keyword_set(autoscale) AND NOT keyword_set(xrange)) then $
       splot_autoscale_x
      if (keyword_set(autoscale) AND NOT keyword_set(yrange)) then $
         splot_autoscale_y

      if (keyword_set(replot)) then begin
         splot_plotall
      endif else begin
         splot_plot1plot, nplot-1
      endelse
   endif else begin
      print, 'Too many calls to SPLOT'
   endelse

   return
end

;------------------------------------------------------------------------------

pro catplot, x, y, _EXTRA=KeywordsForSOPLOT

   common splot_state
   common splot_images
   common splot_pdata

   if (N_params() LT 1) then begin
      print, 'Too few parameters for SPLOT'
      return
   endif

   imagename = ''
   head = ''

   ; If X is a filename, read in the file
   if ((N_params() NE 0) AND (size(x, /tname) EQ 'STRING') ) then begin
      imagename = x
      fits_read, imagename, main_image, head
      x = main_image[*,0] ; Set X equal to the first row of the image ???
   endif

   if (NOT xregistered('splot')) then $
     splot_startup

   state.imagename = imagename
   splot_setheader, head

   if (N_params() EQ 1) then begin
      xplot = lindgen(N_elements(x))
      yplot = x
   endif else begin
      xplot = x
      yplot = y
   endelse

   serase, /norefresh
   soplot, xplot, yplot, /autoscale, /replot, _EXTRA=KeywordsForSOPLOT
   splot_plotparam_refresh

   return
end
;------------------------------------------------------------------------------

pro splot_spinspect_new

    common splot_spinspect_state, spinspect_state, pdata

    i = spinspect_state.fiberid - 1
    if i gt 639 then return

    widget_control, spinspect_state.class_text_id, $
                    set_value = (*pdata[0])[i].class
    widget_control, spinspect_state.subclass_text_id, $
                    set_value = (*pdata[0])[i].subclass
    widget_control, spinspect_state.z_text_id, $
                    set_value = (*pdata[0])[i].z

    widget_control, spinspect_state.manual_class_text_id, $
                    set_value = (*pdata[0])[i].manual_class
    widget_control, spinspect_state.manual_subclass_text_id, $
                    set_value = (*pdata[0])[i].manual_subclass
    widget_control, spinspect_state.manual_z_text_id, $
                    set_value = (*pdata[0])[i].manual_z
    widget_control, spinspect_state.manual_comments_text_id, $
                    set_value = (*pdata[0])[i].manual_comments

    plotspec, spinspect_state.plateid, spinspect_state.fiberid, $
              nsmooth = spinspect_state.nsmooth

    return
end

;------------------------------------------------------------------------------

pro splot_spinspect_update

    common splot_spinspect_state, spinspect_state, pdata

    i = spinspect_state.fiberid - 1
    if i gt 639 then return

    (*pdata[0])[i].inspector =spinspect_state.inspector

    widget_control, spinspect_state.manual_class_text_id, $
                    get_value = manual_class
    (*pdata[0])[i].manual_class = manual_class

    widget_control, spinspect_state.manual_subclass_text_id, $
                    get_value = manual_subclass
    (*pdata[0])[i].manual_subclass = manual_subclass

    widget_control, spinspect_state.manual_z_text_id, $
                    get_value = manual_z
    (*pdata[0])[i].manual_z = manual_z

    widget_control, spinspect_state.manual_comments_text_id, $
                    get_value = manual_comments
    (*pdata[0])[i].manual_comments = manual_comments
  
    return
end

;------------------------------------------------------------------------------

pro splot_spinspect_event, event

   common splot_spinspect_state, spinspect_state, pdata, yanny_structs, $
                                 yanny_hdr

   widget_control, event.id, get_uvalue=uvalue

   case uvalue of
     'spinspect_file_text': begin
         widget_control, spinspect_state.spinspect_file_text_id, $
	   get_value = filename

         spinspect_state.fiberid = 1
         spinspect_state.plateid = long(strmid(filename, 10, 4))
         spinspect_state.mjd = long(strmid(filename, 15, 5))
         inspector = (strmid(filename, 21, strlen(filename) - 21 - 4))[0]
	 spinspect_state.inspector = inspector
	 help, inspector, spinspect_state.inspector_text_id
         widget_control, spinspect_state.inspector_text_id, $
	   set_value = 'Inspector: ' +  inspector

         print, 'Reading file ' + filename
         yanny_read, filename, pdata, structs = yanny_structs, hdr = yanny_hdr
	 splot_spinspect_new
       end

     'nsmooth_text': begin
         widget_control, spinspect_state.nsmooth_text_id, get_value = nsmooth
	 spinspect_state.nsmooth = nsmooth
         plotspec, spinspect_state.plateid, spinspect_state.fiberid, $
                   nsmooth = spinspect_state.nsmooth
       end

     'inspected': begin
         print, 'Inspected Fiber ' + $
	        string(spinspect_state.fiberid, format = '(I3)')
	 splot_spinspect_update
	 print, (*pdata[0])[spinspect_state.fiberid - 1]
	 if spinspect_state.fiberid ge 640 then return
	 spinspect_state.fiberid = spinspect_state.fiberid + 1
	 splot_spinspect_new
       end

     'skip': begin
         print, 'Fiber ' + $
                string(spinspect_state.fiberid, format = '(I3)') + $
	        ' NOT inspected!'
	 print, (*pdata[0])[spinspect_state.fiberid - 1]
	 if spinspect_state.fiberid ge 640 then return
	 spinspect_state.fiberid = spinspect_state.fiberid + 1
	 splot_spinspect_new
       end

     'back': begin
	 if spinspect_state.fiberid le 1 then return
	 spinspect_state.fiberid = spinspect_state.fiberid - 1
         print, 'Moving back to fiber ' + $
                string(spinspect_state.fiberid, format = '(I3)') 
	 splot_spinspect_new
     end

     'goto': begin
         widget_control, spinspect_state.goto_text_id, $
	   get_value = newfiber
	 if newfiber lt 1 or newfiber gt 640 then return
	 spinspect_state.fiberid = newfiber
         print, 'Moving to Fiber ' + string(newfiber, format = '(I3)')
	 splot_spinspect_new
       end
   endcase

   return
end
;------------------------------------------------------------------------------

pro splot_spinspect

   common splot_state
   common splot_spinspect_state, spinspect_state, pdata

   if (NOT xregistered('splot_spinspect')) then begin

     spinspect_state = { $
       plateid: 0L                        , $ ; Current Plate id #
       fiberid: 0L                        , $ ; Current Fiber id #
       mjd: 0L                            , $ ; MJD of plate
       inspector: ""                      , $ ; Name of Inspector
       nsmooth: 5                         , $ ; Pixels to smooth by
       spinspect_file_text_id: 0L         , $ ; ID of spInspect File = widget
       inspector_text_id: 0L              , $ ; ID of Inspector = widget
       class_text_id: 0L                  , $ ; ID of Class = widget
       subclass_text_id: 0L               , $ ; ID of Subclass = widget
       z_text_id: 0L                      , $ ; ID of Z = widget
       manual_class_text_id: 0L           , $ ; ID of Manual Class = widget
       manual_subclass_text_id: 0L        , $ ; ID of Manual Subclass = widget
       manual_z_text_id: 0L               , $ ; ID of Manual Z = widget
       manual_comments_text_id: 0L        , $ ; ID of Manual Comments = widget
       nsmooth_text_id: 0L                , $ ; ID of Smooth By = widget
       goto_text_id: 0L                     $ ; ID Goto Fiber = widget
     }

     spinspect_base = widget_base(/floating, $
        /base_align_left, $
        group_leader = state.base_id, $
        /column, $
        title = 'Splot spInspect Params', $
        uvalue = 'spinspect_base')

     button_base1 =  widget_base(spinspect_base, /row, /base_align_left)
     button_base2 =  widget_base(spinspect_base, /row, /base_align_left)
     button_base3 =  widget_base(spinspect_base, /row, /base_align_left)
     button_base4 =  widget_base(spinspect_base, /row, /base_align_left)
     button_base5 =  widget_base(spinspect_base, /row, /base_align_left)
     button_base6 =  widget_base(spinspect_base, /row, /base_align_left)

     spinspect_state.spinspect_file_text_id = cw_field(button_base1, $
       uvalue = 'spinspect_file_text', /string,  $
       title = ' SpInspect File:', $
       value = '',  $
       xsize = 50, /return_events)

     spinspect_state.inspector_text_id = widget_label(button_base2, $
       uvalue = 'inspector_text',  $
       value = '  Inspector: ' + spinspect_state.inspector, $
       xsize = 150, /align_left )

     spinspect_state.class_text_id = cw_field(button_base3, $
       uvalue = 'class_text',  $
       title = 'Spec1d Class:',  $
       xsize = 8, /string, /noedit)

     spinspect_state.subclass_text_id = cw_field(button_base3, $
       uvalue = 'subclass_text',  $
       title = ' Subclass:',  $
       xsize = 12, /string, /noedit)

     spinspect_state.z_text_id = cw_field(button_base3, $
       uvalue = 'z_text',  $
       title = ' Z:',  $
       xsize = 8, /float, /noedit)

     spinspect_state.manual_class_text_id = cw_field(button_base4, $
       uvalue = 'manual_class_text', /string,  $
       title = 'Manual Class:', $
       value =  '',  $
       xsize = 8)
 
     spinspect_state.manual_subclass_text_id = cw_field(button_base4, $
       uvalue = 'manual_subclass_text', /string,  $
       title = ' Subclass:', $
       value = '',  $
       xsize = 12)
  
     spinspect_state.manual_z_text_id = cw_field(button_base4, $
       uvalue = 'manual_z_text', /float,  $
       title = ' Z:', $
       value = 0.0,  $
       xsize = 8)
 
     spinspect_state.manual_comments_text_id = cw_field(button_base5, $
       uvalue = 'manual_comments_text', /string,  $
       title = ' Comments:', $
       value = ' ',  $
       xsize = 37)

    spinspect_state.nsmooth_text_id = cw_field(button_base5, $
       uvalue = 'nsmooth_text', /int,  $
       title = ' Smooth by:', $
       value = spinspect_state.nsmooth,  $
       /return_events, $
       xsize = 4)

     inspected_button = widget_button(button_base6, $
       value = ' Fiber Inspected ', ysize = 28, $
       uvalue = 'inspected')

     skip_button = widget_button(button_base6, $
       value = ' Skip Fiber ', ysize = 28, $
       uvalue = 'skip')

     back_button = widget_button(button_base6, $
       value = ' Back 1 Fiber ', ysize = 28, $
       uvalue = 'back')

     spinspect_state.goto_text_id = cw_field(button_base6, $
       uvalue = 'goto', /int,  $
       title = ' Go to Fiber:', $
       value = 1,  $
       /return_events, $
       xsize = 4)

     widget_control, spinspect_base, /realize
     xmanager, 'splot_spinspect', spinspect_base, /no_block
   endif

   return
end
