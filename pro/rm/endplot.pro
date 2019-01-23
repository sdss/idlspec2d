;+
; NAME:
;       ENDPLOT 
;
; PURPOSE:
;       Closes and saves the current plotfile.  By default, the file
;       is then spooled to the printer using the unix 'lpr -r'
;       command. Keywords are available to name the plotfile and
;       inhibit spooling.
;
;       See also BEGPLOT and KILLPLOT.
;
; CATEGORY:
;
;       Misc.
;
; CALLING SEQUENCE:
;
;       ENDPLOT
;
;       There are no user defined inputs.
;
; INPUTS:
;
;       None.
;
; KEYWORD PARAMETERS:
;
;       /noprint: Set to prevent spooling. This keyword is now ignored
;       /doprint: send to printer after closing
;       P: Name (string) of printer to spool plotfile to.  'lpr -r
;          -Pprintername' is used.
;
;       /shell: Set to use an intervening shell process for the
;                'spawn, "lpr..."' call.  By default, the 
;                /NOSHELL keyword to spawn is used.
;       /setup: run setupplot after closing device
;       /landfix: fix the landscape so its not flipped (which is the only way
;                 IDL writes landscape mode)
;
; OUTPUTS:
;
;       None.
;
; COMMON BLOCKS:
;
;       PLOTFILE_COM:  See the documentation for BEGPLOT.
;
; EXAMPLE:
;
;       See the documentation for BEGPLOT.
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
;       Based on ENDPLT by Phil Klingner, NIS-3, LANL.
;
;       Modified and renamed ENDPLOT. Michael J. Carter, NIS-1, LANL,
;       March, 1994.
;
;       Added /noshell keyword in call to 'spawn' to speed up printing.
;       Use the /shell keyword to 'endplot' to override.  < mjc - April, 1995 >
;       Made not printing the default. /noprint now quietly ignored. Added
;            /doprint keyword to force printing. 
;         Erin Scott Sheldon UofMich 17-Mar-2001
;      Incorporated setupplot stuff
;         E.S.S. 3-Jul-2002
;      Added /landfix: 15-Aug-2003  Erin Sheldon
;
;-

PRO Endplot, doprint=doprint, noprint=noprint, p=p, shell=shell, help=help, setup=setup, landfix=landfix

;; Procedure to close out and print current or named plot file
   
   COMMON Plotfile_com, plotfile, save_device, spool
   
   IF keyword_set(help) THEN BEGIN
      doc_library, 'endplot'
      return
   ENDIF
   
   sv = size(plotfile)
   IF sv(sv(0)+1) EQ 0 THEN plotfile = 'null'
   sv = size(save_device)
   IF sv(sv(0)+1) EQ 0 THEN save_device = !D.name
   
   IF keyword_set(p) THEN $
    parg = '-P' + strtrim(string(p), 2) $
   ELSE $
     parg = ''
   
   message, 'Closing current plot file ' + plotfile, /continue
   IF !D.name NE 'X' THEN device, /close
   set_plot, save_device

   IF ((keyword_set(doprint)) AND (spool EQ 1)) THEN BEGIN
       IF keyword_set(shell) THEN BEGIN
           spawncmd = 'lpr -h -r' + ' ' + parg + ' ' + plotfile
           spawn, spawncmd
       ENDIF ELSE BEGIN
           IF parg EQ '' THEN $
             spawncmd = ['lpr', '-h', '-r', plotfile] $
           ELSE $
             spawncmd = ['lpr', '-h', '-r', parg, plotfile]
           spawn, spawncmd, /noshell
       ENDELSE 
       message, plotfile+' sent to printer', /continue
   ENDIF
   
   pslayout
   IF keyword_set(setup) OR !pslayout.runsetup THEN setupplot

   IF keyword_set(landfix) THEN pslandfix,plotfile

   plotfile = 'null'

END
