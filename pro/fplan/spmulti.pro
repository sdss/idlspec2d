;+
; NAME:
;   spmulti
;
; PURPOSE:
;   Combine all files in 2dmerge directory
;
; CALLING SEQUENCE:
;   spmulti, mjd=mjd
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   mjd        - MJD's to combine
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   cpbackup
;   idlspec2d_version()
;   idlutils_version()
;   spcoadd_frames
;   splog
;   spreduce
;   yanny_free
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   27-Apr-2000 Written by S Burles, FNAL
;-
;------------------------------------------------------------------------------

pro spmulti, mjd=mjd, extractDir=extractDir, mergeDir=mergeDir

   if (NOT keyword_set(extractDir)) then extractDir = ''
   if (NOT keyword_set(mergeDir)) then mergeDir = ''

   stime0 = systime(1)

   ;----------
   ; Open log files for output

   ;if (keyword_set(logfile)) then begin
   ;   cpbackup, logfile
   ;   splog, filename=logfile
   ;   splog, 'Log file ', logfile, ' opened ', systime()
   ;endif

   ;----------
   ; Combine all red+blue exposures for a given sequence

   mjdstr = 'all'
   outputname = string('spMerge2d-',mjdstr,'-',platestr,'.fits', $
                format='(a,a,a1,a4,a5)')

   for side=1, 2 do begin


       look = string('spSpec2d-?',side,'*fits',format='(a,i1,a)')

       files = findfile(filepath(look, root_dir=extractDir))
       j = where(files NE '', nfile)
       if (nfile GT 0) then files = files[j] else files = 0

       splog, 'Combining ' + strtrim(string(nfile),2) $
         + ' files for side ' + strtrim(string(side),2)

       if (nfile GT 0) then begin
          individualroot  = string('spInd2d-',mjdstr,'-',platestr, $
                format='(a,a,a1,a4)')

          for i=0, nfile-1 do splog, 'Combine file ', files[i]

;               outputroot = string('spMerge2d-',mjd,'-',platestr,'.fits', $
;                format='(a,i5.5,a1,a4,a5)')
;               combine2dout, files, filepath(outputroot, root_dir=combineDir), $
;                side, wavemin=alog10(3750.0d), window=100, $
;                individual=filepath(individualroot, root_dir=combineDir)

               spcoadd_frames, files, $
                filepath(outputname, root_dir=combineDir), $
                wavemin=alog10(3750.0d), window=100
            endif

         endfor

         splog, 'Time to combine sequence ', seqid[iseq], ' = ', $
          systime(1)-stime3, ' seconds', format='(a,i5,a,f6.0,a)'
      endif

   endfor ; End loop for sequence number

   splog, 'Total time for SPALLREDUCE = ', systime(1)-stime0, ' seconds', $
    format='(a,f6.0,a)'
   splog, 'Successful completion of SPALLREDUCE at ', systime()

   ;----------
   ; Close log files

   if (keyword_set(plotfile) AND NOT keyword_set(xdisplay)) then begin
      device, /close
      set_plot, 'x'
   endif

   if (NOT keyword_set(skipsn) AND keyword_set(combineDir)) then begin

     if (NOT keyword_set(snfile) AND keyword_set(plotfile)) then begin
       snfile = plotfile
       strput, snfile, 'spSignal'
     endif

     cpbackup, snfile
     set_plot, 'ps'
     device, filename=snfile, /color, /inches, xs=8.0, ys=10.0, xoff=0.25, $
                 yoff=0.4
     splog, 'Plot file ', snfile, ' opened ', systime()

     cd, combineDir, current=old_dir

     expres=outputname
     dot = strpos(outputname, '.')
     strput, expres, '*', dot

     checksn, expres=expres, title=outputname
     device, /close
     set_plot,'x'

     cd, old_dir
   endif  

   if (keyword_set(logfile)) then splog, /close

   return
end
;------------------------------------------------------------------------------
