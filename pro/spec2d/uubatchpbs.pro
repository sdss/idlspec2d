;+
; NAME:
;   uubatchpbs
;
; PURPOSE:
;   Batch process Spectro-2D and Spectro-1D reductions based upon
;   already-built plan files.
; ADAPTED FROM: batchpbs.pro
;   To function on clusters without node sharing (e.g. University of Utah)
;
; CALLING SEQUENCE:
;   uubatchpbs, [ platenums, topdir=, run2d=, run1d=, platestart=, plateend=, $
;    mjd=, mjstart=, mjend=, upsvers2d=, upsvers1d=, /zcode, queue=, /skip2d, /clobber, $
;    pbsnodes=pbsnodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, pbs_walltime=pbs_walltime, ember=ember]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   platenums  - Plate numbers to reduce; default to '*'
;   topdir     - Optional override value for the environment
;                variable $BOSS_SPECTRO_REDUX.
;   run2d      - Optional override value for the environment variable $RUN2D
;   run1d      - Optional override value for the environment variable $RUN1D
;   platestart - Starting plate number.
;   plateend   - Ending plate number.
;   mjd        - MJD dates to reduce; default to all.
;                Select based upon the MJD of the combine plan file, and
;                reduce data from all nights needed for that combined plate+MJD.
;   mjstart    - Starting MJD dates to reduce.
;   mjend      - Ending MJD dates to reduce.
;   upsvers2d  - If set, then do a "setup idlspec2d $UPSVERS2D" on the
;                remote machine before executing Spectro-2D.  This allows
;                you to batch jobs using a version other than that which
;                is declared current under UPS.
;   upsvers1d  - If set, then do a "setup idlspec2d $UPSVERS2D" on the
;                remote machine before executing Spectro-1D.  This allows
;                you to batch jobs using a version other than that which
;                is declared current under UPS.
;   zcode      - If set, run Zcode in auto mode.
;   queue      - If set, sets the submit queue.
;   skip2d     - If set, then skip the Spectro-2D reductions.
;   clobber    - If set, then reduce all specified plates.  The default is
;                to not reduce plates where the script file already exists.
;   pbs_nodes  - If set, collect the pbs qsub commands into pbs_nodes script files
;                in order to run on clusters without node sharing (ie Utah).
;                default to node sharing, and keep the pbs qsub commands in the
;                individual plate-mjd SCRIPT files.
;   pbs_ppn    - If set, use #PBS -l nodes=1:ppn=pbs_ppn, otherwise
;                default to #PBS -l nodes=1
;   pbs_a      - If set, use #PBS -A pbs_a, otherwise
;                default to none
;   pbs_walltime - If set, use #PBS -l walltime=pbs_walltime, otherwise
;                default to none
;   ember      - If set, then setup the defaults for the ember cluster at the University of Utah:
;                pbs_nodes = 8 (for 8 nodes, without node sharing) 
;                pbs_ppn = 12 (12 processors per node)
;                pbs_a = 'bolton-em' (Bolton's account, limited to 8 nodes: ember253 - ember260)
;                pbs_walltime='48:00:00'
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
;   This is currently written to batch all Spectro-2D and Spectro-1D
;   reductions, without the option of doing only one or the other.
;
;   Script files are generated for partial reductions of a plate,
;   although those jobs are not submitted.  It would be less confusing
;   to not make those script files.
;
; REVISION HISTORY:
;   17-Jan-2010  Written by D. Schlegel, LBL
;   01-Jan-2011  Adapted from batchpbs by Joel R. Brownstein, University of Utah, 
;                to generalize to cluster computers that do not have pbs node sharing
;                by relocating the PBS commands to bundled script files, 
;                in general via the keywords pbs_nodes, pbs_ppn, pbs_a
;                and with University of Utah defaults preset via the keyword ember.
;-
;------------------------------------------------------------------------------
pro uubatchpbs, platenums1, topdir=topdir1, run2d=run2d1, run1d=run1d1, $
 platestart=platestart, plateend=plateend, $
 mjd=mjd, mjstart=mjstart, mjend=mjend, $
 upsvers2d=upsvers2d, upsvers1d=upsvers1d, zcode=zcode, $
 queue=queue, skip2d=skip2d, clobber=clobber, $ 
 pbs_nodes=pbs_nodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, $
 pbs_walltime=pbs_walltime, ember=ember

   if (size(platenums1,/tname) EQ 'STRING') then platenums = platenums1 $
    else if (keyword_set(platenums1)) then $
      platenums = string(platenums1,format='(i4.4)') $
    else platenums = '*'

   ;----------
   ; Determine the top-level of the output directory tree

   if (keyword_set(topdir1)) then topdir = topdir1 $
    else topdir = getenv('BOSS_SPECTRO_REDUX')
   if strpos(topdir,'/',strlen(topdir)-1) lt 0 then topdir+='/'
   splog, 'Setting TOPDIR=', topdir
   if (keyword_set(run2d1)) then run2d = strtrim(run2d1,2) $
    else run2d = getenv('RUN2D')
   splog, 'Setting RUN2D=', run2d
   if (keyword_set(run1d1)) then run1d = strtrim(run1d1,2) $
    else run1d = getenv('RUN1D')
   splog, 'Setting RUN1D=', run1d

   topdir2d = djs_filepath('', root_dir=topdir, subdir=run2d)
   if (keyword_set(run1d)) then run1dstr = ',run1d="'+run1d+'"' $
    else run1dstr = ''
   if (keyword_set(run2d)) then run2dstr = ',run2d="'+run2d+'"' $
    else run2dstr = ''

   if keyword_set(ember) then begin
     if not keyword_set(pbs_nodes) then pbs_nodes=8
     if not keyword_set(pbs_ppn) then pbs_ppn=12
     if not keyword_set(pbs_walltime) then pbs_walltime='48:00:00'
     if not keyword_set(pbs_a) then pbs_a = 'bolton-em'
   endif
      
   ;----------
   ; Create list of plate directories
   ; Limit the list to only those specified by PLATENUMS,PLATESTART,PLATEEND

   platedirs = get_mjd_dir(topdir2d, mjd=platenums, mjstart=platestart, $
    mjend=plateend)

   if (NOT keyword_set(platedirs[0])) then begin
      splog, 'No directories found'
      return
   endif
   ndir = n_elements(platedirs)

   ;----------
   ; In each plate directory, find all 'spPlancomb*.par' files

   for idir=0L, ndir-1L do begin
      planfile = findfile( $
       djs_filepath('spPlancomb*.par', root_dir=topdir2d, $
        subdir=platedirs[idir]), count=nfile)

      for ifile=0, nfile-1 do begin
         yanny_read, planfile[ifile], hdr=hdr
         thismjd = long(yanny_par(hdr, 'MJD'))

         ; Decide if THISMJD is within the bounds specified by MJD,MJSTART,MJEND
         if (mjd_match(thismjd, mjd=mjd, mjstart=mjstart, mjend=mjend)) $
          then begin
            if (keyword_set(platelist)) then begin
               platelist = [platelist, platedirs[idir]]
               planlist = [planlist, planfile[ifile]]
            endif else begin
               platelist = platedirs[idir]
               planlist = planfile[ifile]
            endelse
         endif
      endfor
   endfor

   nplate = n_elements(planlist)
   if (nplate EQ 0) then begin
      splog, 'No plan files found'
      return
   endif

; Do not use spPlancomb files that only have a subset of the MJDs !!!???

   ;----------
   ; For each combine plan file, generate the IDL script files

   fullscriptfile = strarr(nplate)
   plateid = lonarr(nplate)
   mjd_beg = lonarr(nplate)
   mjd_end = lonarr(nplate)

   fq = "'"
   qbatch = bytarr(nplate) + 1B ; default to reduce all plates

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; PBS_NODES:
   ; Setup the bundled script files if pbs_nodes keyword is set
   if keyword_set(pbs_nodes) then begin

     if nplate lt pbs_nodes then pbs_nodes = nplate

     home = getenv('HOME')
     if (home ne '') then begin
      pos0 = strpos(home,'/',/reverse_search)+1
      pos1 = strlen(home)
      userID = strmid(home,pos0,pos1-pos0)
     endif else userID = 'user'
     print, 'UUBATCHBPS: Starting for user:  ',userID
     
     pbs_root_dir = getenv('BOSS_PBS')
     if (pbs_root_dir eq '') then pbs_root_dir = djs_filepath('pbs/'+run2d,root_dir=topdir) $
     else pbs_root_dir = djs_filepath('bossredux/'+run2d,root_dir=pbs_root_dir)
     pbs_dir = djs_filepath(userID,root_dir=pbs_root_dir) + '/'
     if file_test(pbs_dir) then begin
       shift_pbs_dir = djs_filepath(userID+'.*',root_dir=pbs_root_dir) + '/'
       shift_pbs = file_search(shift_pbs_dir, count=nshift_pbs)
       max_shift = -1L
       for i=0,nshift_pbs-1 do begin
         pos0 = strpos(shift_pbs[i],'/'+userID+'.',/reverse_search)+strlen('/'+userID+'.')
         pos1 = strlen(shift_pbs[i])
         next_shift = fix(strmid(shift_pbs[i],pos0,pos1-pos0))
         max_shift = (next_shift gt max_shift) ? next_shift : max_shift
       endfor
       shift_pbs_dir = djs_filepath(userID + '.' + strtrim(max_shift+1,2),root_dir=pbs_root_dir) + '/'
       print, 'UUBATCHBPS: Renaming previous PBS directory to: '+shift_pbs_dir
       file_move, pbs_dir, shift_pbs_dir
       file_mkdir, pbs_dir
     endif else file_mkdir, pbs_dir
         
     pbs_node_index = 'node'+ strtrim(indgen(pbs_nodes),2)
     pbs_node_script = djs_filepath(pbs_node_index+'.pbs',root_dir=pbs_dir)
     pbs_node_lun = intarr(pbs_nodes)
     if keyword_set(pbs_ppn) then begin
       pbs_ppn_index  = '_proc'+ strtrim(indgen(pbs_ppn),2)
       pbs_ppn_script = strarr(pbs_nodes,pbs_ppn)
     endif
     
     for pbs_node = 0, pbs_nodes-1 do begin
       openw, get_lun, pbs_node_script[pbs_node] ,/get_lun
       pbs_node_lun[pbs_node] = get_lun
       printf, pbs_node_lun[pbs_node], '# Auto-generated batch file '+systime()
       if keyword_set(pbs_a) then printf, pbs_node_lun[pbs_node], '#PBS -A '+pbs_a
       if keyword_set(pbs_walltime) then printf, pbs_node_lun[pbs_node], '#PBS -l walltime='+pbs_walltime
       printf, pbs_node_lun[pbs_node], '#PBS -W umask=0022'
       printf, pbs_node_lun[pbs_node], '#PBS -V'
       printf, pbs_node_lun[pbs_node], '#PBS -j oe'
       if (keyword_set(queue)) then printf, pbs_node_lun[pbs_node], '#PBS -q ' + queue
       if keyword_set(pbs_ppn) then begin
         printf, pbs_node_lun[pbs_node], '#PBS -l nodes=1:ppn='+strtrim(pbs_ppn,2)
         pbs_ppn_script[pbs_node,*] = djs_filepath(pbs_node_index[pbs_node] + pbs_ppn_index +'.pbs',root_dir=pbs_dir)
         for pbs_proc = 0, pbs_ppn-1 do printf, pbs_node_lun[pbs_node], 'source '+pbs_ppn_script[pbs_node,pbs_proc] + ' &'
       endif else printf, pbs_node_lun[pbs_node], '#PBS -l nodes=1"
     endfor 

   endif
   pbs_node = 0
   pbs_proc = 0
   pbs_ppn_append = 0
   
   for iplate=0, nplate-1 do begin
      ; Find all relevant 2D plan files
      yanny_read, planlist[iplate], hdr=hdr
      planfile2d = yanny_par(hdr, 'planfile2d')
      plateid[iplate] = yanny_par(hdr, 'plateid')
      mjd = yanny_par(hdr, 'MJD')
      platemjd = string(plateid[iplate],format='(i4.4)') + '-' $
       + string(mjd,format='(i5.5)')
      platefile = 'spPlate-'+platemjd+'.fits'

      ; Track the beginning and ending MJD going into this plate
      mjd_beg[iplate] = min(strmid(planfile2d,14,5))
      mjd_end[iplate] = max(strmid(planfile2d,14,5))

      ; Split the combine plan file name into a directory and file name
      planfilecomb = fileandpath(planlist[iplate], path=pathcomb)

      ; Construct the name of the batch file
      fullscriptfile[iplate] = djs_filepath('script-'+platemjd, $
       root_dir=pathcomb)
      if (keyword_set(skip2d)) then fullscriptfile[iplate] += '-' + run1d

      ; Write the batch file
      if (keyword_set(clobber) EQ 0) then $
       qbatch[iplate] = file_test(fullscriptfile[iplate]) EQ 0
       
      if (qbatch[iplate]) then begin
         openw, olun, fullscriptfile[iplate], /get_lun
         printf, olun, '# Auto-generated batch file '+systime()
         if not keyword_set(pbs_nodes) then begin
           if keyword_set(pbs_ppn) then printf, olun, '#PBS -l nodes=1:ppn='+strtrim(pbs_ppn,2) $
            else printf, olun, '#PBS -l nodes=1"
           if keyword_set(pbs_a) then printf, olun, '#PBS -A '+pbs_a
           printf, olun, '#PBS -l walltime=48:00:00'
           printf, olun, '#PBS -W umask=0022'
           printf, olun, '#PBS -V'
           printf, olun, '#PBS -j oe'
           ; set queue if asked
           if (keyword_set(queue)) then printf, olun, '#PBS -q ' + queue
           printf, olun, 'cd $PBS_O_WORKDIR'
         endif else printf, olun, 'cd '+pathcomb


         ; Echo commands to make debugging easier
         printf, olun, 'set -o verbose'

         if (keyword_set(skip2d) EQ 0) then begin
            ; Set up requested code version
            if (keyword_set(upsvers2d)) then $
             printf, olun, 'setup idlspec2d '+upsvers2d

            ; Create sorted photoPlate files
            for i=0, n_elements(planfile2d)-1 do $
             printf, olun, 'echo '+fq+'sdss_plate_sort,"'+planfile2d[i]+'"'+fq+' | idl'

            ; Run Spectro-2D
            for i=0, n_elements(planfile2d)-1 do $
             printf, olun, 'echo '+fq+'spreduce2d,"'+planfile2d[i]+'"'+fq+' | idl'
            printf, olun, 'echo '+fq+'spcombine_v5,"'+planfilecomb+'"'+fq+' | idl'
         endif

         ; Run Spectro-1D
         if (keyword_set(upsvers1d)) then $
          printf, olun, 'setup idlspec2d '+upsvers1d
         printf, olun, 'echo '+fq+'spreduce1d,"'+platefile+'"'+run1dstr+fq+' | idl'

         ; Run Zcode
         if (keyword_set(zcode)) then begin
            printf, olun, 'setup runz'
            printf, olun, 'runz_BOSS.sh ' + platefile +' -a'
         endif
         close, olun
         free_lun, olun
         
         ;----------
         ; Do not reduce any plan files that are only partial reductions
         ; (If we tried, then we would get multiple instances of SPREDUCE2D
         ; running on the same data.)
         ; Make this decision regardless of the values of CLOBBER.
         
         indx = where(plateid EQ plateid[iplate] AND mjd_end GT mjd_end[iplate] $
          AND mjd_beg LE mjd_end[iplate], ct)
         if (ct GT 0) then qbatch[iplate] = 0B

         ; Run the scriptfile from a bundled PBS script file (node) 
         ; in the absence of node sharing
         if (keyword_set(pbs_nodes)) then begin
           script_cmd = (ct GT 0) ? '#skip ' : 'source ' 
           if keyword_set(pbs_ppn) then begin
             openw, pbs_ppn_lun, pbs_ppn_script[pbs_node,pbs_proc], append=pbs_ppn_append, /get_lun
             printf, pbs_ppn_lun, script_cmd+fullscriptfile[iplate]+' > '+fullscriptfile[iplate]+'.o'
             close, pbs_ppn_lun
             free_lun, pbs_ppn_lun
             pbs_node += 1
             if pbs_node ge pbs_nodes then begin
               pbs_node = 0
               pbs_proc += 1
               if pbs_proc ge pbs_ppn then begin
                 pbs_proc = 0 
                 if not pbs_ppn_append then pbs_ppn_append = 1
               endif
             endif
           endif else begin
             printf, pbs_node_lun[pbs_node], script_cmd+fullscriptfile[iplate]+' &
             pbs_node += 1
             if pbs_node ge pbs_nodes then pbs_node = 0           
           endelse
         endif
      endif

   endfor
   
   if not pbs_ppn_append then begin
     pbs_node_i = pbs_node
     pbs_proc_i = pbs_proc
     for pbs_proc = pbs_proc_i,pbs_ppn-1 do begin
       for pbs_node = pbs_node_i,pbs_nodes-1 do begin
         openw, pbs_ppn_lun, pbs_ppn_script[pbs_node,pbs_proc], /get_lun
         printf, pbs_ppn_lun, "#done"
         close, pbs_ppn_lun
         free_lun, pbs_ppn_lun
       endfor
       pbs_node_i = 0
     endfor
   endif else print, 'completed'
   
   ; Close the bundled script files if pbs_nodes keyword is set
   if keyword_set(pbs_nodes) then begin
     for pbs_node = 0, pbs_nodes-1 do begin
       printf, pbs_node_lun[pbs_node], 'wait'
       close, pbs_node_lun[pbs_node]
       free_lun, pbs_node_lun[pbs_node]
     endfor
   endif

   ;----------
   ; Trim the plate list to only those needing reductions

   ibatch = where(qbatch, nbatch)
   if (nbatch EQ 0) then begin
      splog, 'All plates have been reduced'
      return
   endif

   ;----------
   ; Submit jobs to the PBS queue

   if (not keyword_set(pbs_nodes)) then begin
     for i=0L, nbatch-1L do begin
        thisfile = fileandpath(fullscriptfile[ibatch[i]], path=thispath)
        if (keyword_set(thispath)) then cd, thispath
        spawn, 'qsub '+thisfile
     endfor
   endif else begin
     cd, pbs_dir
     for i=0L, pbs_nodes-1 do spawn, 'qsub '+pbs_node_script[i]
   endelse

   return
end
;------------------------------------------------------------------------------
