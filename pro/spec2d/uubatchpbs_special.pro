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
;    mjd=, mjstart=, mjend=, upsvers2d=, upsvers1d=, upsversutils=, rawdata_dir=, $
;    boss_spectro_redux=, scratchdir=, /zcode, /galaxy, upsversgalaxy=, boss_galaxy_redux=, boss_galaxy_scratch=, $
;    pbsdir=, /verbose, queue=, qos=, ebossvers=ebossvers, /daily, /skip2d, $
;    /skip_granada_fsps, /skip_portsmouth_stellarmass, /skip_portsmouth_emlinekin, /skip_wisconsin_pca, $
;   /clobber, /nosubmit, /test, /slurm, $
;    pbsnodes=pbsnodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, pbs_walltime=pbs_walltime, /pbs_batch, /riemann, /ember, /kingspeak]
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
;   upsvers2d  - If set, then do a "module switch idlspec2d idlspec2d$UPSVERS2D" on the
;                remote machine before executing Spectro-2D.  This allows
;                you to batch jobs using a version other than that which
;                is declared current under UPS.
;   upsvers1d  - If set, then do a "module switch idlspec2d idlspec2d$UPSVERS2D" on the
;                remote machine before executing Spectro-1D.  This allows
;                you to batch jobs using a version other than that which
;                is declared current under UPS.
;   upsversutils - If set, then do a "module switch idlutils idlutils/$IDLUTILS" on the
;                remote machine.
;   scratchdir   - If set, then treat this as topdir until the computation is complete
;   zcode      - If set, run Zcode in auto mode.
;   galaxy     - If set, run Galaxy (Granada, Portsmouth, Wisconsin) Suite of Products.
;   upsversgalaxy  - If set, then do a "module switch galaxy galaxy/$GALAXY" on the
;                remote machine.
;   boss_galaxy_redux   - Optional override value for the environment variable $BOSS_GALAXY_REDUX,
;   boss_galaxy_scratch - Optional override value for the environment variable $GALAXY_SCRATCH_DIR.
;   pbsdir     - Optional override value for the environment variable $PBS_SCRATCH_DIR.
;   verbose    - If set, then add "set -o verbose" for easier debugging.
;   queue      - If set, sets the submit queue.
;   qos        - If set, sets the submit qos.
;   ebossvers  - If set, switch to this version of the eboss module
;   daily      - If set, switch to the daily module (overrides ebossvers)
;   skip2d     - If set, then skip the Spectro-2D reductions.
;   skip_wisconsin_pca   - If galaxy set and if not set [default], then run wisconsin_pca code
;   skip_granada_fsps - If galaxy set and if not set [default], then run granada_fsps code
;   skip_portsmouth_stellarmass   - If galaxy set and if not set [default], then run portsmouth_stellarmass 
;   skip_portsmouth_emlinekin   - If galaxy set and if not set [default], then run portsmouth_emlinekin
;   clobber    - If set, then reduce all specified plates.  The default is
;                to not reduce plates where the script file already exists.
;   pbs_batch  - If set, collect the pbs_nodes into an PBS array
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
;                pbs_nodes = 12 (for 12 nodes, without node sharing) 
;                pbs_ppn = 12 (12 processors per node)
;                pbs_a = 'bolton-em' (Bolton's account, limited to 12 nodes: ember253-260,377-380)
;                pbs_walltime='336:00:00'
;                slurm=1
;   kingspeak  -If set, then setup the defaults for the kingspeak cluster at the University of Utah:
;                pbs_nodes = 27 (for 27 nodes, without node sharing)
;                pbs_ppn = 16 (16 processors per node)
;                pbs_a = 'sdss-kp' (SDSS account, limited to 27 nodes)
;                pbs_walltime='336:00:00'
;                slurm=1
;   riemann    - If set, then setup the defaults for the riemann cluster at LBL:
;                pbs_nodes = 12 (for 12 nodes, without node sharing) 
;                pbs_ppn = 8 (8 processors per node)
;                pbs_walltime='336:00:00'
;                slurm=0
;   nosubmit   - If set, generate script file but don't submit to queue
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
;                and with LBL defaults preset via the keyword riemann.
;                and with University of Utah defaults preset via the keyword ember or kingspeak.
;   02-Apr-2015  Modified by Joel R. Brownstein, University of Utah, in order to refactor the
;                #PBS directives into a method (uubatchpbs_directives) in order to add a keyword
;                to switch directives to #SBATCH for slurm (Simple Linux Utility for Resource Management).
;   05-Sep-2015  6-digit plate number armaggedon fixed by Julian Bautista, University of Utah
;   30-Jun-2020  Corrected printed file names by Sean Morrison
;-
;------------------------------------------------------------------------------

pro uubatchpbs_directives, pbs_batch_lun=pbs_batch_lun, slurm=slurm, pbs_dir=pbs_dir, pbs_nodes=pbs_nodes, pbs_ppn=pbs_ppn, pbs_mem_per_cpu=pbs_mem_per_cpu, pbs_a=pbs_a, pbs_walltime=pbs_walltime, qos=qos, batch_array=batch_array,pbs_cpus=pbs_cpus,pbs_share=pbs_share, nplate=nplate, full_run=full_run

   if keyword_set(slurm) then begin
        printf, pbs_batch_lun, '#!/bin/bash'
        printf, pbs_batch_lun, '# Auto-generated by uubatchpbs.pro '+systime()
        if keyword_set(pbs_a) then begin
           account = pbs_a
           if keyword_set(qos) then account = account + strtrim((qos eq 'sdss-fast') ? '-fast' : qos,2)
           printf, pbs_batch_lun, '#SBATCH --account='+account
           if keyword_set(pbs_share) then begin 
             temp_strin=strsplit(pbs_a, /EXTRACT, '-') 
             partition = temp_strin[0]+'-shared-'+temp_strin[1]
           endif else begin 
             partition = pbs_a 
           endelse 
           printf, pbs_batch_lun, '#SBATCH --partition='+partition
        endif
        printf, pbs_batch_lun, '#SBATCH --nodes=1'
        if keyword_set(pbs_share) then begin
           if keyword_set(pbs_ppn) then printf, pbs_batch_lun, '#SBATCH --ntasks='+strtrim(min([nplate,pbs_ppn]),2)
        endif else begin
           if keyword_set(pbs_ppn) then printf, pbs_batch_lun, '#SBATCH --ntasks='+strtrim(pbs_ppn,2)
        endelse
        if keyword_set(pbs_mem_per_cpu) then  printf, pbs_batch_lun, '#SBATCH --mem-per-cpu='+strtrim(pbs_mem_per_cpu,2)
        ;printf, pbs_batch_lun, '#SBATCH --job-name=uubatch'
        printf, pbs_batch_lun, '#SBATCH --job-name=uubatch'
        if keyword_set(pbs_walltime) then printf, pbs_batch_lun, '#SBATCH --time='+pbs_walltime
        if keyword_set(batch_array) then begin
           printf, pbs_batch_lun, '#SBATCH --array=1-'+strtrim(pbs_nodes,2)
           printf, pbs_batch_lun, '#SBATCH --output=uubatch_%A[%a].out'
           printf, pbs_batch_lun, '#SBATCH --err=uubatch_%A[%a].err'
        endif else begin
           printf, pbs_batch_lun, '#SBATCH --output=uubatch_%j.out'
           printf, pbs_batch_lun, '#SBATCH --err=uubatch_%j.err'
           if keyword_set(full_run) then begin
             printf, pbs_batch_lun, 'export OMP_NUM_THREADS=1'
             printf, pbs_batch_lun, 'export OPENBLAS_NUM_THREADS=1'
             printf, pbs_batch_lun, 'export MKL_NUM_THREADS=1'
             printf, pbs_batch_lun, 'export VECLIB_MAXIMUM_THREADS=1'
             printf, pbs_batch_lun, 'export NUMEXPR_NUM_THREADS=1'
             printf, pbs_batch_lun, 'export TPOOL_NTHREADS = 1'
           endif
        endelse
        printf, pbs_batch_lun, 'umask 0027'
        if (keyword_set(daily)) then printf, pbs_batch_lun, 'module switch eboss eboss/daily' $
        else if (keyword_set(ebossvers)) then printf, pbs_batch_lun, 'module switch eboss eboss/'+strtrim(ebossvers,2)
        if keyword_set(riemann) and keyword_set(galaxy) and not keyword_set(skip_portsmouth_stellarmass) then printf, pbs_batch_lun, 'source /home/boss/.intel64'
        if keyword_set(batch_array) then begin
           printf, pbs_batch_lun, 'SBATCH_NODE=$( printf "%02d\n" "$SLURM_ARRAY_TASK_ID" )'
           printf, pbs_batch_lun, 'source '+pbs_dir+'node${SBATCH_NODE}.pbs'
        endif
   endif else begin
        printf, pbs_batch_lun, '#!/bin/bash'
        printf, pbs_batch_lun, '# Auto-generated by uubatchpbs.pro '+systime()
        if keyword_set(pbs_a) then printf, pbs_batch_lun, '#PBS -A '+pbs_a
        if keyword_set(pbs_walltime) then printf, pbs_batch_lun, '#PBS -l walltime='+pbs_walltime
        printf, pbs_batch_lun, '#PBS -W umask=0027'
        printf, pbs_batch_lun, '#PBS -V'
        printf, pbs_batch_lun, '#PBS -j oe'
        printf, pbs_batch_lun, '#PBS -t 1-'+strtrim(pbs_nodes,2)
        printf, pbs_batch_lun, '#PBS -N uubatch'
        if (keyword_set(queue)) then printf, pbs_batch_lun, '#PBS -q ' + queue
        if (keyword_set(qos)) then printf, pbs_batch_lun, '#PBS -l qos=' + qos
        if keyword_set(pbs_ppn) then printf, pbs_batch_lun, '#PBS -l nodes=1:ppn='+strtrim(pbs_ppn,2) $
        else printf, pbs_batch_lun, '#PBS -l nodes=1'
        if (keyword_set(daily)) then printf, pbs_batch_lun, 'module switch eboss eboss/daily' $
        else if (keyword_set(ebossvers)) then printf, pbs_batch_lun, 'module switch eboss eboss/'+strtrim(ebossvers,2)
        if keyword_set(riemann) and keyword_set(galaxy) and not keyword_set(skip_portsmouth_stellarmass) then printf, pbs_batch_lun, 'source /home/boss/.intel64'
        if keyword_set(batch_array) then begin
           printf, pbs_batch_lun, 'PBS_JOBID=$( printf "%02d\n" "$PBS_ARRAYID" )'
           printf, pbs_batch_lun, 'source '+pbs_dir+'node${PBS_JOBID}.pbs'
        endif
   endelse
end

;------------------------------------------------------------------------------

pro trim_keywords, keywords
   if not strmatch(strmid(keywords,keywords.LastIndexOf(','),strlen(keywords)-keywords.LastIndexOf(',')),'*[a-z]*') then begin
      keywords=strmid(keywords,0,keywords.LastIndexOf(','))
   endif
   keywords=keywords
end

;------------------------------------------------------------------------------

pro uubatchpbs_special, platenums1, topdir=topdir1, run2d=run2d1, run1d=run1d1, $
 platestart=platestart, plateend=plateend, $
 mjd=mjd, mjstart=mjstart, mjend=mjend, $
 upsvers2d=upsvers2d, upsvers1d=upsvers1d, upsversutils=upsversutils, $
 rawdata_dir=rawdata_dir, slurm=slurm, $
 boss_spectro_redux=boss_spectro_redux, scratchdir=scratchdir, $
 zcode=zcode, galaxy=galaxy, upsversgalaxy=upsversgalaxy, pbsdir=pbsdir, $
 boss_galaxy_redux=boss_galaxy_redux, boss_galaxy_scratch=boss_galaxy_scratch, $
 verbose=verbose, queue=queue, qos=qos, ebossvers=ebossvers, daily=daily, skip2d=skip2d, clobber=clobber, $
 nosubmit=nosubmit, test=test, $
 skip_granada_fsps=skip_granada_fsps, skip_portsmouth_stellarmass=skip_portsmouth_stellarmass, $
 skip_portsmouth_emlinekin=skip_portsmouth_emlinekin, skip_wisconsin_pca=skip_wisconsin_pca,  $
 pbs_nodes=pbs_nodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, pbs_batch=pbs_batch, MWM_fluxer=MWM_fluxer, $
 noxcsao=noxcsao, $
 pbs_walltime=pbs_walltime, lco=lco, plate_s=plate_s, legacy=legacy, full_run=full_run, _EXTRA=Extra
;, riemann=riemann, ember=ember, kingspeak=kingspeak


   if (size(platenums1,/tname) EQ 'STRING') then platenums = platenums1 $
    else if (keyword_set(platenums1)) then $
      platenums = plate_to_string(platenums1) $
    else platenums = '*'
   
   ;----------
   ; Determine the top-level of the output directory tree
   if keyword_set(lco) then begin
     obsdir='LCO'
   endif else begin
     obsdir='APO'
   endelse
   obsdir='';coment this line for the final version HJIM
   if (keyword_set(topdir1)) then topdir = topdir1 $
   else begin
     topdir = getenv('BOSS_SPECTRO_REDUX')
     ;topdir=concat_dir(topdir, obsdir)
     if strpos(topdir,'/',strlen(topdir)-1) lt 0 then topdir+='/'
     if keyword_set(test) and not (strlen(topdir)-rstrpos(dir,'/test/') eq strlen('/test/')) then topdir=djs_filepath('',root_dir=topdir, subdir='test')
   endelse
   if strpos(topdir,'/',strlen(topdir)-1) lt 0 then topdir+='/'
   splog, 'Setting TOPDIR=', topdir

   
   if keyword_set(galaxy) then begin
       if (keyword_set(boss_galaxy_redux)) then boss_galaxy_redux = strtrim(boss_galaxy_redux,2) else begin
         boss_galaxy_redux = getenv('BOSS_GALAXY_REDUX')
         if strpos(boss_galaxy_redux,'/',strlen(boss_galaxy_redux)-1) lt 0 then boss_galaxy_redux+='/'
         if keyword_set(test) and not (strlen(boss_galaxy_redux)-rstrpos(boss_galaxy_redux,'/test/') eq strlen('/test/')) then boss_galaxy_redux=djs_filepath('',root_dir=boss_galaxy_redux, subdir='test')
       endelse
       if strpos(boss_galaxy_redux,'/',strlen(boss_galaxy_redux)-1) lt 0 then boss_galaxy_redux+='/'
       splog, 'Setting BOSS_GALAXY_REDUX=', boss_galaxy_redux

       if (keyword_set(boss_galaxy_scratch)) then boss_galaxy_scratch = strtrim(boss_galaxy_scratch,2) else begin
         boss_galaxy_scratch = getenv('GALAXY_SCRATCH_DIR')
         if strpos(boss_galaxy_scratch,'/',strlen(boss_galaxy_scratch)-1) lt 0 then boss_galaxy_scratch+='/'
         if keyword_set(test) and not (strlen(boss_galaxy_scratch)-rstrpos(boss_galaxy_scratch,'/test/') eq strlen('/test/')) then boss_galaxy_scratch=djs_filepath('',root_dir=boss_galaxy_scratch, subdir='test')
       endelse
       if strpos(boss_galaxy_scratch,'/',strlen(boss_galaxy_scratch)-1) lt 0 then boss_galaxy_scratch+='/'
       splog, 'Setting GALAXY_SCRATCH_DIR=', boss_galaxy_scratch
   endif

   if (keyword_set(run2d1)) then run2d = strtrim(run2d1,2) $
    else run2d = getenv('RUN2D')
   splog, 'Setting RUN2D=', run2d
   if (keyword_set(run1d1)) then run1d = strtrim(run1d1,2) $
    else run1d = getenv('RUN1D')
   splog, 'Setting RUN1D=', run1d
   if (keyword_set(upsvers1d)) then splog, 'Setting IDLSPEC2D=', upsvers1d $
   else if (keyword_set(upsvers2d)) then splog, 'Setting IDLSPEC2D=', upsvers2d
   if (keyword_set(upsversutils)) then splog, 'Setting IDLUTILS=', upsversutils
   if (keyword_set(upsversgalaxy)) then splog, 'Setting GALAXY=', upsversgalaxy
   
    if (keyword_set(pbsdir)) then pbsdir = strtrim(pbsdir,2) else begin
        ;pbsdir = getenv('PBS_SCRATCH_DIR')
        pbsdir = getenv('SLURM_SCRATCH_DIR')
        if strpos(pbsdir,'/',strlen(pbsdir)-1) lt 0 then pbsdir+='/'
        if keyword_set(test) and not (strlen(pbsdir)-rstrpos(pbsdir,'/test/') eq strlen('/test/')) then pbsdir=djs_filepath('',root_dir=pbsdir, subdir='test')
    endelse
    if strpos(pbsdir,'/',strlen(pbsdir)-1) lt 0 then pbsdir+='/'
    splog, 'Setting PBS_SCRATCH_DIR=', pbsdir

   topdir2d = djs_filepath('', root_dir=topdir, subdir=run2d)

   if (keyword_set(run1d)) then begin
        run1dstr = ', run1d="'+run1d+'"'
        run1dstr_py = '--run1d "'+run1d+'"'
   endif else begin
        run1dstr = ''
        run1dstr_py = ''
   endelse
   if (keyword_set(run2d)) then run2dstr = ',run2d="'+run2d+'"' $
    else run2dstr = ''

   ;if keyword_set(riemann) then begin
    ; if not keyword_set(pbs_nodes) then pbs_nodes=28
    ; if not keyword_set(pbs_ppn) then pbs_ppn=8
    ; if not keyword_set(pbs_walltime) then pbs_walltime='336:00:00'
    ; slurm = 1
   ;endif else if keyword_set(ember) then begin
   ;  if not keyword_set(pbs_nodes) then pbs_nodes=12
   ;  if not keyword_set(pbs_ppn) then pbs_ppn=12
   ;  if not keyword_set(pbs_walltime) then pbs_walltime='336:00:00'
   ;  if not keyword_set(pbs_a) then pbs_a = 'bolton-em'
   ;  slurm = 1
   ;endif else if keyword_set(kingspeak) then begin
   ;  if not keyword_set(pbs_nodes) then pbs_nodes=getenv('SLURM_NODES')
   ;  if not keyword_set(pbs_ppn) then pbs_ppn=getenv('SLURM_PPN')
   ;  if not keyword_set(pbs_walltime) then pbs_walltime='336:00:00'
   ;  if not keyword_set(pbs_a) then pbs_a = getenv('SLURM_ALLOC')
   ;  slurm = 1
   pbs_share=0
   if keyword_set(slurm) then begin
;     if not keyword_set(pbs_nodes) then pbs_nodes=getenv('SLURM_NODES')
     if not keyword_set(pbs_nodes) then pbs_nodes=getenv('SLURM_NODES')
     if not keyword_set(pbs_ppn) then pbs_ppn=getenv('SLURM_PPN')
     if not keyword_set(pbs_walltime) then pbs_walltime='336:00:00'
     if not keyword_set(pbs_a) then pbs_a = getenv('SLURM_ALLOC')
     if not keyword_set(pbs_mem_per_cpu) then pbs_mem_per_cpu=getenv('SLURM_MEM_PER_CPU')
     if not keyword_set(pbs_share) then begin
       if getenv('SLURM_VERS') eq 'notchpeak' then pbs_share=1
     endif
   endif

   sub = keyword_set(slurm) ? 'sbatch' : 'qsub'

   if ((keyword_set(slurm)) and (pbs_nodes gt 1)) then pbs_batch = 1L

   ;----------
   ; Create list of plate directories
   ; Limit the list to only those specified by PLATENUMS,PLATESTART,PLATEEND
;   if keyword_set(plate_s) then begin
;     platedirs = get_mjd_dir(topdir2d, mjd=platenums, mjstart=platestart, $
;       mjend=plateend,/alldirs)
;     for ili=0, n_elements(platedirs)-1 do begin
;       ;print,strmid(strtrim(platedirs[ili],2),5,1)
;       ;print,strmid(strtrim(platedirs[ili],2),4,1)
;       if strmid(strtrim(platedirs[ili],2),4,1) ne 'p' then begin
;         if strmid(strtrim(platedirs[ili],2),5,1) ne 'p' then begin
;           platedirs[ili]=''
;         endif
;       endif
;     endfor
;     ii = where(platedirs NE '', ct)
;     if (ct EQ 0) then begin
;       splog, 'No plate directories found'
;       return
;     endif else begin
;       platedirs = platedirs[ii]
;     endelse
;   endif else begin

    platedirs = get_mjd_dir(topdir2d, mjd=platenums, mjstart=platestart, $
      mjend=plateend)
;   endelse
  
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
         temp_plan=strsplit(planfile[ifile],'/',/extract)
         temp_mjd=strsplit(temp_plan[n_elements(temp_plan)-1],'-',/extract)
         temp_mjd=long(repstr(temp_mjd[n_elements(temp_mjd)-1],'.par'))
         if keyword_set(plate_s) and not keyword_set(legacy) then begin
           min_mjd=59030;LIMIT the use of single spectrograph after mjd 59005
           max_mjd=70000; Needs to change
         endif else begin
           if keyword_set(legacy) then begin
             min_mjd=0
             max_mjd=59030
           endif else begin
             min_mjd=0
             max_mjd=70000
           endelse
         endelse
         if temp_mjd ge min_mjd and temp_mjd lt max_mjd then begin
         ;print,planfile[ifile]
         yanny_read, planfile[ifile], hdr=hdr
         thismjd = long(yanny_par(hdr, 'MJD'))
         ; Decide if THISMJD is within the bounds specified by MJD,MJSTART,MJEND
         ;if (mjd_match(thismjd, mjd=mjd, mjstart=mjstart, mjend=mjend)) then begin
            if (keyword_set(platelist)) then begin
               platelist = [platelist, platedirs[idir]]
               planlist = [planlist, planfile[ifile]]
            endif else begin
               platelist = platedirs[idir]
               planlist = planfile[ifile]
            endelse
         ;endif
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

     if keyword_set(pbs_ppn) then nodes_required = ceil(float(nplate)/pbs_ppn) else nodes_required = nplate
     if nodes_required lt pbs_nodes then pbs_nodes = nodes_required
     ncycle = ceil(float(nplate)/(pbs_nodes*pbs_ppn))
     if pbs_share eq 1 then begin
       if nplate le pbs_ppn then pbs_ncpus=nplate else pbs_ncpus = pbs_ppn
     endif else begin
       pbs_ncpus = pbs_ppn
     endelse

     home = getenv('HOME')
     if (home ne '') then begin
      pos0 = strpos(home,'/',/reverse_search)+1
      pos1 = strlen(home)
      userID = strmid(home,pos0,pos1-pos0)
     endif else userID = 'user'
     
     date = strtrim(bin_date(),2)
     for i=1,n_elements(date)-1 do $
        if (strlen(date[i]) eq 1) then date[i] = '0'+date[i] 

     splog, 'Starting for user:  ',userID
     userID+='_'+string(date, format='(A4,A2,A2,A2,A2,A2)')

     if (pbsdir eq '') then begin
        if keyword_set(scratchdir) then pbsdir = djs_filepath('pbs/'+run2d,root_dir=scratchdir) $
        else pbsdir = djs_filepath('pbs/'+run2d,root_dir=topdir)
     endif else pbs_dir = djs_filepath('bossredux/'+run2d,root_dir=pbsdir)
     pbs_dir = djs_filepath('',root_dir=pbsdir,subdir=run2d+'/'+userID)
     if file_test(pbs_dir) then begin
       shift_pbs_dir = djs_filepath('',root_dir=pbsdir,subdir=run2d+'/'+userID+'.*')
       shift_pbs = file_search(shift_pbs_dir, count=nshift_pbs)
       max_shift = -1L
       for i=0,nshift_pbs-1 do begin
         pos0 = strpos(shift_pbs[i],'/'+userID+'.',/reverse_search)+strlen('/'+userID+'.')
         pos1 = strlen(shift_pbs[i])
         next_shift = fix(strmid(shift_pbs[i],pos0,pos1-pos0))
         max_shift = (next_shift gt max_shift) ? next_shift : max_shift
       endfor
       shift_pbs_dir = djs_filepath('',root_dir=pbsdir,subdir=run2d+'/'+userID+ '.' + strtrim(max_shift+1,2))
       splog, 'Renaming previous PBS directory to: '+shift_pbs_dir
       file_move, pbs_dir, shift_pbs_dir
       file_mkdir, pbs_dir
     endif else file_mkdir, pbs_dir
     splog, "cd "+pbs_dir
         
     pbs_node_index = 'node'+ string(indgen(pbs_nodes)+1,format='(i2.2)')
     pbs_node_script = djs_filepath(pbs_node_index+'.pbs',root_dir=pbs_dir)
     pbs_node_lun = intarr(pbs_nodes)
     if keyword_set(pbs_ppn) then begin
       pbs_ppn_index  = '_proc'+ string(indgen(pbs_ppn)+1,format='(i2.2)')
       pbs_ppn_script = strarr(pbs_nodes,pbs_ppn)
     endif
     
     if keyword_set(pbs_ppn) then splog, 'Preparing to '+sub+' '+strtrim(pbs_nodes,2)+' nodes ('+strtrim(pbs_ppn,2)+' processors per node) for '+strtrim(nplate,2)+' plates.' $
     else splog, 'Preparing to '+sub+' '+strtrim(nplate,2)+' plates.'

     for pbs_node = 0, pbs_nodes-1 do begin
       openw, get_lun, pbs_node_script[pbs_node] ,/get_lun
       pbs_node_lun[pbs_node] = get_lun
       uubatchpbs_directives, pbs_batch_lun=pbs_node_lun[pbs_node], slurm=slurm, pbs_nodes=pbs_nodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, pbs_walltime=pbs_walltime, qos=qos, $
       pbs_cpus=pbs_cpus, pbs_share=pbs_share,full_run=full_run, nplate=nplate
       if keyword_set(pbs_ppn) then begin
         pbs_ppn_script[pbs_node,*] = djs_filepath(pbs_node_index[pbs_node] + pbs_ppn_index +'.pbs',root_dir=pbs_dir)
         for pbs_proc = 0, pbs_ppn-1 do printf, pbs_node_lun[pbs_node], 'source '+pbs_ppn_script[pbs_node,pbs_proc] + ' &'
       endif
       close, pbs_node_lun[pbs_node]
     endfor 

     if keyword_set(pbs_batch) then begin
        pbs_batch_script = djs_filepath('uubatch.pbs',root_dir=pbs_dir)
        openw, pbs_batch_lun, pbs_batch_script, /get_lun
        uubatchpbs_directives, pbs_batch_lun=pbs_batch_lun, slurm=slurm, pbs_dir=pbs_dir, pbs_nodes=pbs_nodes, pbs_ppn=pbs_ppn, pbs_a=pbs_a, pbs_walltime=pbs_walltime, $
        pbs_share=pbs_share, qos=qos, nplate=nplate,full_run=full_run, /batch_array
        close, pbs_batch_lun
        free_lun, pbs_batch_lun
     endif
     
   endif
   pbs_node = 0
   pbs_proc = 0
   pbs_ppn_append = 0
   
   if keyword_set(galaxy) then begin
       n_redux = 5
       galaxy_redux = replicate({group:'', product:'', file:'', counter:0L, done:0L, skip:0B, keyword:''},n_redux)
       galaxy_redux.group = ['granada','portsmouth','portsmouth','utah','wisconsin']
       galaxy_redux.product = ['fsps','stellarmass','emlinekin','bells','pca']
       galaxy_redux.skip = [keyword_set(skip_granada_fsps),keyword_set(skip_portsmouth_stellarmass),keyword_set(skip_portsmouth_emlinekin),keyword_set(skip_utah_bells),keyword_set(skip_wisconsin_pca)]
       galaxy_redux.keyword = ['/skip_granada_fsps','/skip_portsmouth_stellarmass','/skip_portsmouth_emlinekin','/skip_utah_bells','/skip_wisconsin_pca']
   endif

   cycle = 0
   for iplate=0, nplate-1 do begin
      ; Find all relevant 2D plan files
      yanny_read, planlist[iplate], hdr=hdr
      planfile2d = yanny_par(hdr, 'planfile2d')
      ;print,planfile2d
      if keyword_set(plate_s) then begin
         plateid[iplate] = yanny_par(hdr, 'plateid')
      endif else begin
         plateid[iplate] = yanny_par(hdr, 'fieldid')
      endelse
      mjd = yanny_par(hdr, 'MJD')
      platemjd = field_to_string(plateid[iplate]) + '-' $
       + string(mjd,format='(i5.5)')
      platefile = 'spField-'+platemjd+'.fits'

      ; Track the beginning and ending MJD going into this plate
      mjd_beg[iplate] = min(strmid(planfile2d,14,5))
      mjd_end[iplate] = max(strmid(planfile2d,14,5))

      ; Split the combine plan file name into a directory and file name
      planfilecomb = fileandpath(planlist[iplate], path=pathcomb)
           
      if keyword_set(scratchdir) then begin
        scratchdir2d = djs_filepath(field_to_string(plateid[iplate]), root_dir=scratchdir, subdir=run2d)
        scratchdir1d = djs_filepath('', root_dir=scratchdir2d, subdir=run1d)
        fullscriptfile[iplate] = djs_filepath('redux-'+platemjd, root_dir=scratchdir2d)
        redux_file = djs_filepath('redux-'+platemjd, root_dir=topdir2d,subdir=field_to_string(plateid[iplate]))
      endif else begin
        fullscriptfile[iplate] = djs_filepath('redux-'+platemjd, root_dir=pathcomb)
        redux_file = fullscriptfile[iplate]
      endelse 
      if (keyword_set(skip2d)) then fullscriptfile[iplate] += '-' + run1d
      
      if keyword_set(clobber) then file_delete, redux_file, /quiet, /allow_nonexistent $
      else qbatch[iplate] = file_test(redux_file) ? 0B : 1B

      if keyword_set(galaxy) then begin
        galaxy_outdir =  djs_filepath('', root_dir=boss_galaxy_redux, subdir=run2d+'/'+field_to_string(plateid[iplate])+'/'+run1d)

        galaxy_redux_file = strarr(n_redux)
        for r=0,n_redux-1 do galaxy_redux_file[r] = djs_filepath(galaxy_redux[r].group + '_' + galaxy_redux[r].product + '_redux-'+ platemjd, root_dir=galaxy_outdir, subdir=galaxy_redux[r].group + '/' + galaxy_redux[r].product)

        keywords = ''
        for r=0,n_redux-1 do begin
          if keyword_set(galaxy_redux[r].skip) then begin
              keywords += ', '+galaxy_redux[r].keyword
              if file_test(galaxy_redux_file[r]) then galaxy_redux[r].done += 1L
          endif else begin
              galaxy_redux[r].counter += 1L
              if keyword_set(clobber) then file_delete, galaxy_redux_file[r], /quiet, /allow_nonexistent
          endelse
        endfor
        if keyword_set(test) then keywords += ', /test'
      endif

      if (qbatch[iplate]) then begin

        if keyword_set(scratchdir) then begin
        
          ; Construct run2d and run1d directories for each plate within scratchdir
          if (not file_test(scratchdir2d)) then file_mkdir, scratchdir2d
          
          ; cp the plan files to scratch if needed:
          file_copy, planlist[iplate], scratchdir2d, /over
          planfile2d_source = file_search(djs_filepath(planfile2d,root_dir=topdir2d,subdir=field_to_string(plateid[iplate])),count=has_plan2d)
          if keyword_set(has_plan2d) then file_copy, planfile2d_source, scratchdir2d, /over
        endif

         openw, olun, fullscriptfile[iplate], /get_lun
         printf, olun, '# Auto-generated batch file '+systime()
         if not keyword_set(pbs_nodes) then begin
           if keyword_set(pbs_ppn) then printf, olun, '#PBS -l nodes=1:ppn='+strtrim(pbs_ppn,2) $
            else printf, olun, '#PBS -l nodes=1'
           if keyword_set(pbs_a) then printf, olun, '#PBS -A '+pbs_a
           if keyword_set(pbs_walltime) then printf, olun, '#PBS -l walltime='+pbs_walltime
           printf, olun, '#PBS -W umask=0027'
           printf, olun, '#PBS -V'
           printf, olun, '#PBS -j oe'
           ; set queue if asked
           if (keyword_set(queue)) then printf, olun, '#PBS -q ' + queue
           printf, olun, 'cd $PBS_O_WORKDIR'
         endif else begin
            if keyword_set(scratchdir) then  printf, olun, 'cd '+scratchdir2d $
            else printf, olun, 'cd '+pathcomb
         endelse 

         ; Define the observatory data
         if (keyword_set(lco)) then begin
             printf, olun, 'export BOSS_SPECTRO_DATA=$BOSS_SPECTRO_DATA_S'
         endif else begin
             printf, olun, 'export BOSS_SPECTRO_DATA=$BOSS_SPECTRO_DATA_N'
         endelse
         ; Override environment variables if requested
         if (keyword_set(rawdata_dir)) then begin
             printf, olun, 'export BOSS_SPECTRO_DATA='+rawdata_dir
         endif
         if (keyword_set(boss_spectro_redux)) then begin
             printf, olun, 'export BOSS_SPECTRO_REDUX='+boss_spectro_redux
         endif
         if (keyword_set(run2d)) then begin
             printf, olun, 'export RUN2D='+run2d
         endif
         if (keyword_set(run1d)) then begin
             printf, olun, 'export RUN1D='+run1d
         endif
         if keyword_set(verbose) then begin
             printf, olun, ''
             printf, olun, '#- Echo commands to make debugging easier'
             printf, olun, 'set -o verbose'
             printf, olun, ''
             printf, olun, '#- The real work'
         endif

        ; Delete old files in clobber keyword is set - Added by TH 5 Aug 2015
         if keyword_set(clobber) then begin
            file_delete, 'spec2d-'+platemjd+'.started', /quiet, /allow_nonexistent
            file_delete, 'spec2d-'+platemjd+'.done', /quiet, /allow_nonexistent
            file_delete, 'spec1d-'+platemjd+'.started', /quiet, /allow_nonexistent
            file_delete, 'spec1d-'+platemjd+'.done', /quiet, /allow_nonexistent
            file_delete, 'redmonster-'+platemjd+'.started', /quiet, /allow_nonexistent
            file_delete, 'redmonster-'+platemjd+'.done', /quiet, /allow_nonexistent
         endif

         if (keyword_set(skip2d) EQ 0) then begin
            ; Set up requested code version
            if (keyword_set(upsvers2d)) then $
             printf, olun, 'module switch idlspec2d idlspec2d/'+upsvers2d
            if (keyword_set(upsversutils)) then printf, olun, 'module switch idlutils idlutils/'+upsversutils

            ; Create sorted photoPlate files
            ;for i=0, n_elements(planfile2d)-1 do begin
            ; if keyword_set(plate_s) then begin
            ;   printf, olun, 'echo '+fq+'sdss_plate_sort,"'+planfile2d[i]+'"'+fq+' | idl'
            ; endif else begin
            ;   printf, olun, 'echo '+fq+'sdss_field_sort,"'+planfile2d[i]+'"'+fq+' | idl'
            ; endelse
            ;endfor
            ; Run Spectro-2D
            
            ;spreduce2d_keys=''
            rm_combine_keys='/xyfit, /loaddesi,'
            if keyword_set(plate_s) then begin
                if keyword_set(legacy) then begin
                    ;spreduce2d_keys = spreduce2d_keys +' /legacy,'
                    rm_combine_keys = rm_combine_keys +' /legacy,'
                endif else begin
                    ;spreduce2d_keys = spreduce2d_keys +' /plates,'
                    rm_combine_keys = rm_combine_keys +' /plates,'
                    if keyword_set(MWM_fluxer) then $
                        ;spreduce2d_keys = spreduce2d_keys +' /MWM_fluxer,'
                        rm_combine_keys = rm_combine_keys +' /MWM_fluxer,'
                endelse
            endif
            
            ;for i=0, n_elements(planfile2d)-1 do begin
            ;    pmjd=STRJOIN((STRSPLIT((STRSPLIT(planfile2d[i],'.',/EXTRACT))[0],'-',/extract))[1:*],'-')
            ;    ;printf, olun, 'touch spec2d-'+platemjd+'.started'       ; Added TH 4 Aug 2015
            ;    printf, olun, 'touch spec2d-'+pmjd+'.started'       ; Added TH 4 Aug 2015 ; updated SM 40 Jun 2020
            ;    printf, olun, 'echo '+fq+'spreduce2d, '+spreduce2d_keys+' "'+planfile2d[i]+'"'+fq+' | idl'
            ;    ;printf, olun, 'touch spec2d-'+platemjd+'.done'          ; Added TH 4 Aug 2015
            ;    printf, olun, 'touch spec2d-'+pmjd+'.done'       ; Added TH 4 Aug 2015 ; updated SM 40 Jun 2020
            ;endfor
            ;printf, olun, 'echo '+fq+'spcombine_v5,"'+planfilecomb+'"'+fq+' | idl'
            printf, olun, 'touch specombine-'+platemjd+'.started'       ; Added HI 21 Nov 2018
            ;printf, olun, 'echo '+fq+'spcombine_v5,"'+planfilecomb+'",minsn2=0.0'+fq+' | idl'
            ;printf, olun, 'echo '+fq+'spcombine_v5,"'+planfilecomb+'",minsn2=0.0'+fq+' | idl'
            printf, olun, 'echo '+fq+'rm_combine_script,"'+planfilecomb+'", '+rm_combine_keys+' bscore=0.01, run2d="'+run2d+'"'+fq+' | idl'
            printf, olun, 'touch specombine-'+platemjd+'.done'       ; Added HI 21 Nov 2018
         endif

         ; Run Spectro-1D
         if (keyword_set(upsvers1d)) then $
          printf, olun, 'module switch idlspec2d idlspec2d/'+upsvers1d
         if (keyword_set(upsversutils)) then printf, olun, 'module switch idlutils idlutils/'+upsversutils
         printf, olun, 'touch spec1d-'+platemjd+'.started'              ; Added TH 4 Aug 2015
         ;printf, olun, 'echo '+fq+'spreduce1d,"'+platefile+'"'+run1dstr+fq+' | idl'
         printf, olun, 'echo '+fq+'spreduce1d_empca,"'+platefile+'"'+run1dstr+fq+' | idl'
         if not keyword_set(noxcsao) then printf, olun, 'run_PyXCSAO.py '+platefile+' '+run1dstr_py
         printf, olun, 'touch spec1d-'+platemjd+'.done'                 ; Added TH 4 Aug 2015
         printf, olun, ''
         printf, olun, '#- Make final spectra files'
         
         reformat_keywords=''
         conflist_keywords='/create, '
         fieldmerge_keywords='/include_bad, '
         if keyword_set(plate_s) then begin
            if keyword_set(legacy) then begin
                reformat_keywords=reformat_keywords+'/legacy, '
                conflist_keywords=conflist_keywords+'/legacy, '
                fieldmerge_keywords=fieldmerge_keywords+'/legacy, '
            endif else begin
                reformat_keywords=reformat_keywords+'/plates, '
                conflist_keywords=conflist_keywords+'/plates, '
                fieldmerge_keywords=fieldmerge_keywords+'/plates, '
            endelse
         endif
         
         if not keyword_set(noxcsao) then begin
            reformat_keywords=reformat_keywords+'/XCSAO, '
            fieldmerge_keywords=fieldmerge_keywords+'/XCSAO, '
         endif
         trim_keywords, fieldmerge_keywords
         fieldmerge_keywords=fieldmerge_keywords+' '+run2dstr+', programs="*eFEDS*", '
         
         trim_keywords, reformat_keywords
         trim_keywords, conflist_keywords
         trim_keywords, fieldmerge_keywords

         printf, olun, 'echo '+fq+'reformat_spec, "'+platefile+'", '+reformat_keywords +run1dstr+' '+run2dstr+fq+' | idl'
         printf, olun, 'echo '+fq+'conflist, '+ conflist_keywords+' '+run1dstr+' '+run2dstr+fq+' | idl'
         printf, olun, 'echo '+fq+'fieldmerge, field='+strtrim(string(plateid[iplate]),2)+', mjd='+strtrim(string(mjd),2)+', '+fieldmerge_keywords+fq+' | idl'
         printf, olun, 'echo '+fq+'reformat_spec, "'+platefile+'", /lite, '+reformat_keywords +run1dstr+' '+run2dstr+fq+' | idl'

;         if keyword_set(plate_s) then begin
;           if keyword_set(legacy) then begin
;             printf, olun, 'echo '+fq+'reformat_spec,"'+platefile+'", /legacy '+run1dstr+fq+' | idl'
;             printf, olun, 'echo '+fq+'conflist, /legacy, /create '+fq+' | idl'
;             printf, olun, 'echo '+fq+'fieldmerge, field='+strtrim(string(plateid[iplate]),2)+', mjd='+strtrim(string(mjd),2)+', /legacy, /include_bad '+fq+' | idl'
;             printf, olun, 'echo '+fq+'reformat_spec,"'+platefile+'", /lite , /legacy '+run1dstr+fq+' | idl'
;           endif else begin
;             printf, olun, 'echo '+fq+'reformat_spec,"'+platefile+'", /plates '+run1dstr+fq+' | idl'
;             printf, olun, 'echo '+fq+'conflist, /plates, /create '+fq+' | idl'
;             printf, olun, 'echo '+fq+'fieldmerge, field='+strtrim(string(plateid[iplate]),2)+', mjd='+strtrim(string(mjd),2)+', /plates, /include_bad '+fq+' | idl'
;             printf, olun, 'echo '+fq+'reformat_spec,"'+platefile+'", /lite , /plates '+run1dstr+fq+' | idl'
;           endelse
;         endif else begin
;           printf, olun, 'echo '+fq+'reformat_spec,"'+platefile+'"'+run1dstr+fq+' | idl'
;           printf, olun, 'echo '+fq+'conflist, /create '+fq+' | idl'
;           printf, olun, 'echo '+fq+'fieldmerge, field='+strtrim(string(plateid[iplate]),2)+', mjd='+strtrim(string(mjd),2)+', /include_bad '+fq+' | idl'
;           printf, olun, 'echo '+fq+'reformat_spec,"'+platefile+'", /lite '+run1dstr+fq+' | idl'
;         endelse
         

         ; Run Zcode
         if (keyword_set(zcode)) then begin
            printf, olun, ''
            printf, olun, 'module load runz'
            printf, olun, 'runz_BOSS.sh ' + platefile +' -a'
            printf, olun, 'runz_BOSS.sh ' + platefile +' -a -G -t GAL'
         endif
         
         ; Run Galaxy Suite of Products
         if (keyword_set(galaxy)) then begin
             if (keyword_set(upsversgalaxy)) then printf, olun, 'module switch galaxy galaxy/'+upsversgalaxy
             printf, olun, 'export BOSS_GALAXY_REDUX='+boss_galaxy_redux
             printf, olun, 'export GALAXY_SCRATCH_DIR='+boss_galaxy_scratch
             skip_keywords = ''
             if keyword_set(skip_granada_fsps) then skip_keywords += ', /skip_granada_fsps'
             if keyword_set(skip_wisconsin_pca) then skip_keywords += ', /skip_wisconsin_pca'
             if keyword_set(skip_portsmouth_stellarmass) then skip_keywords += ', /skip_portsmouth_stellarmass'
             if keyword_set(skip_portsmouth_emlinekin) then skip_keywords += ', /skip_portsmouth_emlinekin'
             for i=0, n_elements(planfile2d)-1 do $
             printf, olun, 'echo '+fq+'galaxy_pipeline,"'+planfile2d[i]+'"'+skip_keywords+fq+' | idl'
         endif


         ; Run reformat_spectra for spec files, full and lite
         




         ; splog, "run1d is ", run1d
         ; splog, "run2d is ", run2d
         
         ; Make pretty pictures
         ;- post-DR9, no longer supported; use spectrawebapp or plotspec instead
          idlcmd  = "plate_spec_image, " + field_to_string(plateid[iplate]);string(plateid[iplate],format='(i4.4)') 
          idlcmd += ", mjd=" + string(mjd,format='(i5.5)')
          idlcmd += ", run1d='" + run1d + "'"
          idlcmd += ", run2d='" + run2d + "'"
          idlcmd += ", /silent"
          if keyword_set(plate_s) then begin
           if keyword_set(legacy) then begin
              idlcmd += ", /legacy"
           endif else begin
              idlcmd += ", /plates"
           endelse
          endif
          printf, olun, ''
          printf, olun, '#- Make pretty pictures'
          printf, olun, 'idl -e "' + idlcmd + '"'
          
          printf, olun, '#- update spAll file'
          printf, olun, 'echo '+fq+'fieldmerge, '+fieldmerge_keywords+fq+' | idl'

;          if keyword_set(plate_s) then begin
;           if keyword_set(legacy) then begin
;             printf, olun, 'echo '+fq+'fieldmerge, /legacy, /include_bad '+fq+' | idl'
;           endif else begin
;             printf, olun, 'echo '+fq+'fieldmerge, /plates, /include_bad '+fq+' | idl'
;           endelse
;          endif else begin
;           printf, olun, 'echo '+fq+'fieldmerge, /include_bad '+fq+' | idl'
;          endelse
          
          ;printf, olun, '#- Make the healpix links'
          ;printf, olun, 'module load sas/main'
          ;printf, olun, 'module switch python miniconda'
          ;iprintf, olun, 'module load miniconda/3.7.7'
          ;printf, olun, 'module load apogee_drp/master'
          ;printf, olun, 'module load dlnpyutils/master'
          ;printf, olun, 'sas_mwm_healpix --spectro boss --mjd '+strtrim(string(mjd),2)+' --telescope apo25m --drpver '+run2d+' -v'
         
         ; If using scratchdir, uubatchcp (selected) final reductions to topdir
         if (keyword_set(scratchdir)) then begin
            for i=0, n_elements(planfile2d)-1 do $
            printf, olun, 'echo '+fq+'uubatchcp,"'+planfile2d[i]+'", topdir="'+topdir+'", run2d="'+run2d+'", run1d="'+run1d+'", scratchdir="'+scratchdir+'"'+fq+' | idl'
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
           if not keyword_set(pbs_node) and not keyword_set(pbs_proc) then begin
             cycle += 1
             splog, "Preparing node cycle "+string(cycle,format='(i2)')+'/'+strtrim(ncycle,2)
           endif
           if keyword_set(pbs_ppn) then begin
             openw, pbs_ppn_lun, pbs_ppn_script[pbs_node,pbs_proc], append=pbs_ppn_append, /get_lun
             printf, pbs_ppn_lun, script_cmd+fullscriptfile[iplate]+' > '+fullscriptfile[iplate]+'.o' $
                                                                 + ' 2> '+fullscriptfile[iplate]+'.e'
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
             openw, pbs_node_lun[pbs_node], pbs_node_script[pbs_node], /append
             printf, pbs_node_lun[pbs_node], script_cmd+fullscriptfile[iplate]+' &'
             close, pbs_node_lun[pbs_node]
             pbs_node += 1
             if pbs_node ge pbs_nodes then pbs_node = 0           
           endelse
         endif
      endif

   endfor
   
   if keyword_set(pbs_ppn) and not keyword_set(pbs_ppn_append) then begin
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
       openw, pbs_node_lun[pbs_node], pbs_node_script[pbs_node], /append
       printf, pbs_node_lun[pbs_node], 'wait'
       printf, pbs_node_lun[pbs_node], 'echo "DONE"'
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

   if keyword_set(pbs_batch) and keyword_set(galaxy) then begin
       galaxy_redux_file = djs_filepath('galaxy_redux.fits',root_dir=pbs_dir)
       splog, "CREATE: "+galaxy_redux_file
       mwrfits, galaxy_redux, galaxy_redux_file, /create, /silent
   endif

   void = where(qbatch eq 1B,n_todo)
   void = where(qbatch eq 0B,n_done)
   message = string("boss_redux",format='(a22)')+': '
   message += "#BOSS PLATES DONE = " + string(n_done,format='(i4)') + "       #BOSS PLATES TO DO = " + string(n_todo,format='(i4)')
   if not keyword_set(clobber) then splog, message $
   else splog, message + ".  [/clobber]"
   if keyword_set(galaxy) then begin
     for r=0,n_redux-1 do begin
       message = string(strtrim(galaxy_redux[r].group,2) + "_" + strtrim(galaxy_redux[r].product,2),format='(a22)')+': '
       message += "#GALAXY JOBS DONE = " + string(galaxy_redux[r].done,format='(i4)') + "       #GALAXY JOBS TO DO = " + string(galaxy_redux[r].counter,format='(i4)')
       if not galaxy_redux[r].skip then splog, message $
       else splog, message + ".  ["+galaxy_redux[r].keyword+"]"
     endfor
   endif

   ;----------
   ; Submit jobs to the PBS queue
   ;exit
   if (not keyword_set(pbs_nodes)) then begin
      for i=0L, nbatch-1L do begin
         thisfile = fileandpath(fullscriptfile[ibatch[i]], path=thispath)
         if (keyword_set(thispath)) then cd, thispath, current=origdir
         if keyword_set(nosubmit) then begin
            splog, 'Generated '+thisfile+' but not submitting to queue'
         endif else begin
            splog, sub+' '+thisfile
            spawn, sub+' '+thisfile, identifier
            splog, strupcase(sub)+'> '+identifier
         endelse
      endfor
      cd, origdir
   endif else begin
     cd, pbs_dir
     if keyword_set(pbs_batch) then begin
        if keyword_set(nosubmit) then begin
            splog, 'Generated '+pbs_batch_script+' but not submitting to queue'
        endif else begin
           splog, 'Submitting '+pbs_batch_script
           spawn, sub+' '+pbs_batch_script
        endelse
     endif else begin
       for i=0L, pbs_nodes-1 do begin
          if keyword_set(nosubmit) then begin
              splog, 'Generated '+pbs_node_script[i]+' but not submitting to queue'
          endif else begin
             splog, 'Submitting '+pbs_node_script[i]
             spawn, sub+' '+pbs_node_script[i]
          endelse
       endfor
     endelse
   endelse

   return
end
;------------------------------------------------------------------------------
