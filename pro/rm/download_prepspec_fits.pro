; download the fits from prepspec

pro download_prepspec_fits, suffix=suffix, rmid=rmid, topurl=topurl, outdir=outdir, subdir=subdir

   if ~keyword_set(topurl) then topurl='http://star-www.st-and.ac.uk/~kdh1/pub/sdss/2014b/'
   if ~keyword_set(outdir) then outdir='/data3/yshen/ftp/sdssrm/collab/prepspec/2014b/'
   if ~keyword_set(subdir) then subdir='/ACBFJ/'  ;'/ACBFJ/','/ACBHYGF/'

   if n_elements(rmid) eq 0 then begin
       ; rmid=indgen(849L)
       ; default is the first 100 objects
       file = getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
       target=mrdfits(file,1)
       target=target[0:848]
       ind=sort(target.zfinal)
       rmid=ind ;[0:99]
   endif


   if ~keyword_set(suffix) then begin
     ;suffix=['avgrms.ps','ACF_fit.ps','ACBHYGF_fit.ps', $ 
     ;        'c5100.dat', 'blr_t.dat','ha_t.dat', 'hb_t.dat', 'he2_t.dat', 'mg2_t.dat', $
     ;        'avg.dat', 'rms.dat']
     suffix=['c5100.dat', 'c3000.dat', 'c2500.dat', 'c1700.dat', 'c1350.dat', $
        'avg_w.dat', 'rms_w.dat','nlr_w.dat', $
        'ha_t.dat','hb_t.dat', 'he2_4686_t.dat','mg2_t.dat', 'c3_t.dat', 'c4_t.dat', 'si4_t.dat', $
        'n5_t.dat', 'lya_t.dat', 'p0_t.dat', $
        'vnlr.dat', 'vblr.dat', 'vvblr.dat', $
        'blr_w.dat', 'ha_w.dat', 'hb_w.dat', 'he2_4686_w.dat', 'mg2_w.dat', $
        'c3_w.dat', 'c4_w.dat', 'si4_w.dat', 'n5_w.dat', 'lya_w.dat', $
        'c5100_stats.dat', 'c3000_stats.dat', 'c2500_stats.dat', 'c1700_stats.dat', 'c1350_stats.dat', $ ; stats
        'ha_t_stats.dat','hb_t_stats.dat', 'he2_4686_t_stats.dat','mg2_t_stats.dat', 'c3_t_stats.dat', 'c4_t_stats.dat', 'si4_t_stats.dat', $
        'n5_t_stats.dat', 'lya_t_stats.dat' ]
   endif
   ntag=n_elements(suffix)

   nobj=n_elements(rmid)
   for i=0L, nobj-1 do begin
      rmtag0=string(rmid[i],format='(i3.3)')
      rmtag = 'rm'+string(rmid[i],format='(i3.3)')
      outdir1 = outdir + rmtag + '/'
      tmp = file_test(outdir1,/dir)
      if tmp eq 0 then spawn, 'mkdir ' + outdir1
      cd, outdir1

      for j=0L, ntag-1 do begin
         file = topurl + rmtag0 + subdir + rmtag + '_' + suffix[j]
         spawn, 'wget ' + file
      endfor
   endfor
   
end
