; Script to test running DESIGN_PLATE
;------------------------------------------------------------------------------
pro design_random

   blankobj = create_struct( $
    'holetype', '', $
    'objtype' , '', $
    'mag'     , fltarr(5)+25, $
    'ra'      , 0.d, $
    'dec'     , 0.d )

   nguide = 50
   nstd = 30
   nsky = 1000
   nobj = 2000
   ntot = nguide + nstd + nsky + nobj

   racen = 50.0
   deccen = 30.0
   objs = replicate(blankobj, ntot)
   objs.ra = randomu(12345L, ntot)*4 + racen - 2
   objs.dec = randomu(67890L, ntot)*4 + deccen - 2

   objs[0:nguide-1].holetype = 'GUIDE'

   objs[nguide:nguide+nstd-1].holetype = 'OBJECT'
   objs[nguide:nguide+nstd-1].objtype = 'SPECTROPHOTO_STD'

   objs[nguide+nstd:nguide+nstd+nsky-1].holetype = 'OBJECT'
   objs[nguide+nstd:nguide+nstd+nsky-1].objtype = 'SKY'

   objs[nguide+nstd+nsky:nguide+nstd+nsky+nobj-1].holetype = 'OBJECT'
   objs[nguide+nstd+nsky:nguide+nstd+nsky+nobj-1].objtype = 'SERENDIPITY_MANUAL'

   design_plate, objs, racen=racen, deccen=deccen, tilenum=666, nminsky=200

;   plug = yanny_readone('plPlugMapT-0666.par')
;   splot, plug.ra, plug.dec, psym=4
;   isky = where(strmatch(plug.holetype,'COHERENT_SKY*'))
;   soplot, plug[isky].ra, plug[isky].dec, $
;    psym=4, color='green', symsize=1.5

   return
end
;------------------------------------------------------------------------------
