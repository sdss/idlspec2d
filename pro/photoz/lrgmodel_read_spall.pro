;------------------------------------------------------------------------------
function lrgmodel_read_spall, public=public, regenerate=regenerate

   if (size(public,/tname) EQ 'STRING') then pstring = '-'+public[0] $
    else if (keyword_set(public)) then pstring = '-public' $
    else pstring = ''
   trimfile = 'spAll-photoz' + pstring + '.fits'

   foo = findfile(trimfile, count=ct)
   if (ct EQ 0 OR keyword_set(regenerate)) then begin

      ; Read the spAll file
      spfile = filepath('spAll'+pstring+'.fits', $
       root_dir=getenv('BOSS_SPECTRO_REDUX'))
      foo = findfile(spfile, count=ct)
      if (ct EQ 0) then message, 'Failed to find file ' + spfile

      columns = ['PROGNAME', 'PLATEQUALITY', 'PLATE', 'FIBERID', 'MJD', $
       'RUN', 'RERUN', 'CAMCOL', 'FIELD', 'ID', 'COLC', $
       'CLASS', 'SPECPRIMARY', 'PRIMTARGET', $
       'Z', 'Z_ERR', 'ZWARNING', $
       'PETROFLUX', 'PSFFLUX', 'DEVFLUX', 'PETROR50', 'R_DEV', $
       'MODELFLUX', 'MODELFLUX_IVAR', 'EXTINCTION']
      spall = hogg_mrdfits(spfile, 1, columns=columns, $
       nrowchunk=10000L)
;       nrowchunk=10000L, range=[0,25000L]) ; ??? Test

      ; Trim to LRGs
      itrim = where( (strmatch(spall.platequality,'good*') $
        OR strmatch(spall.platequality,'marginal*')) $
       AND strmatch(spall.progname,'main*') $
       AND strmatch(spall.class,'GALAXY*') $
;       AND total(spall.modelflux_ivar GT 0,1) EQ 5 $ ; photom in all bands
       AND total(spall.modelflux_ivar[1:3] GT 0,1) EQ 3 $ ; photom in 3 bands
       AND (spall.primtarget AND (2L^5+2L^26)) NE 0 $ ; Targetted as LRG
       AND spall.zwarning EQ 0 $
       AND spall.specprimary EQ 1, ntrim) ; Best spectroscopic observations
      spall = spall[itrim]

      mwrfits, spall, trimfile, /create
   endif

   if (NOT keyword_set(spall)) then spall = mrdfits(trimfile, 1)
   return, spall
end
;------------------------------------------------------------------------------
