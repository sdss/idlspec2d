; List of number of galaxies and quasars per plate, the number of redshifts
; that differ from prior reductions, and the fraction of various 2D flags set.

pro redmonsterlist, newdir=newdir, olddir=olddir

;setenv,'SPECTRO_DATA=/scr0/2d_v4_8'
   if (NOT keyword_set(olddir)) then olddir = '/u/dss/spectro_v4_7'
   if (NOT keyword_set(newdir)) then newdir = getenv('SPECTRO_DATA')

   splog, filename='monsterlist.log'

   platelist, plist=plist

   splog, 'PLATE MJD   %Gal   %QSO   Ndiff    Whoppr Reject Scattr X-Talk BrtSky BadFlx BadChi Monstr', /noname

   splog, '----- ----- ------ ------ ------   ------ ------ ------ ------ ------ ------ ------ ------', /noname

   for i=0, n_elements(plist)-1 do begin
      if (strtrim(plist[i].status2d) EQ 'Done') then begin

         setenv, 'SPECTRO_DATA='+newdir
         readspec, plist[i].plate, mjd=plist[i].mjd, $
          plug=plug, ormask=ormask, zans=zans, /silent

         setenv, 'SPECTRO_DATA='+olddir
         readspec, plist[i].plate, mjd=plist[i].mjd, zans=oldz, /silent

         if (keyword_set(zans) AND keyword_set(oldz)) then begin
            zdiff = 3.e5 * abs(zans.z - oldz.z)
            junk = where(zdiff GT 500. AND zans.zwarning EQ 0 $
             and oldz.zwarning EQ 0, ndiff)
         endif else begin
            ndiff = 0
         endelse

         igal = where((plug.primtarget AND $
          (2L^6 + 2L^7 + 2L^8 + 2L^5 + 2L^26)) NE 0, ngal)
         iqso = where((plug.primtarget AND $
          (2L^0 + 2L^1 + 2L^2 + 2L^3 + 2L^4 + 2L^25)) NE 0, nqso)

         if (ngal GT 0 AND keyword_set(zans)) then begin
            junk = where(zans[igal].zwarning EQ 0 $
             AND strtrim(zans[igal].class) EQ 'GALAXY', nggal)
            fgal = 100 * float(nggal) / ngal
         endif else begin
            fgal = 0.
         endelse
         if (nqso GT 0 AND keyword_set(zans)) then begin
            junk = where(zans[iqso].zwarning EQ 0 $
             AND strtrim(zans[iqso].class) EQ 'QSO', ngqso)
            fqso = 100 * float(ngqso) / nqso
         endif else begin
            fqso = 0.
         endelse

         junk = where((ormask AND (2L^8+2L^9)) NE 0, nwhopper)
         junk = where((ormask AND (2L^18+2L^19)) NE 0, nreject)
         junk = where((ormask AND 2L^20) NE 0, nscatter)
         junk = where((ormask AND 2L^21) NE 0, nxtalk)
         junk = where((ormask AND 2L^23) NE 0, nbsky)
         junk = where((ormask AND 2L^26) NE 0, nbadflux)
         junk = where((ormask AND 2L^27) NE 0, nbadchi)
         junk = where((ormask AND 2L^28) NE 0, nmonster)

         npix = n_elements(ormask)
         fwhopper = 100 * float(nwhopper) / npix
         freject = 100 * float(nreject) / npix
         fscatter = 100 * float(nscatter) / npix
         fxtalk = 100 * float(nxtalk) / npix
         fbsky = 100 * float(nbsky) / npix
         fbadflux = 100 * float(nbadflux) / npix
         fbadchi = 100 * float(nbadchi) / npix
         fmonster = 100 * float(nmonster) / npix

         splog, string(plist[i].plate, plist[i].mjd, fgal, fqso, ndiff, $
          fwhopper, freject, fscatter, fxtalk, fbsky, fbadflux, $
          fbadchi, fmonster, $
          format='(i5,i6,2f7.1,i7,"  ",8f7.2)'), /noname

      endif
   endfor

   splog, /close
end
