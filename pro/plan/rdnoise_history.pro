;+
; NAME:
;   rdnoise_history
;
; PURPOSE:
;   Compute the noise history for the BOSS CCDs
;
; CALLING SEQUENCE:
;   rdnoise_history, [docams=, mjd= ]
;
; INPUTS:
;
; REQUIRED KEYWORDS:
;
; OPTIONAL KEYWORDS:
;   docams     - Camera(s); default to ['b1','b2','r1','r2']
;   mjd        - Search for all files in $BOSS_SPECTRO_DATA/MJD;
;                default to '55???'
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Output files with 12 columns: file name, MJD, TAI, exposure time (seconds),
;   read noise in DN in each of the 4 quadrants, standard deviation of the
;   bias drift in DN (row-to-row) in each of the 4 quadrants.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   djs_avsigclip()
;   djsig()
;   fileandpath()
;   sdssproc
;   splog
;
; REVISION HISTORY:
;   21-Oct-2009  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro rdnoise_history, docams=docams, mjd=mjd1

   if (NOT keyword_set(docams)) then docams = ['b1','b2','r1','r2']
   if (n_elements(docams) GT 1) then begin
      for i=0, n_elements(docams)-1 do $
       rdnoise_history, docams=docams[i], mjd=mjd1
      return
   endif

   if (keyword_set(mjd1)) then mjdstr = strtrim(mjd1,2) $
    else mjdstr = '55???'

   mjdlist = ''
   for i=0, n_elements(mjdstr)-1 do begin
      spawn, '\ls -d '+getenv('BOSS_SPECTRO_DATA') + '/' + mjdstr[i], mjdlist1
      if (keyword_set(mjdlist1)) then $
       mjdlist = keyword_set(mjdlist) ? [mjdlist,mjdlist1] : mjdlist1
   endfor
   if (NOT keyword_set(mjdlist)) then begin
      print, 'No files found'
      return
   endif
   nmjd = n_elements(mjdlist)

   for imjd=0, nmjd-1 do begin
      thisname = filepath('sdR-'+docams+'-????????.fit*', $
       root_dir=mjdlist[imjd])
      thisname = findfile(thisname, count=ct)
      if (ct GT 0) then begin
         if (NOT keyword_set(framename)) then framename = thisname $
          else framename = [framename, thisname]
      endif
   endfor

   camcolor = strmid(docams,0,1)
   case camcolor of
   'b': begin
      xbias0 = [  0, 4284,   0,  4284]
      xbias1 = [ 67, 4340,  67,  4340]
      y0 = [  56,   56, 2112,  2112]
      y1 = [2111, 2111, 4167,  4167]
      end
   'r': begin
      xbias0 = [  0, 4241,  0,  4241]
      xbias1 = [110, 4351,110,  4351]
      y0 = [ 100,  100, 2112,  2112]
      y1 = [2111, 2111, 4123,  4123]
      end
   endcase

   splog, filename='rdnoise-'+docams+'.log', /noname, $
    '# File MJD TAI EXPTIME Rdnoise0 Rdnoise1 Rdnoise3 Drift0 Drift1 Drift2 Drift3'
   splog, /noname, $
    '# All noise values in DN'

   nfile = n_elements(framename)
   mjd = lonarr(nfile)
   tai = dblarr(nfile)
   exptime = fltarr(nfile)
   rdnoise = fltarr(nfile,4)
   rddrift = fltarr(nfile,4)
   for ifile=0, nfile-1 do begin
      print, 'Reading file ', ifile, ' of ', nfile, ': ' + framename[ifile]
      sdssproc, framename[ifile], hdr=hdr, /silent
      exptime[ifile] = sxpar(hdr, 'EXPTIME')
      mjd[ifile] = sxpar(hdr, 'MJD')
      tai[ifile] = sxpar(hdr, 'TAI-BEG')
      img = mrdfits(framename[ifile],/fscale)
      if (keyword_set(img)) then begin
         if (mjd[ifile] LT 55113) then img = rotate(img,2)
         for j=0,3 do begin
            biasimg = img[xbias0[j]:xbias1[j],y0[j]:y1[j]]
            rdnoise[ifile,j] = djsig(biasimg) * 1.015
            bvec = djs_avsigclip(biasimg, 1)
            rddrift[ifile,j] = stddev(bvec)
         endfor
      endif
      splog, fileandpath(framename[ifile]), $
       mjd[ifile], tai[ifile], (exptime[ifile]>0)<99999, $
       rdnoise[ifile,*], rddrift[ifile,*], $
       format='(a,i6,f12.0,f8.1,8f10.3)', /noname
   endfor

   splog, /close

   return
end
;------------------------------------------------------------------------------
