;+
; NAME:
;   sdss_plate_sort
;
; PURPOSE:
;   Resort the photoPlate files to the fiber numbering of a plugmap
;
; CALLING SEQUENCE:
;   sdss_plate_sort, [ planfile ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name(s) of output plan file; default to looping through
;                all plan files matching 'spPlan2d*.par'
;
; OUTPUT:
;
; COMMENTS:
;   The following files are input:
;     $PHOTOPLATE_DIR/$PLATE/photoMatchPlate-$PLATE.fits
;     $PHOTOPLATE_DIR/$PLATE/photoPlate-$PLATE.fits
;     $PHOTOPLATE_DIR/$PLATE/photoPosPlate-$PLATE.fits
;   The following files are output in the same directories as
;   the plan files:
;     photoMatchPlate-$PLATE-$MJD.fits
;     photoPlate-$PLATE-$MJD.fits
;     photoPosPlate-$PLATE-$MJD.fits
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Apr-2010  Written by David Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro sdss_plate_sort, planfile

   if (NOT keyword_set(planfile)) then planfile = findfile('spPlan2d*.par')

   ;----------
   ; If multiple plan files exist, then call this script recursively
   ; for each such plan file.

   nplan = n_elements(planfile)
   if (nplan GT 1) then begin
      for i=0, nplan-1 do $
       sdss_plate_sort, planfile[i]
      return
   endif

   ;----------
   ; Strip path from plan file name, and change to that directory

   thisplan = fileandpath(planfile[0], path=thispath)
   cd, thispath, current=origdir
   if (NOT keyword_set(thispath)) then cd, origdir

   ;----------
   ; Find the SPEXP structure

   allseq = yanny_readone(thisplan, 'SPEXP', hdr=hdr, /anon)
   if (N_elements(allseq) EQ 0) then begin
      splog, 'ABORT: No SPEXP structures in plan file ' + thisplan
      cd, origdir
      return
   endif

   ;----------
   ; Find keywords from the header

   sortstr = string(allseq.plateid) + ' ' + string(allseq.mjd)
   ilist = uniq(sortstr,uniq(sortstr))
   for i=0, n_elements(ilist)-1 do begin
      platestr = string(allseq[ilist[i]].plateid,format='(i4.4)')
      mjdstr = string(allseq[ilist[i]].mjd,format='(i5.5)')
      plugfile = 'plPlugMapM-'+allseq[ilist[i]].mapname+'.par'
      plugdir = getenv('SPECLOG_DIR')+'/'+mjdstr
      matchfile = (findfile(djs_filepath('photoMatchPlate-'+platestr+'.fits*', $
       root_dir=getenv('PHOTOPLATE_DIR'), subdir=platestr)))[0]
      infile1 = (findfile(djs_filepath('photoPlate-'+platestr+'.fits*', $
       root_dir=getenv('PHOTOPLATE_DIR'), subdir=platestr)))[0]
      infile2 = (findfile(djs_filepath('photoPosPlate-'+platestr+'.fits*', $
       root_dir=getenv('PHOTOPLATE_DIR'), subdir=platestr)))[0]
      outfile0 = 'photoMatchPlate-'+platestr+'-'+mjdstr+'.fits'
      outfile1 = 'photoPlate-'+platestr+'-'+mjdstr+'.fits'
      outfile2 = 'photoPosPlate-'+platestr+'-'+mjdstr+'.fits'

      plugmap = readplugmap(plugfile, plugdir=plugdir)
      if (keyword_set(matchfile)) then $
       matchdat = mrdfits(matchfile, 1, hdr0) $
      else matchdat = 0
      if (keyword_set(infile1)) then $
       objdat1 = mrdfits(infile1, 1, hdr1) $
      else objdat1 = 0
      if (keyword_set(infile2)) then $
       objdat2 = mrdfits(infile2, 1, hdr2) $
      else objdat2 = 0

      qplug_exist = keyword_set(plugmap)
      qobj_exist = keyword_set(matchdat) * keyword_set(objdat1) $
       * keyword_set(objdat2)

      if (qplug_exist EQ 0) then $
       splog, 'WARNING: Missing plPlugMapM files!'
      if (qobj_exist EQ 0) then $
       splog, 'WARNING: Missing photoPlate files!'

      if (qplug_exist AND qobj_exist) then begin

         splog, 'Input photoMatchPlate file = '+matchfile
         splog, 'Input photoPlate file = '+infile1
         splog, 'Input photoPosPlate file = '+infile2

         spherematch, plugmap.ra, plugmap.dec, $
          ((matchdat.match_ra+360.0) MOD 360), matchdat.match_dec, $
          1./3600, i1, i2, d12
         nfiber = n_elements(plugmap)

         if (n_elements(i1) NE nfiber) then $
          message, 'Failure matching all objects!'
         if (total(i1[sort(i1)] EQ lindgen(nfiber)) LT nfiber) then $
          message, 'Double-matching of some objects!'
         if (total(i2[sort(i2)] EQ lindgen(nfiber)) LT nfiber) then $
          message, 'Double-matching of some objects!'
         isort = lonarr(nfiber)
         isort[i1] = lindgen(nfiber)

         splog, 'Writing '+outfile0
         mwrfits, matchdat[i2[isort]], outfile0, hdr0, /create
         splog, 'Writing '+outfile1
         mwrfits, objdat1[i2[isort]], outfile1, hdr1, /create
         splog, 'Writing '+outfile2
         mwrfits, objdat2[i2[isort]], outfile2, hdr2, /create
      endif
   endfor

   ; Change back to original directory
   cd, origdir
   return
end
;------------------------------------------------------------------------------
