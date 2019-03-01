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
;   09-Apr-2014  Fix upstream bug where failed position matches could
;                still have thing_id>=0.  SJB
;                Also replace failed match ra,dec with plug_ra,plug_dec
;-
;------------------------------------------------------------------------------
pro sdss_field_sort, planfile

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

   sortstr = string(allseq.fieldid) + ' ' + string(allseq.mjd)
   ilist = uniq(sortstr,uniq(sortstr))
   for i=0, n_elements(ilist)-1 do begin
      fieldstr = field_to_string(allseq[ilist[i]].fieldid) ;- JEB plate number problem OK
      mjdstr = string(allseq[ilist[i]].mjd,format='(i5.5)')
      calobjobssfile = 'obsSummary-'+allseq[ilist[i]].mapname+'.par'
      plugdir = getenv('SDSSCORE')+'/'+mjdstr
      matchfile = (findfile(djs_filepath('photoMatchField-'+fieldstr+'*.fits*', $
       root_dir=getenv('SDSSCORE')+'/photofield/', subdir=fieldstr)))[0]
      infile1 = (findfile(djs_filepath('photoField-'+fieldstr+'*.fits*', $
       root_dir=getenv('SDSSCORE')+'/photofield/', subdir=fieldstr)))[0]
      infile2 = (findfile(djs_filepath('photoPosField-'+fieldstr+'*.fits*', $
       root_dir=getenv('SDSSCORE')+'/photofield/', subdir=fieldstr)))[0]
      outfile0 = 'photoMatchField-'+fieldstr+'-'+mjdstr+'.fits'
      outfile1 = 'photoField-'+fieldstr+'-'+mjdstr+'.fits'
      outfile2 = 'photoPosField-'+fieldstr+'-'+mjdstr+'.fits'
            
      plugmap = readobssummary(calobjobssfile, plugdir=plugdir)
      
      if (keyword_set(matchfile)) then $
       matchdat = mrdfits(matchfile, 1, hdr0) $
      else matchdat = 0
      if (keyword_set(infile1)) then $
       objdat1 = mrdfits(infile1, 1, hdr1) $
      else objdat1 = 0
      if (keyword_set(infile2)) then $
       objdat2 = mrdfits(infile2, 1, hdr2) $
      else objdat2 = 0
      ;print, 'photoPosField-'+fieldstr+'*.fits*'
      ;print, getenv('SDSSCORE')+'/photofield/'+fieldstr
      ;----------
      ; Fix bug where position match failures could still have thing_id >= 0

      ibad = where((objdat2.ra EQ 0.0) AND (objdat2.dec EQ 0.0) $
            AND (objdat2.thing_id GE 0), nbad)
      if (nbad GT 0) then begin
          splog, "Correcting " + strtrim(string(nbad),2) $
              + " failed matches that still had thing_id>=0"
          objdat2[ibad].thing_id = -1
      endif

      qplug_exist = keyword_set(plugmap)
      qobj_exist = keyword_set(matchdat) * keyword_set(objdat1) $
       * keyword_set(objdat2)

      if (qplug_exist EQ 0) then $
       splog, 'WARNING: Missing obsSummary files!'
      if (qobj_exist EQ 0) then $
       splog, 'WARNING: Missing photoField files!'

      if (qplug_exist AND qobj_exist) then begin

         splog, 'Input photoMatchField file = '+matchfile
         splog, 'Input photoField file = '+infile1
         splog, 'Input photoPosField file = '+infile2

         ;----------
         ; spherematch plugmap and match file

         spherematch, plugmap.ra, plugmap.dec, $
          ((matchdat.match_ra+360.0) MOD 360), matchdat.match_dec, $
          1./3600, i1, i2, d12
         nfiber = n_elements(plugmap)

         if (n_elements(i1) NE nfiber) then $
          message, 'ERROR: Failure matching all objects!'
         if (total(i1[sort(i1)] EQ lindgen(nfiber)) LT nfiber) then $
          message, 'ERROR: Double-matching of some objects!'
         if (total(i2[sort(i2)] EQ lindgen(nfiber)) LT nfiber) then $
          message, 'ERROR: Double-matching of some objects!'
         isort = lonarr(nfiber)
         isort[i1] = lindgen(nfiber)

         ;----------
         ; Reorder match arrays to plugmap fiber order

         matchdat = matchdat[i2[isort]]
         objdat1 = objdat1[i2[isort]]
         objdat2 = objdat2[i2[isort]]

         ;----------
         ; Replace failed matches with plug_ra,plug_dec
         ; Photometric information is already 0, and thing_id=-1

         nomatch = where(objdat2.thing_id LT 0, nmiss)
         if (nmiss GT 0) then begin
            splog, "Updating " + strtrim(string(nmiss),2) $
                + " failed matches to use plug_ra,plug_dec"
            objdat2[nomatch].ra  = plugmap[nomatch].ra
            objdat2[nomatch].dec = plugmap[nomatch].dec
         endif

         ;----------
         ; Write output files

         splog, 'Writing '+outfile0
         mwrfits, matchdat, outfile0, hdr0, /create
         splog, 'Writing '+outfile1
         mwrfits, objdat1, outfile1, hdr1, /create
         splog, 'Writing '+outfile2
         mwrfits, objdat2, outfile2, hdr2, /create
      endif
   endfor

   ; Change back to original directory
   cd, origdir
   return
end
;------------------------------------------------------------------------------
