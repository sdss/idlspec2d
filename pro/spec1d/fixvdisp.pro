; * Fix the MJD in the output files spZbest, spZall to be consistent
;   with the file name.
; * Compute VDISP,VDISP_ERR for galaxies and add to structures.
;------------------------------------------------------------------------------
pro fixvdisp

   ;---------------------------------------------------------------------------
   ; Fix spZbest files

   filename = findfile('0*/spZbest*.fits', count=nfile)

   for ifile=0, nfile-1 do begin
      print, 'Reading ' + filename[ifile]

      zans = mrdfits(filename[ifile], 1, zhdr)
      spawn, 'mv -f ' + filename[ifile] + ' ' + filename[ifile]+'.OLD'

      mjdname = long( strmid( fileandpath(filename[ifile]), 13, 5) )

      ; Compute the velocity dispersions only if there are galaxies
      ; in this file, and VDISP doesn't yet exist in the structure.
      igal = where(strtrim(zans.class,2) EQ 'GALAXY')
      qdo = (where(tag_names(zans) EQ 'VDISP'))[0] EQ -1
      if (qdo EQ 1 AND igal[0] NE -1) then begin
         thisplatefile = repstr(filename[ifile], 'spZbest', 'spPlate')
         objflux = mrdfits(thisplatefile, 0, objhdr)
         objivar = mrdfits(thisplatefile, 1)

         vdispfit, objflux[*,igal], objivar[*,igal], hdr=objhdr, $
          zobj=zans[igal].z, sigma=sigma, sigerr=sigerr

         newans = struct_addtags(zans, $
          replicate(create_struct('VDISP', 0.0, 'VDISP_ERR', 0.0), 640))
         newans[igal].vdisp = sigma
         newans[igal].vdisp_err = sigerr

         sxaddpar, zhdr, 'MJD', mjdname ; Force the MJD to be correct
         mwrfits, 0, filename[ifile], zhdr, /create
         mwrfits, newans, filename[ifile]
      endif
   endfor
end
;------------------------------------------------------------------------------
