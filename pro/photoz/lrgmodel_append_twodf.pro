; Input an existing SPALL structure, and append data from 2dF files.
;------------------------------------------------------------------------------
function lrgmodel_append_twodf, spall, prefix, nadd=nadd, colorcut=colorcut

   blankobj = spall[0]
   struct_assign, {junk:0}, blankobj

   tdf_dir = getenv('SDSS_2DF_DIR')
   file1 = filepath('catalogue'+prefix+'.fits', root_dir=tdf_dir, subdir='data')
   file2 = filepath('calibObj-'+prefix+'.fits', root_dir=tdf_dir, subdir='data')
   dat1 = mrdfits(file1, 1)
   dat2 = mrdfits(file2, 1)
   ndat = n_elements(dat1)
   if (n_elements(dat2) NE ndat OR NOT keyword_set(dat1)) then $
    message, 'Problem reading the 2dF data files'

   ; Trim to only good redshifts where we have matched photometry
   qgood = dat1.z_final GT 0.01 AND dat2.modelflux[2] GT 0

   ; Make Nikhil's color cuts if so requested
   if (keyword_set(colorcut)) then begin
      grcolor = - 2.5*alog10((dat2.modelflux[1]>0.1) / (dat2.modelflux[2]>0.1))
      ricolor = - 2.5*alog10((dat2.modelflux[2]>0.1) / (dat2.modelflux[3]>0.1))
      ourcolor = ricolor - grcolor/8.
      qgood = qgood AND ourcolor GT 0.55 AND grcolor GT 1.4
   endif

   itrim = where(qgood, ntrim)
   moredat = replicate(blankobj, ntrim)
   copy_struct, dat2[itrim], moredat
   moredat.z = dat1[itrim].z_final

   ; Nikhil's calibObj files omit RUN,RERUN,CAMCOL,FIELD, which must
   ; be found in the 2dF file.
   moredat.run = dat1[itrim].run
   moredat.rerun = dat1[itrim].rerun
   moredat.camcol = dat1[itrim].camcol
   moredat.field = dat1[itrim].field

   nadd = ntrim
   return, [spall, moredat]
end
;------------------------------------------------------------------------------

