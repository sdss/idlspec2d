;------------------------------------------------------------------------------
function lrgmodel_append_twodf, spall, prefix, nadd=nadd

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
   itrim = where(dat1.z_final GT 0.01 AND dat2.modelflux[2] GT 0, ntrim)

   moredat = replicate(blankobj, ntrim)
   copy_struct, dat2[itrim], moredat
   moredat.z = dat1[itrim].z_final

   nadd = ntrim
   return, [spall, moredat]
end
;------------------------------------------------------------------------------

