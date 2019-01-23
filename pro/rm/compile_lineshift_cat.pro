; Combine the fitting results for the line shifts study

pro compile_lineshift_cat

   ; this is the results from the global fit
   ; file='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/qso_prop-0000-56837_lineshift.fits'
   file = '/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/evenmore_lines/civ_3gauss/qso_prop-0000-56837.fits'
   result=mrdfits(file,1)

   ; now add tags for NeV3426, OII3728 and CaII3934
   file='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/local_fit/qso_prop-0000-56837_lineshift.fits'
   result1=mrdfits(file,1)
   alltag=tag_names(result1)
   ind=where(strmatch(alltag, 'CAII3934*') eq 1 or $
             strmatch(alltag, 'OII3728*') eq 1  or $
             strmatch(alltag, 'NEV3426*') eq 1, comp=indd)
   tag_rej=alltag[indd]
   remove_tags, result1, tag_rej, result2

   result=struct_addtags(result, result2)

   ; output the fits catalog
   ; outfile='/data3/yshen/work/lineshifts/lineshift.fits'
   outfile = '/data3/yshen/work/lineshifts/lineshift_civ3gauss.fits'
   mwrfits, result, outfile, /create

end

; the routine below will compile the results for the results on different Sn
pro compile_lineshift_cat_sn

err_scal_arr=[2.,4.,6.,8.,10.]
nerr=n_elements(err_scal_arr)
topdir='/data3/yshen/spectro/bossredux/v5_7_1/0000/qsofit/wh_skysub/'

for i=0L, nerr - 1 do begin

   ; flesh out the output structure
   junk=temporary(result)
   
   ; this is the results from the global fit
   file=topdir+'evenmore_lines/err_scal_' $
     + string(err_scal_arr[i],format='(f0.1)') + '/' + 'qso_prop-0000-56837.fits'
   if file_test(file) eq 1 then $
     result=mrdfits(file,1)

   ; now add tags for NeV3426, OII3728 and CaII3934
   file=topdir+'local_fit/err_scal_' $
     + string(err_scal_arr[i],format='(f0.1)') + '/' + 'qso_prop-0000-56837.fits'
   result1=mrdfits(file,1)
   alltag=tag_names(result1)
   ind=where(strmatch(alltag, 'CAII3934*') eq 1 or $
             strmatch(alltag, 'OII3728*') eq 1  or $
             strmatch(alltag, 'NEV3426*') eq 1, comp=indd)
   tag_rej=alltag[indd]
   remove_tags, result1, tag_rej, result2

   if n_elements(result) ne 0 then result=struct_addtags(result, result2) else $
     result=result2

   ; output the fits catalog
   outfile='/data3/yshen/work/lineshifts/lineshift_err_scal_'+ string(err_scal_arr[i],format='(f0.1)') + '.fits'
   mwrfits, result, outfile, /create

endfor

end
