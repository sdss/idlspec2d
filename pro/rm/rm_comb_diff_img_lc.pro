; combine the light curves from Karen's difference imaging work

pro rm_comb_diff_img_lc, specfile=specfile, outdir=outdir

  file=getenv('IDLRM_DIR')+'/etc/target_fibermap.fits'
  target=mrdfits(file,1)
  epoch_info=mrdfits(file,2)

  ; diff LC dir
  topdir='/data3/yshen/ftp/sdssrm/collab/photo_lc/img_diff_lc/'

  topdir_std_bok='/data3/yshen/ftp/sdssrm/collab/photo_lc/std/bok/'

  ; output path for the combined LCs
  if ~keyword_set(outdir) then outdir=topdir

  ; read in synthetic flux from spectra
  if ~keyword_set(specfile) then $
     specfile='/data3/yshen/ftp/sdssrm/collab/photo_lc/spec_syn_prepspec_nov11_2015.fits'
  specflux=mrdfits(specfile,1)

  ;rmid=indgen(849)
  ;nnn=n_elements(rmid)
  ;nep=200L ; reserved for 200 photo epochs
  ;result={RMID:0L, zfinal:0.D, objc_type:0L, MJD:dblarr(nep), $
  ;        g_mag_cfht:dblarr(nep), g_mag_cfht_err:dblarr(nep), $
  ;        g_fnu_cfht:dblarr(nep), g_fnu_cfht_err:dblarr(nep), $
  ;        i_mag_cfht:dblarr(nep), i_mag_cfht_err:dblarr(nep), $
  ;        i_fnu_cfht:dblarr(nep), i_fnu_cfht_err:dblarr(nep), $
  ;        g_mag_bok:dblarr(nep), g_mag_bok_err:dblarr(nep), $
  ;        g_fnu_bok:dblarr(nep), g_fnu_bok_err:dblarr(nep), $
  ;        i_mag_bok:dblarr(nep), i_mag_bok_err:dblarr(nep), $
  ;        i_fnu_bok:dblarr(nep), i_fnu_bok_err:dblarr(nep)}
  ;tags=tag_names(result)
  ;result=replicate(result,nnn)
  ;result.zfinal=target[rmid].zfinal
  ;result.objc_type=target[rmid].objc_type

  ; get the LCs
  band=['g', 'i']
  for iband=0,1 do begin
    ; for CFHT data
    fld=['a','b','c','d','e','f','g','h','i']
    nfld=n_elements(fld)
    ; loop over all fields
    for idith=1,2 do begin
      for ifld=0,nfld-1 do begin
        fldtag=fld[ifld]+string(idith,format='(i0)')+band[iband]
        datdir=topdir + fldtag +'/'
        orig0='cfht_'+fld[ifld]+string(idith,format='(i0)')
        ; get ccd#
        ccdfile=datdir+fldtag+'.dat'
        readcol,ccdfile,format='a,x,x,x,x,x,x,a',/silent, rmid_ccd, ccd_no
        rmid_ccd=long(strmid(rmid_ccd,2))

        files=file_search(datdir+'rm*.dat')
        len=strlen(datdir)
        rmid=long(strmid(strmid(files,len),2, 3))
        nobj=n_elements(rmid) 
        for jj=0, nobj-1 do begin
          readcol, files[jj], format='d,d,d',mjd,flux,flerr,/silent
          ; find the ccd number for each rmid
          ind_ccd=where(rmid_ccd eq rmid[jj])
          if ind_ccd[0] ne -1 then begin
             orig=orig0 + '_' + ccd_no[ind_ccd[0]]
          endif else splog, 'obj not found in dat file: rm', rmid[jj],fldtag+'.dat'
          ntmp=n_elements(mjd)
          if n_elements(mjd_arr) eq 0 then mjd_arr=mjd else mjd_arr=[mjd_arr,mjd]
          if n_elements(flux_arr) eq 0 then flux_arr=flux else flux_arr=[flux_arr,flux]
          if n_elements(flerr_arr) eq 0 then flerr_arr=flerr else flerr_arr=[flerr_arr,flerr]
          if n_elements(rmid_arr) eq 0 then rmid_arr=replicate(rmid[jj],ntmp) else $
             rmid_arr=[rmid_arr, replicate(rmid[jj],ntmp)]
          if n_elements(band_arr) eq 0 then band_arr=replicate(band[iband],ntmp) else $
             band_arr=[band_arr, replicate(band[iband],ntmp)]
          if n_elements(orig_arr) eq 0 then orig_arr=replicate(orig,ntmp) else $
             orig_arr=[orig_arr, replicate(orig,ntmp)]
        endfor
        print, 'band,fld,pos finished:',band[iband],fld[ifld],idith
      endfor
    endfor

    ; for bok data
    fld=string(indgen(18),format='(i2.2)') & nfld=18
    for ipos=1,4 do begin
      for ifld=0,nfld-1 do begin
         datdir=topdir+'bokdata_'+band[iband]+'/'+strupcase(band[iband])+fld[ifld] $
           +'_ccd'+string(ipos,format='(i0)')+'/'
         orig='bok_'+fld[ifld]+'_ccd'+string(ipos,format='(i0)')
         files=file_search(datdir+'rm*.dat')
        len=strlen(datdir)
        rmid=long(strmid(strmid(files,len),2, 3))
        nobj=n_elements(rmid) 
        for jj=0, nobj-1 do begin
          readcol, files[jj], format='d,d,d',mjd,flux,flerr,/silent
          ntmp=n_elements(mjd)
          mjd_arr=[mjd_arr,mjd]
          flux_arr=[flux_arr,flux]
          flerr_arr=[flerr_arr,flerr]
          rmid_arr=[rmid_arr, replicate(rmid[jj],ntmp)]
          band_arr=[band_arr, replicate(band[iband],ntmp)]
          orig_arr=[orig_arr, replicate(orig,ntmp)]
        endfor
        ; now do it for the 70 std stars -- note they are in a different dir
        datdir=topdir_std_bok+strupcase(band[iband])+fld[ifld] $
           +'stnd_ccd'+string(ipos,format='(i0)')+'/'
        files=file_search(datdir+'rm*.dat',count=nfound)
        if nfound gt 0 then begin ; only if a std is found
          len=strlen(datdir)
          rmid=long(strmid(strmid(files,len),2, 3))
          nobj=n_elements(rmid)
          for jj=0, nobj-1 do begin
            readcol, files[jj], format='d,d,d',mjd,flux,flerr,/silent
            ntmp=n_elements(mjd)
            mjd_arr=[mjd_arr,mjd]
            flux_arr=[flux_arr,flux]
            flerr_arr=[flerr_arr,flerr]
            rmid_arr=[rmid_arr, replicate(rmid[jj],ntmp)]
            band_arr=[band_arr, replicate(band[iband],ntmp)]
            orig_arr=[orig_arr, replicate(orig,ntmp)]
          endfor
        endif
        print, 'band,fld,pos finished:',band[iband],fld[ifld],ipos
      endfor
    endfor
  endfor
  ; flip the sign of the flux for the original diff img LC
  flux_arr = -flux_arr

  ; for spec-synthetic
  nep_spec=n_elements(specflux[0].mjd)
  for iobj=0, max(specflux.rmid) do begin
     ; g-band, prepspec flux
     rmid_arr=[rmid_arr, replicate(specflux[iobj].rmid, nep_spec) ]
     band_arr=[band_arr, replicate('g', nep_spec) ]
     orig_arr=[orig_arr, replicate('spec', nep_spec) ]
     mjd_arr=[mjd_arr, specflux[iobj].mjd]
     flux_arr=[flux_arr, reform((specflux[iobj].gri_fnu)[0,*]) ]
     flerr_arr=[flerr_arr, reform((specflux[iobj].gri_fnu_err)[0,*]) ]
     ; g-band, original flux
     rmid_arr=[rmid_arr, replicate(specflux[iobj].rmid, nep_spec) ]
     band_arr=[band_arr, replicate('g', nep_spec) ]
     orig_arr=[orig_arr, replicate('spec_orig', nep_spec) ]
     mjd_arr=[mjd_arr, specflux[iobj].mjd]
     flux_arr=[flux_arr, reform((specflux[iobj].gri_fnu_ori)[0,*]) ]
     flerr_arr=[flerr_arr, reform((specflux[iobj].gri_fnu_ori_err)[0,*]) ]
     ; i-band, prepspec flux
     rmid_arr=[rmid_arr, replicate(specflux[iobj].rmid, nep_spec) ]
     band_arr=[band_arr, replicate('i', nep_spec) ]
     orig_arr=[orig_arr, replicate('spec', nep_spec) ]
     mjd_arr=[mjd_arr, specflux[iobj].mjd]
     flux_arr=[flux_arr, reform((specflux[iobj].gri_fnu)[2,*]) ]
     flerr_arr=[flerr_arr, reform((specflux[iobj].gri_fnu_err)[2,*]) ]
     ; i-band, original flux
     rmid_arr=[rmid_arr, replicate(specflux[iobj].rmid, nep_spec) ]
     band_arr=[band_arr, replicate('i', nep_spec) ]
     orig_arr=[orig_arr, replicate('spec_orig', nep_spec) ]
     mjd_arr=[mjd_arr, specflux[iobj].mjd]
     flux_arr=[flux_arr, reform((specflux[iobj].gri_fnu_ori)[2,*]) ]
     flerr_arr=[flerr_arr, reform((specflux[iobj].gri_fnu_ori_err)[2,*]) ]
     splog, 'Spec RMID finished: ', iobj
  endfor 

  ntot=n_elements(rmid_arr)
  alldata={rmid:0L, mjd:0.D, flux:0.D, flerr:0.D, band:'', orig:''}
  alldata=replicate(alldata, ntot)
  alldata.rmid=rmid_arr & alldata.mjd=mjd_arr & alldata.flux=flux_arr
  alldata.flerr=flerr_arr & alldata.band=band_arr & alldata.orig=orig_arr
  outfile1=outdir + 'alldata.fits'
  mwrfits, alldata, outfile1, /create

  ; now output individual LC for each object
  nnn=max(rmid_arr)+1
  for i=0L, nnn-1 do begin
     ind=where(rmid_arr eq i, nfound)
     if nfound gt 0 then begin
        outfile=outdir+'by_rmid/rm'+string(i,format='(i3.3)')+'.fits'
        output=alldata[ind]
        indd=sort(output.mjd)
        output=output[indd]
        mwrfits, output, outfile, /create,/silent
     endif
  endfor

end


; untar all the img files
pro untar_img_file, topdir=topdir, tag=tag

band=['i','g']
for iband=0,1 do begin  
   if ~keyword_set(topdir) then $
       dir='/data3/yshen/ftp/sdssrm/collab/photo_lc/img_diff_lc/bokdata_' + band[iband] + '/' $
    else dir=topdir
   cd, dir
   for ifld=0,17 do begin
      for iccd=1,4 do begin
         if keyword_set(tag) then name=strupcase(band[iband])+string(ifld,format='(i2.2)')+tag+'_ccd'+string(iccd,format='(i0)') $
           else name=strupcase(band[iband])+string(ifld,format='(i2.2)')+'_ccd'+string(iccd,format='(i0)')
         spawn, 'mkdir ' + name
         cd, name
         spawn, 'mv ../'+name+'.tar ./'
         spawn, 'tar xfv ' + name+'.tar' 
         cd, dir
      endfor
   endfor
endfor
end
