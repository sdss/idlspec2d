function spectra_sample, samples, ver1Ddir, list = list

  if (keyword_set(list)) then begin
   inputlist = samples 
   readcol, inputlist, samples, format='(a)'
  endif

  nfiles=n_elements(samples)

  for j=0,nfiles-1 do begin
 
     sample=mrdfits(samples[j],1)

     good=where(sample.plate NE 0)
     ngals=n_elements(good)
     plate=sample[good[0]].plate
     mjd=sample[good[0]].mjd

     i=10
     ff = string('/usr/sdss/data05/spectro/'+ver1Ddir+'/' $
           +string(plate,format='(i4.4)')+$
        '/1d/spSpec-'+string(mjd,format='(i5.5)')+'-'$
         +string(plate,format='(i4.4)')+'-'$
         +string(sample[good[i]].fiberid,format='(i3.3)')+'.fit')
     spec = mrdfits(ff,0,hdr)
     npix= (size(spec[*,0]))[1]
     galspec=fltarr(npix,4,640)
     npix1=npix - 1

     for i=0,ngals-1 do begin

         ff = string('/usr/sdss/data05/spectro/'+ver1Ddir+'/' $
           +string(plate,format='(i4.4)')+$
        '/1d/spSpec-'+string(mjd,format='(i5.5)')+'-'$
         +string(plate,format='(i4.4)')+'-'$
         +string(sample[good[i]].fiberid,format='(i3.3)')+'.fit')

         spec = mrdfits(ff,0,hdr)
        
        galspec[0:npix1,0,good[i]]=spec[0:npix1,0]
        galspec[0:npix1,1,good[i]]=spec[0:npix1,2]
        zp=sxpar(hdr,'CRVAL1')
        galspec[0:npix1,2,good[i]]=zp+findgen(npix)*0.0001
        galspec[0:npix1,3,good[i]]=spec[0:npix1,3]

	print,i,good[i]
     endfor

     namesample=string('spectra_'+string(mjd,format='(i5.5)')+'_' $
                    +string(plate, format='(i4.4)')+'.fit')
     mwrfits, galspec, namesample, /create
  
     endfor
     
return, galspec

end