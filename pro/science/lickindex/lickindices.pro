function lickindices, samples, ver1Ddir, contype, test = test, list = list

; contype=0  without continuum
; contype=1  continuum obtained splining only the pseudo continuum regions
; contype=2  continuum obtained splining all the wavelength range
; contype=3  continuum obtained with a straight line

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

     for i=0,ngals-1 do begin

         ff = string('/usr/sdss/data05/spectro/'+ver1Ddir+'/' $
           +string(plate,format='(i4.4)')+$
        '/1d/spSpec-'+string(mjd,format='(i5.5)')+'-'$
         +string(plate,format='(i4.4)')+'-'$
         +string(sample[good[i]].fiberid,format='(i3.3)')+'.fit')

         spec = mrdfits(ff,0,hdr)

         npix= (size(spec[*,0]))[1]
	 wave=sxpar(hdr,'COEFF0')+ findgen(npix)*0.0001

         one_indices, spec[*,0], wave, spec[*,2], sample[good[i]].z_1D, $
                     contype, abslines, abslineserr, sn

         sample[good[i]].lickindex=abslines
         sample[good[i]].lickindexerr=abslineserr

    endfor   

    if (keyword_set(test)) then $
    namesample=string('spec_param_'+string(mjd,format='(i5.5)')+'_' $
                    +string(plate, format='(i4.4)') $
                    +string(test, format='(a)')+'.fit') $
    else $
    namesample=string('spec_param_'+string(mjd,format='(i5.5)')+'_' $
                    +string(plate, format='(i4.4)')+'.fit')  
    mwrfits, sample, namesample, /create

  endfor

return, sample

end