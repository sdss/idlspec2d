; Stupid routine for reading 1D outputs at Princeton

function read1spec, plate, fiber, fluxerr=fluxerr, wave=wave

   filename = '/u/dss/data/spectro/*' + strtrim(string(plate),2) $
    + '/spSpec-*' + string(plate,format='(i4.4)') $
    + '-' + string(fiber,format='(i3.3)') + '.fit'
   filename = (findfile(filename))[0]

   data = readfits(filename, hdr)
   flux = data[*,0]
   fluxerr = data[*,2]

   coeff0 = sxpar(hdr, 'COEFF0')
   coeff1 = sxpar(hdr, 'COEFF1')
   wave = 10^( coeff0 + coeff1 * findgen(n_elements(flux)) )

   return, flux
end
