
function fitvacset, xpeak, lambda, wset, xset, ncoeff=ncoeff

      if (NOT keyword_set(ncoeff)) then ncoeff=5

      ;------------------
      ; First convert lambda, and skywaves to log10 vacuum

      splog, 'Converting wavelengths to vacuum'
      vaclambda = lambda
      airtovac, vaclambda
      vacloglam = alog10(vaclambda)

      splog, 'Tweaking to sky lines'

      if (size(xset,/tname) EQ 'UNDEFINED') then begin
         splog, 'WARNING: Sky lines are too noisy! No shifting!'
         xshift = xpeak*0.0
      endif else begin
	 traceset2xy, xset, transpose(xpeak), xshift 
         xshift = transpose(xshift)

         ; Move this QA plot elsewhere ???
         plot, xpeak, xshift, ps=3, xtitle='Arc line position', /ynozero, $
           ytitle = 'Offset to Sky Line [pix]', $
           title = 'Offset to Sky Lines From Wavelength-Solution'

      endelse

      vacset = wset
      nfiber = (size(xpeak))[1]

      xy2traceset, transpose(double(xpeak+xshift)), $
                   vacloglam # (dblarr(nfiber)+1), $
                   vacset, ncoeff=ncoeff, xmin=wset.xmin, xmax=wset.xmax

      return, vacset
end



