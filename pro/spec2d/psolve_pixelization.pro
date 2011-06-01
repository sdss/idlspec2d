;------------------------------------------------------------------------------
function psolve_pixelization, pradius=pradius, rradius=rradius

   nmap = n_elements(pradius)
   nhalf = ceil(max([pradius,rradius]))
   nside = 2 * nhalf + 1
   rmap = shift(dist(nside,nside), nhalf, nhalf)

   ntot = 0L
   psfpix = lonarr(nside,nside,nmap)
   for imap=0, nmap-1 do begin
      psfpix1 = lonarr(nside,nside) - 1

      indx1 = where(rmap LT pradius[imap], n1)
      if (n1 GT 0) then $
       psfpix1[indx1] = lindgen(n1) + ntot
      ntot = ntot + n1

      indx2 = where(rmap GE pradius[imap] AND rmap LT rradius[imap], n2)
      if (n2 GT 0) then begin
         thisval = floor(rmap[indx2] - pradius[imap])
         psfpix1[indx2] = thisval + ntot
         ntot = ntot + max(thisval) + 1
      endif

      psfpix[*,*,imap] = psfpix1
   endfor

   return, psfpix
end
;------------------------------------------------------------------------------
