;+
; NAME:
;   boss_galaxy_select.pro
; PURPOSE:
;   Version 0.00001 of the galaxy selection algorithm for BOSS
; CALLING SEQUENCE:
;   ilist=boss_galaxy_select(obj)
; INPUTS:
;   obj   : A calibObj structure
; OPTIONAL INPUTS:
;   
; KEYWORDS:
;   
; OUTPUTS:
;   ilist : A mask array
;     ilist AND 2  : Cut I
;     ilist AND 4  : Cut II
;     ilist AND 8  : Cut III
; OPTIONAL OUTPUTS:
;   
; RESTRICTIONS:
;   
; EXAMPLES:
;   
; COMMENTS:
;   
; REVISION HISTORY:
;  Written by Nikhil Padmanabhan, LBL
;----------------------------------------------------------------------
function lrg_dperp, modelmag

   grcolor = modelmag[1,*]-modelmag[2,*]
   ricolor = modelmag[2,*]-modelmag[3,*]
   dperp = (ricolor) - (grcolor)/8.d0
   return, dperp

end


function lrg_cpllel, modelmag

   grcolor = modelmag[1,*]-modelmag[2,*]
   ricolor = modelmag[2,*]-modelmag[3,*]
   c= 0.7
   cpllel = c*(grcolor) + (1.d0-c)*4.0d0*((ricolor) - 0.18)
   return, cpllel

end

function lrg_cperp, modelmag

   grcolor = modelmag[1,*]-modelmag[2,*]   
   ricolor = modelmag[2,*]-modelmag[3,*]
   cperp = (ricolor) - (grcolor/4.d0)-0.18
   return, cperp
end


function boss_galaxy_select, calibobj

     ; Get the number of galaxies
     ngal = n_elements(calibobj)

     ; Set some reasonable defaults
     maglim=[13.4d0, 19.5d0]  ; Cut I parameters
     dperp0=0.55              ; Cut II dperp
     dperp1=0.75              ; Cut III dperp
     ilimit=19.5              ; Cut II faint limit
     ilow=17.5                ; Cut II, III bright limit
     ihi=19.9                ; Cut III faint limit

     ; Check to see if calibobj has any elements
     if (ngal EQ 0) then $
       message, 'ERROR : calibobj has no elements'

     ; Extract the relevant magnitudes
     ; We really should have used cmodelmags, and will later.
     modelmag = 22.5 - 2.5*alog10(calibobj.modelflux > 0.001)
     devmag = 22.5 - 2.5*alog10(calibobj.devflux > 0.001)
     expmag = 22.5 - 2.5*alog10(calibobj.expflux > 0.001)
     cmodelflux = calibobj.devflux*calibobj.fracpsf + calibobj.expflux*(1.0-calibobj.fracpsf)
     cmodelmag = 22.5-2.5*alog10(cmodelflux > 0.001)
     fibermag = 22.5 - 2.5*alog10(calibobj.fiberflux > 0.001)
     psfmag = 22.5 - 2.5*alog10(calibobj.psfflux > 0.001)
     ; Extinction correct
     modelmag = modelmag - calibobj.extinction
     psfmag = psfmag - calibobj.extinction
     cmodelmag = cmodelmag - calibobj.extinction
     fibermag = fibermag - calibobj.extinction
     
     ; Compute cperp and cpllel
     cperp = lrg_cperp(modelmag)
     dperp = lrg_dperp(modelmag)
     cpllel = lrg_cpllel(modelmag)     

     ; Cut I
     icut1 = cmodelmag[2,*] LT (maglim[0] + cpllel/0.3d0)
     icut1 = icut1 AND (abs(cperp) LT 0.2d0)
     icut1 = icut1 AND (cmodelmag[2,*] LT maglim[1]) 
     ilrg = icut1 * 2L^1

     ; Cut II
     icut2 = (cmodelmag[3,*] GT ilow) AND (cmodelmag[3,*] LT ilimit)
     icut2 = icut2 AND ((modelmag[2,*] - modelmag[3,*]) LT 2.0d0)
     icut2 = icut2 AND (dperp GT dperp0) AND (dperp LE dperp1)
     icut2 = icut2 AND (fibermag[3,*] LT 21.3)
     ilrg  = ilrg + icut2*2L^2

     ; Cut III
     icut3 = (cmodelmag[3,*] GT ilow) AND (cmodelmag[3,*] LT ihi)
     icut3 = icut3 AND ((modelmag[2,*] - modelmag[3,*]) LT 2.0d0)
     icut3 = icut3 AND (dperp GT dperp1)
     icut3 = icut3 AND (fibermag[3,*] LT 21.3)
     ilrg  = ilrg + icut3*2L^3

     ; Only work with those objects that photo calls a galaxy, and don't
     ; have processing flags thrown.
     icut3 = reform(icut2*0b)
     indx = sdss_selectobj(calibobj, objtype=3, /trim) 
     if (indx[0] GT -1) then $      
       icut3[indx] = 1
     ilrg = icut3*reform(ilrg) 

     return, ilrg
end
