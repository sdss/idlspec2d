function configuration::init, mjd
  self.mjd=mjd
  self.sdss3_start_mjd=55000
  return,1
end

function configuration::getMJD
  return, self.mjd
end


function configuration::getNumberFibersPerSpectrograph
  return, self->getNumberBundles() * self->getNumberFibersPerBundle()
end

function configuration::getNumberFibersTotal
  return, self->getNumberFibersPerSpectrograph() * self->getNumberSpectrographs()
end

pro configuration::toString
  print, self.mjd, self.sdss3_start_mjd
end

function configuration::cleanup
  return,1
end

;
;  Methods that contain parameters that are expected to be static for
;  the lifetime of this code
;
function configuration::getNumberFibersPerBundle
  return, 20
end

function configuration::getNumberSpectrographs
  return,2
end

;
; The methods that contain the special configurations that must
; be specified
;


; l395 need new keyword in fitdispersion
function configuration::getNumberBundles
  if self->isSDSS2() then return, 16
  return, 25
end

;+
;return the detector format.  blue is assumed for anything that
;starts with b so 'blue', 'b1' or 'b2'.  likewise for red.
;-
function configuration::getDetectorFormat, color
  if self->isSDSS2() then return, [2048,2048]
  if self->isLaurenSimulation_v1() then return, [4096,4096]
  
  if (strmid(color,0,1) EQ 'b') then return, [4112, 4096]
  if (strmid(color,0,1) EQ 'r') then return,  [4128, 4114]
end

;
; Changes in spcalib.pro, tag v5_4_0
;

;l191 xy2traceset keyword
function configuration::spcalib_xy2traceset_ncoeff
  if self->isSDSS2() then return,7
  return, 7 ; ASB: 4 too few for BOSS
end


;l194
function configuration::spcalib_rejecttheshold
  if self->isSDSS2() then return,10000
  return,  1000000
end

;l228
;l395
function configuration::spcalib_sigmaguess
  if self->isSDSS2() then return,1.0
  return,  1.15 ; Compromise value, ASB 2010aug
end

;l236 extract_image
;l244
function configuration::spcalib_extract_image_proftype
  if self->isSDSS2() then return,3
  return,  1
end

;l241 fitflatwidth
function configuration::spcalib_fitflatwidth_ncoeff
  if self->isSDSS2() then return,5
  return,  3
end

;l241
function configuration::spcalib_arccoeff
  if self->isSDSS2() then return,5
  return,  6
end

;l393
function configuration::spcalib_ncoeff, color
  if self->isSDSS2() then return,color EQ 'red' ? 4 : 3
  return,  2
end

;l519
function configuration::spcalib_doscatter_maxkernr, color
  if self->isSDSS2() then begin
    if color EQ 'r' then return, 256
    if color EQ 'b' then return, 64
  endif else begin
    if color EQ 'r' then return, 100*4
    if color EQ 'b' then return, 100
  endelse
end
;l519
function configuration::spcalib_doscatter_addrows, color
  if self->isSDSS2() then begin
    if color EQ 'r' then return, 64
    if color EQ 'b' then return, 32
  endif else begin
    if color EQ 'r' then return, 100
    if color EQ 'b' then return, 50
  endelse
end

;l555 set keyword for fiberflat
function configuration::spcalib_fiberflat_minval,flux
  if self->isSDSS2() then return, 0.03 * median(flux)
  return, 0.05 * median(flux)
end

;l555 set keyword for fiberflat
function configuration::spcalib_fiberflat_badflatfracthresh
  if self->isSDSS2() then return, 0.7
  return, 0.5
end

;l241 set keyword for fitflatwidth
function configuration::spcalib_fitflatwidth_mask,flux,fluxivar
  if self->isSDSS2() then return, (flux GT 0) * (fluxivar GT 0)
  return, flux * sqrt(fluxivar) GT 10
end

;l241 set keyword for fitflatwidth
function configuration::spcalib_fitflatwidth_inmask,flux,fluxivar
  if self->isSDSS2() then return, 0
  return, reform(self->spcalib_fitflatwidth_mask(flux,fluxivar),(size(flux,/dimen))[0],20*self->getNumberBundles())
end


;l175 set keyword for reject_flact
function configuration::spcalib_reject_calib_percent80thresh
  if self->isSDSS2() then return, 1000.
  if self->isFaintflat() then return, 500.
  return, 900.
end


;l175 set keyword for reject_flact
function configuration::spcalib_mthreshFactor
  if self->isFaintflat() then return, .2
  return, .5
end

;l186 set keyword for trace320crude
function configuration::spcalib_trace320crude_padding
  if self->isSDSS2() then return, 0
  return, 10
end


;l186 set keyword for trace320crude
; within trace320crude l90 for trace320cen
function configuration::spcalib_trace320cen_deltax, color, spectrographid
  if self->isSDSS2() then return, 6.25 ; SDSS-I

  ; For BOSS, return the fiber spacing per bundle
  nbundle = 25
  if (color EQ 'blue' AND spectrographid EQ 1) then $
;   return, 6.605 + 0.05*abs(2*(findgen(nbundle)-0.5*nbundle)/nbundle)
   return, 6.545 + 0.12*abs(2*(findgen(nbundle)-0.5*nbundle)/nbundle)
  if (color EQ 'blue' AND spectrographid EQ 2) then $
   return, 6.543 + 0.11*abs(2*(findgen(nbundle)-0.5*nbundle)/nbundle)
  if (color EQ 'red' AND spectrographid EQ 1) then begin
     ; The CCD was moved and scale changed on MJD 55168 (Dec 4, 2009)
     if (self.mjd lt 55168) then $
      return, 6.653 + 0.15*abs(2*(findgen(nbundle)-0.5*nbundle)/nbundle)
     return, 6.551 + 0.15*abs(2*(findgen(nbundle)-0.5*nbundle)/nbundle)
  endif
  if (color EQ 'red' AND spectrographid EQ 2) then $
   return, 6.595 + 0.16*abs(2*(findgen(nbundle)-0.5*nbundle)/nbundle)
end

function configuration::spcalib_trace320cen_xstart, color, spectrographid
  if self->isSDSS2() then return, 0. ; SDSS-I

  if (color EQ 'blue' AND spectrographid EQ 1) then return, 266.
  if (color EQ 'blue' AND spectrographid EQ 2) then return, 288.
  if (color EQ 'red' AND spectrographid EQ 1) then begin
    if (self.mjd lt 55168) then return, 245.
    return, 271.
  endif
  if (color EQ 'red' AND spectrographid EQ 2) then return, 259.
end


;l186 set keyword for trace320crude
; within trace320crude l90 for trace320cen
function configuration::spcalib_trace320cen_bundlebreakrange
  if self->isSDSS2() then return, [1,1.6]
  return, [1,2.7]
end

;l186 set keyword for trace320crude
; within trace320crude l90 for trace320cen
function configuration::spcalib_trace320cen_nextpeakrange
  if self->isSDSS2() then return, [0.35,0.25]
  return, [0.45,0.25]
end

;l186 set keyword for trace320crude
; within trace320crude l90 for trace320cen
function configuration::spcalib_trace320cen_bundlegap
  if self->isSDSS2() then return, 0.28
  return, 1.40
end

;l363,385 set keyword for fitarcimage
; within fitarcimage l90 for arcfit_guess
function configuration::spcalib_arcfitguess_acoeff,color
  if (self->isSDSS2()) then begin
    if (color EQ 'blue') then return, [3.6930, -0.1044, -0.0041, 0.00016]
    if (color EQ 'red') then return, [ 3.8700, 0.1008, -0.0044, -0.00022]
    return,-1
  endif
  if (self->isLaurenSimulation_v1() or self->isLaurenSimulation_v2()) then begin
    if (color EQ 'blue') then return, [3.6842,0.2002,-0.0242,0.0012]
    if (color EQ 'red') then return, [ 3.8991,0.1613,-0.0178,0.0012]
    return,-1
  endif
  ; all other cases
; The following are for the BOSS commissioning data w/out rotating the images
;  if (color EQ 'blue') then return, [3.6846659,-0.19959101, -0.023845067, -0.0012433996]
;  if (color EQ 'red') then return, [3.8995092,-0.15889163, -0.017838774, -0.00075860025]
  if (color EQ 'blue') then return, [3.6846659,0.19959101, -0.023845067, 0.0012433996]
  if (color EQ 'red') then return, [3.8927011, 0.16506676, -0.018754040, 0.0011861915, -1.1795306e-05, -2.4502617e-05]
  return,-1
end

;l363,385 set keyword for fitarcimage
; within fitarcimage l90 for arcfit_guess
function configuration::spcalib_arcfitguess_dcoeff,color
  if self->isSDSS2() then return, 0
  return, [0.15, 0.008, 0.0003, 0.0001, 1e-5, 1e-5]
end

;l363,385 set keyword for fitarcimage
;function configuration::spcalib_fitarcimage_trimrange,arc,color
;  if self->isSDSS2() then begin
;    npix= (size(arc, /dim))[0]
;    return, [1,npix-2]
;  endif
;  if (self->isLaurenSimulation_v1() or self->isLaurenSimulation_v2()) then begin
;    if (color EQ 'blue') then return, [1075,3110]
;    if (color EQ 'red') then return, [490,2850]
;  endif
;  if (color EQ 'blue') then return, [870,3350]
;  if (color EQ 'red') then return, [360,3600]
;end

;l363,385 set keyword for fitarcimage
function configuration::spcalib_fitarcimage_wrange,color
  if self->isSDSS2() then return,0
  if (color EQ 'blue') then return, [3500.,6540.]
  if (color EQ 'red') then return, [5400.,10500.]
end


function configuration::extract_object_skyfit2dbspline
  if self->isSDSS2() then return,1
  return, 0
end

;l237 extract_object
function configuration::getscatter, camera, flatimg, flativar, wset, sigma=sigma
  if self->isSDSS2() then return, doscatter(camera, flatimg, flativar, wset, sigma=sigma)
  return, 0.*flatimg
end

;
function configuration::extract_object_skysubtract_alternative
  if self->isSDSS2() then return, 0
  return, 1
end

;l188,189
function configuration::extract_object_fixcolumns
  if self->isSDSS2() then return, 1
  return, 0
end

;l609 sdssproc
function configuration::sdssproc_skipprocessrawdata
  if self->isLaurenSimulation_v1() then return,1
  return,  0
end


function configuration::spflux_v5_minfracthresh
  if self->isSDSS2() then return,0.80
  return,  0.50
end

function configuration::crosstalk_ampgrid, camera
  if self->isSDSS2() then return, [1,2]
  return, [2,2]
end

function configuration::crosstalk_crosstalk, camera

  if self->isSDSS2() then return, 0
  ;
  ;
  ; The cross-talk parameters for each camera and amp
  ; the formation is xtalk[i,j] where i is the stimulus and j is the affected amp
  ;
  ; r1
  r1crosstalk = fltarr(4,4)
  r1crosstalk[0,*] = [1, 0.000959340, 0.000680067, 7.48254e-05]
  r1crosstalk[1,*] = [5.30889e-05, 1, 4.16112e-06, 0.000394995]
  r1crosstalk[2,*] = [9.45941e-06, -7.76818e-07, 1, -3.52740e-06]
  r1crosstalk[3,*] = [-0.000120624, -0.000259071, 8.34742e-06, 1]
  
  ;r2
  r2crosstalk = fltarr(4,4)
  r2crosstalk[0,*] = [1, 0.00127615, 0.000823262, 4.27848e-05]
  r2crosstalk[1,*] = [-7.01562e-05, 1, -1.80703e-05, 0.000393628]
  r2crosstalk[2,*] = [5.04834e-05, -4.20647e-05, 1, 2.25005e-05]
  r2crosstalk[3,*] = [0.000170907, 9.27249e-05, -6.91713e-05, 1]
  
  ;b1
  b1crosstalk = fltarr(4,4)
  b1crosstalk[0,*] = [1, 0.000126933,0.000103371,1.32322e-05]
  b1crosstalk[1,*] = [5.07551e-05, 1, 1.39503e-05, 0.000108620]
  b1crosstalk[2,*] = [0.000120932, 1.30320e-05, 1, 5.67557e-05]
  b1crosstalk[3,*] = [4.44612e-06, 6.58372e-05, 0.000105030, 1]
  
  ;b2
  b2crosstalk = fltarr(4,4)
  b2crosstalk[0,*] = [1, 0.000189313, 0.000970775, 3.42132e-05]
  b2crosstalk[1,*] = [8.45262e-05, 1,4.05977e-05, 0.00295630]
  b2crosstalk[2,*] = [0.000140322,3.57650e-06, 1, 8.10733e-05]
  b2crosstalk[3,*] = [6.28033e-06,0.000173701, 0.000161124, 1]
  
  if (camera eq 'r1') then return, r1crosstalk
  if (camera eq 'r2') then return, r2crosstalk
  if (camera eq 'b1') then return, b1crosstalk
  if (camera eq 'b2') then return, b2crosstalk
  
end


;
; put the different possibilities here
;
function configuration::isSDSS2
  return, self.mjd lt self.sdss3_start_mjd
end

function configuration::isLaurenSimulation_v1
  return, self.mjd eq 55025
end

function configuration::isLaurenSimulation_v2
  return, self.mjd eq 55026
end

function configuration::isFaintFlat
  return, (self.mjd eq 55052 or self.mjd eq 55070)
end

function configuration::isSDSS3
end

pro configuration__define
  void={configuration, mjd:0L, sdss3_start_mjd:0.D}
  return
end
