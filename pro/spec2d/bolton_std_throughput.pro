;+
;
; NAME:
;  bolton_std_throughput
;
; PURPOSE:
;  Compute BOSS end-to-end throughput estimate from
;  spectrophotometric standard stars
;
; USAGE:
;   tput_struc = bolton_std_throughput(fileid, $
;    plate=plate, mjd=mjd [, radius=radius, plug=plug])
;
; ARGUMENTS:
;  Required:
;   fileid: a file identifier of the 'b1-00124616' variety
;           (no 'sdR-', no '.fit*') for a science frame.
; (one of the following two is mandatory:)
;   plate: plate, to give subdir of spCFrame file
;   mjd: mjd, to give subdir of sdR file
;  Optional:
;   radius: pixel radius for boxcar extraction
;
; RETURNS:
;  a structure containing wavelength and end-to-end
;  throughput information calculated for all standard stars
;  in the observation.
;
; Optional output:
;  plug: plugmap structure for the standard stars.
;
; NOTES:
;  This procedure uses sdR files and spCFrame files.
;  It needs to find them via the environment variables
;    $BOSS_SPECTRO_DATA
;    $BOSS_SPECTRO_REDUX
;    $RUN2D
;
;  Calculated throughput is total end-to-end, incl. atmosphere
;  through CCD QE, also including any seeing-related fiber loss.
;  Telescope is 2.5m with central area of 27% obscured.
;
; WRITTEN:
;  A. Bolton, U. of Utah, 2011 Feb.
;
;-

function bolton_std_throughput, fileid, plate=plate, mjd=mjd, $
 radius=radius, plug=plug

if (not keyword_set(radius)) then radius = 3.0

; Find the files.
; If plate passed, get MJD from spCFrame header.
; if MJD passed, get plate from sdR header.
if keyword_set(plate) then begin
    ; spCFrame:
    cframe = getenv('BOSS_SPECTRO_REDUX') + '/' $
      + getenv('RUN2D') + '/' + string(plate, format='(i4)') + $
      '/spCFrame-' + fileid + '.fit*'
    cframe = (file_search(cframe))[0]
    if (cframe eq '') then begin
        splog, 'spCFrame file not found for ' + fileid
        return, 0
    endif
    mjd = sxpar(headfits(cframe), 'MJD')
    ; sdR:
    sciframe = getenv('BOSS_SPECTRO_DATA') + '/' $
      + string(mjd, format='(i5)') + $
      '/sdR-' + fileid + '.fit*'
    sciframe = (file_search(sciframe))[0]
    if (sciframe eq '') then begin
        splog, 'sdR file not found for ' + fileid
        return, 0
    endif
endif else if keyword_set(mjd) then begin
    ; sdR:
    sciframe = getenv('BOSS_SPECTRO_DATA') + '/' $
      + string(mjd, format='(i5)') + $
      '/sdR-' + fileid + '.fit*'
    sciframe = (file_search(sciframe))[0]
    if (sciframe eq '') then begin
        splog, 'sdR file not found for ' + fileid
        return, 0
    endif
    plate = sxpar(headfits(sciframe), 'PLATEID')
    ; spCFrame:
    cframe = getenv('BOSS_SPECTRO_REDUX') + '/' $
      + getenv('RUN2D') + '/' + string(plate, format='(i4)') + $
      '/spCFrame-' + fileid + '.fit*'
    cframe = (file_search(cframe))[0]
    if (cframe eq '') then begin
        splog, 'spCFrame file not found for ' + fileid
        return, 0
    endif
endif else begin
    splog, 'Need either plate or MJD to find file!'
    return, 0
endelse

; Get the stuff we need from the sdR file:
sdssproc, sciframe, image, invvar, hdr=i_hdr, /applybias, applypixflat=0

; Get the stuff we need from the spCFrame file:
flux = mrdfits(cframe, 0, c_hdr)
ivar = mrdfits(cframe, 1)
loglam = mrdfits(cframe, 3)
plug = mrdfits(cframe, 5)
sky = mrdfits(cframe, 6)
xpos = mrdfits(cframe, 7)

; Downselect to standard stars.
; Also get rid of blue-side spCFrame row padding, if necessary.
npix = (size(image))[2]
wh_std = where(strtrim(plug.objtype,2) eq 'SPECTROPHOTO_STD', n_std)
plug = plug[wh_std]
flux = flux[0:npix-1,wh_std]
ivar = ivar[0:npix-1,wh_std]
loglam = loglam[0:npix-1,wh_std]
sky = sky[0:npix-1,wh_std]
xpos = xpos[0:npix-1,wh_std]
ypos = findgen(npix) # replicate(1., n_std)

; Extract raw counts:
c_extract = extract_boxcar(image * (invvar gt 0.), xpos, ypos, radius=radius)

; Convert fluxed spectrum from flux per angstrom to flux per pixel:
wave = 10.^loglam
dwave = 0.5 * (wave[2:npix-1,*] - wave[0:npix-3,*])
dwave = [dwave[0,*], dwave, dwave[npix-3,*]]
flux_per_pix = (flux + sky) * dwave

; Collecting area of SDSS telescope in cm^2
; (27% central obscuration of area as per project book):
t_area = (1. - 0.27) * !pi * (250. / 2.)^2
; Exposure time:
exptime = sxpar(c_hdr, 'EXPTIME')

; Convert fluxed spectrum into ergs per pixel:
energy_per_pix = flux_per_pix * t_area * exptime * 1.e-17

; Convert energy into photons:
h_cgs = 6.626e-27
c_cgs = 2.998e10
lambda_cgs = wave * 1.e-8
energy_per_phot = h_cgs * c_cgs / lambda_cgs
photons_per_pix = energy_per_pix / energy_per_phot

; Compute throughput:
tput = (photons_per_pix gt 0.) * (ivar gt 0.) * c_extract / (photons_per_pix > 1.)

; Pack up and return:
ostruc = {plate: plate, fileid: fileid, exptime: exptime, wave: wave, tput: tput}

return, ostruc
end
