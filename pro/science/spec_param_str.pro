;+
; NAME:
;   spec_param_str
; PURPOSE:
;   Build structure SPEC_PARAM for parameters derived from analysis of
;     galaxy spectra.
; CALLING SEQUENCE:
;   str= spec_param_str(N)
; INPUTS:
;   tsObjfile  - tsObj file name
;   ver1Ddir  - 1D directory name (including version number)
; OPTIONAL INPUTS:
; KEYWORDS
; OUTPUTS:
;   str     - array of structures of type "SPEC_PARAM" with N elements
; COMMENTS:
; BUGS: 
;   heliocentric correction will have to be removed when the 2D bug is fixed
; REVISION HISTORY:
;   2000-Sep-07  written by Mariangela Bernardi
;-
;------------------------------------------------------------------------------
function spec_param_str, tsObjfile, ver1Ddir, list = list

  if (keyword_set(list)) then begin
   inputlist = tsObjfile 
   readcol, inputlist, tsObjfile, format='(a)'
  endif
  
  nfiles=n_elements(tsObjfile)

  for j=0,nfiles-1 do begin

  photo = mrdfits(tsObjfile[j],1) 
 
  nchar=strlen(tsObjfile[j]) - 19
  splate =strmid(tsObjfile[j], nchar+12, 3)
  smjd = strmid(tsObjfile[j], nchar+6, 5)
  plate= long(splate)
  mjd = long(smjd)

  ngals=640
  ntemp= 32

  tempearly = {SPEC_PARAM, vers_spec: ' ', $
                       plate: 0L, mjd: 0L, fiberid: 0L, $
                       run:0L, rerun: 0L, camcol: 0L, field: 0L, id: 0L, $
                       parent: 0L, $ 
                       specobjid: lonarr(5), $
                       ra: 0.0, dec: 0.0, $  
                       l: 0.0, b: 0.0, $        
                       objc_rowc: 0.0, objc_colc: 0.0, $
                       psfcounts: fltarr(5), psfcountserr: fltarr(5), $
                       fracpsf: fltarr(5), $
                       fibercounts: fltarr(5), fibercountserr: fltarr(5), $ 
                       petrocounts: fltarr(5), petrocountserr: fltarr(5), $
                       petrorad: fltarr(5), petroraderr: fltarr(5), $ 
                       petror50: fltarr(5), petror50err: fltarr(5), $ 
                       petror90: fltarr(5), petror90err: fltarr(5), $  
                       counts_model: fltarr(5), counts_modelerr: fltarr(5), $
                       counts_dev: fltarr(5), counts_deverr: fltarr(5), $
                       counts_exp: fltarr(5), counts_experr: fltarr(5), $
                       dev_l:  fltarr(5), exp_l: fltarr(5), $
                       r_dev: fltarr(5), r_deverr: fltarr(5), $
                       r_exp: fltarr(5), r_experr: fltarr(5), $
                       ab_dev: fltarr(5), ab_deverr: fltarr(5), $
                       ab_exp: fltarr(5), ab_experr: fltarr(5), $
                       reddening: fltarr(5), $
                       primtarget:0L, sectarget:0L, $
                       flags: lonarr(5), flags2: lonarr(5), $ 
                       firstint: 0.0, $
		       type: ' ', $  	
                       primtar: 0L, sectar: 0L, $
                       z_1Dclass: 0L, z_1Dstatus: 0L, $
                       z_1Dwarning: 0L,  $
                       z_1D: 0.0, z_1Derr: 0.0, z_1Dconf: 0.0, $
                       sn_g: 0.0, sn_r: 0.0, sn_i: 0.0, $
                       z: fltarr(ntemp), z_err: fltarr(ntemp), zconf: fltarr(ntemp), $
                       sigma_cc: intarr(ntemp), sigma_ccerr: intarr(ntemp), $
                       sigma_diff: intarr(ntemp), sigma_differr: intarr(ntemp), $ 
                       sigma_quotient: intarr(ntemp), sigma_quotienterr: intarr(ntemp), $
                       sigma_quot_diff: intarr(ntemp), sigma_quot_differr: intarr(ntemp), $
                       sigma_pc: intarr(ntemp), sigma_pcerr: intarr(ntemp), $
                       lickindex: fltarr(21), lickindexerr: fltarr(21) }
 

  sample = replicate(tempearly, ngals)

  id=lonarr(3)
  teta=fltarr(640)

  for i=0,ngals-1 do begin

         ii=i+1
         ff = string('/usr/sdss/data05/spectro/'+ver1Ddir+'/'+string(plate,$
         format='(i4.4)')+$
        '/1d/spSpec-'+string(mjd,format='(i5.5)')+'-'$
         +string(plate,format='(i4.4)')+'-'$
         +string(ii,format='(i3.3)')+'.fit')

         spec = mrdfits(ff,0,hdr)          

         if ((size(spec))[2] GT 1) then begin
   
            type=sxpar(hdr,'OBJTYPE')
            id[0]= sxpar(hdr,'PLATEID')
            id[1]= sxpar(hdr,'MJD')
            id[2]= sxpar(hdr,'FIBERID')
            sobjid= sxpar(hdr,'OBJID')
            sobjid=str_sep(sobjid,'         ')
            nobjid=n_elements(sobjid)-1
            objid=long(sobjid[1:nobjid])	  
            ra=sxpar(hdr,'RAOBJ')
            dec=sxpar(hdr,'DECOBJ')
            primtar= sxpar(hdr,'PRIMTARG')
            sectar= sxpar(hdr,'SECTARG')
            zclass=sxpar(hdr,'SPEC_CLN')
            z=sxpar(hdr,'Z')
            zerr= sxpar(hdr,'Z_ERR')
            zconf= sxpar(hdr,'Z_CONF')
            zstat= sxpar(hdr,'Z_STATUS')
            zwarn= sxpar(hdr,'Z_WARNIN')
            sng=sxpar(hdr,'SN_G')
            snr=sxpar(hdr,'SN_R')
            sni= sxpar(hdr,'SN_I')
            helio2= sxpar(hdr,'HELIO_RV') 
            helio1= sxpar(hdr,'EARTH_RV')

            rar=double(ra * !PI/180.0)
            decr=double(dec * !PI/180.0)
          
         for kk=0,ngals-1 do begin         
            raphoto=double(photo[kk].ra * !PI/180.0)
            decphoto=double(photo[kk].dec * !PI/180.0)
            teta[kk]=acos(sin(decr)*sin(decphoto)+cos(decr)*cos(decphoto)*$
                 cos(rar-raphoto))*180/!PI*60*60.
         endfor

            good=where(teta LE 1)
 
           if (good[0,0] EQ -1) then begin
           print,'warning: no objects identified' 
           endif else begin
    
            if (n_elements(good) GT 1) then begin 
            print,'warning: two objects identified; choosing closest'
            junk=min(teta[good],best)
            good=good[best]  
            endif

;;;;;  WARNING: hel. redshift correction
 
            c=299792
            z=z+(helio1-2*helio2)/c 

            sample[i].vers_spec = ver1Ddir
            sample[i].plate=id[0]
            sample[i].mjd=id[1]
	    sample[i].fiberid=id[2]
            sample[i].run = photo[good].run
            sample[i].rerun = photo[good].rerun
            sample[i].camcol = photo[good].camcol
            sample[i].field = photo[good].field
            sample[i].id = photo[good].id
            sample[i].parent = photo[good].parent
            sample[i].specobjid= objid[*]
            sample[i].ra=ra
            sample[i].dec=dec
            sample[i].l= photo[good].l
            sample[i].b=photo[good].b
            sample[i].objc_rowc = photo[good].objc_rowc
            sample[i].objc_colc = photo[good].objc_colc
            sample[i].psfcounts = photo[good].psfcounts
            sample[i].psfcountserr = photo[good].psfcountserr
            sample[i].fracpsf=photo[good].fracpsf
            sample[i].fibercounts = photo[good].fibercounts
            sample[i].fibercountserr = photo[good].fibercountserr
            sample[i].petrocounts = photo[good].petrocounts
            sample[i].petrocountserr = photo[good].petrocountserr
            sample[i].petrorad = photo[good].petrorad
            sample[i].petroraderr = photo[good].petroraderr
            sample[i].petror50 = photo[good].petror50
            sample[i].petror50err = photo[good].petror50err
            sample[i].petror90 = photo[good].petror90
            sample[i].petror90err = photo[good].petror90err
            sample[i].counts_model = photo[good].counts_model
            sample[i].counts_modelerr = photo[good].counts_modelerr
            sample[i].counts_dev = photo[good].counts_dev
            sample[i].counts_deverr = photo[good].counts_deverr
            sample[i].counts_exp = photo[good].counts_exp
            sample[i].counts_experr = photo[good].counts_experr
            sample[i].dev_l = photo[good].dev_l
            sample[i].exp_l = photo[good].exp_l
            sample[i].r_dev = photo[good].r_dev
            sample[i].r_deverr = photo[good].r_deverr
            sample[i].r_dev = photo[good].r_exp
            sample[i].r_deverr = photo[good].r_experr   
            sample[i].ab_dev = photo[good].ab_dev
            sample[i].ab_deverr = photo[good].ab_deverr
	    sample[i].ab_exp = photo[good].ab_exp
            sample[i].ab_experr = photo[good].ab_experr
            sample[i].reddening = photo[good].reddening
            sample[i].primtarget = photo[good].primtarget
            sample[i].sectarget = photo[good].sectarget 
            sample[i].flags = photo[good].flags
            sample[i].flags2 = photo[good].flags2 
            sample[i].firstint = photo[good].firstint
            sample[i].type=type 
            sample[i].primtar= primtar
            sample[i].sectar= sectar
            sample[i].z_1Dclass=zclass
            sample[i].z_1Dstatus=zstat
            sample[i].z_1Dwarning=zwarn
            sample[i].z_1D=z
            sample[i].z_1Derr=zerr 
            sample[i].z_1Dconf=zconf
            sample[i].sn_g=sng
            sample[i].sn_r=snr 
            sample[i].sn_i=sni
            
       endelse

    endif          

  endfor

  namesample=string('spec_param_'+string(mjd,format='(i5.5)')+'_' $
                    +string(plate, format='(i4.4)')+'.fit')
  mwrfits, sample, namesample, /create

    if (size(fullsample,/tname) EQ 'STRUCT') then begin
         fullsample = [fullsample, sample]
    endif else begin
         fullsample = sample
    endelse

  endfor

return, fullsample

end
