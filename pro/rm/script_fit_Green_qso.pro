; fit Paul Green's MD fields QSO

pro script_fit_Green_qso, tag=tag, ntrial=ntrial, comp_only=comp_only, _extra=extra

; set /comp_only to complie the fitting results

topdir = '/data3/yshen/work/md_qso/'
file = topdir + 'var_qso.fits'
varqso=mrdfits(file,1)

; set up the qsofit file
emparfile=getenv('IDLRM_DIR')+'/etc/qsoline_qsovar.par'

linename=['Halpha_br', 'Hbeta_br', 'MgII_br', 'CIII_br', 'CIV_br', 'OIII5007', 'OIII5007c', $
             'OIII5007w', 'SiIII1892', 'AlIII1857', 'HeII1640', 'OIII1663']
contiwave=[1350.,1700.,3000.,5100.]

;tag='MD01','MD03', 'MD_DR10', 'S82_DR10'
if tag eq 'MD07' then begin 
   ; this is the RM field
   coadd_mjd = 56783
   outdir=topdir + tag + '/'
   outfile='qso_prop-' + tag + '.fits'
   figfile='QA-'+ tag + '-0000-' + string(coadd_mjd,format='(i5.5)') + '.ps'
   if ~keyword_set(comp_only) then begin
       rm_fitplate,0,coadd_mjd,/rm_plate,emparfile=emparfile,outdir=outdir,/coadd, $
        linename=linename,contiwave=contiwave,/add_bhmass,outfile=outfile,figfile=figfile,ntrial=ntrial
   endif else $
       rm_compile_qsofit,0,coadd_mjd, /rm_plate, outdir_in=outdir,outfile=outfile, $
        linename=linename,contiwave=contiwave,/add_bhmass
endif else begin
   ind=where(strtrim(varqso.sample) eq tag)
   outdir=topdir + tag + '/'
   outfile='qso_prop-' + tag + '.fits'
   if ~keyword_set(comp_only) then begin
       fitqso_batch,varqso[ind].plate,varqso[ind].fiber,varqso[ind].mjd, $
          zfit=varqso[ind].z_pipe,outdir=outdir,tag=tag,ra=varqso[ind].ra,dec=varqso[ind].dec, $
          emparfile=emparfile,ntrial=ntrial, _extra=extra
   endif else begin
       rm_compile_qsofit,varqso[ind].plate,varqso[ind].mjd, input_fiber=varqso[ind].fiber, $
         outdir_in=outdir,outfile=outfile,rm_plate=0, $
         linename=linename,contiwave=contiwave,/add_bhmass
   endelse
endelse

end
