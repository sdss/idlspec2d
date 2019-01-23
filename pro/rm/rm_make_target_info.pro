;+
; NAME:
;   rm_make_target_info
;
;
; PURPOSE:
;   compile some basic quasar properties for the RM targets
;

pro rm_make_target_info

   target_file = getenv('IDLRM_DIR') + '/etc/target_fibermap.fits'
   fibermap = mrdfits(target_file,1,/silent)
   ; keep only RM targets
   target=fibermap[0:848]
   nobj=n_elements(target)

   struct=replicate({imag:0.D, Mi_z2:0.D, logL5100:0.D, tau_obs:-1.D}, nobj)
   target=struct_addtags(target,struct)

   target.imag=reform((target.psfmag)[3,*])
   mi_z2=get_abs_mag(target.imag, target.zfinal)
   target.mi_z2=mi_z2
   d10pc = 3.08d19
   logL2500 = -0.4*(mi_z2 + 48.6 + 2.5*alog10(3.)) + alog10(4d*!PI*d10pc^2) + alog10(3d18/2500d)
   logL5100 = logL2500 - 0.5*alog10(5100./2500.)
   target.logL5100 = logL5100

   ; assign BLR size (hbeta BLR), in days
   lag = 10.D^( -21.3 + 0.519*target.logL5100 ) ; + randomn(seed, nobj)*0.15 )
   target.tau_obs = lag*(1. + target.zfinal)

   outfile='/data3/quasar/yshen/work/composite_lag/target_info.fits'
   mwrfits, target, outfile, /create

end
