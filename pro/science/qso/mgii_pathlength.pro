
function mgii_pathlength, summary, mgspec, ew_limits=ew_limits, sigma=sigma, $
       binlnz=binlnz

;
;   These ew_limits are in **PIXELS**, for SDSS and MgII these
;   correspond to rest EW of 0.32, 0.64, 0.97 Ang
;
    if NOT keyword_set(ew_limits) then ew_limits = [0.5, 1.0, 1.5, 2.0]
    number_ew = n_elements(ew_limits)

    if NOT keyword_set(sigma) then sigma=5.0
    if NOT keyword_set(binlnz) then binlnz = 0.05
    al = alog(10)

    template = { minlnz : 0.0, maxlnz : 0.0, $
                 n_abs : ew_limits*0.0, pathlength : ew_limits*0.0 } 

    lnrest = mgspec.loglam * al - alog(2803.531)
    maxwave = max(lnrest * (mgspec.mask))
    minwave = min(lnrest + (9999.0 * (1 - mgspec.mask)))

    minlnz = fix(minwave / binlnz) * binlnz
    nlnz   = fix(maxwave / binlnz) - fix(minwave / binlnz) + 1
    maxlnz = minlnz + nlnz*binlnz

    results = replicate(template, nlnz)
    results.minlnz = minlnz + findgen(nlnz)*binlnz
    results.maxlnz = minlnz + (findgen(nlnz)+1.0)*binlnz

    summary_lnz = alog(1.0+summary.zabs)
    for i=0,nlnz-1 do begin
       
      ; find good pixels :
      pass2803 = where(lnrest GE results[i].minlnz AND $
                       lnrest LT results[i].maxlnz)
      if pass2803[0] EQ -1 then continue

      ewinv = (mgspec.finv)[pass2803] < (mgspec.finv)[pass2803-11]

      ewtest = sqrt(ewinv) # (ew_limits/sigma) 

;
;	Sum total pathlength in each spectrum in units of ln(1.0+z)
;       Each spectrum can contribute at most binlnz in pathlength
;
      results[i].pathlength = total(ewtest GE 1.0, 1) * 1.0e-4 * al 

      nabstest = (((summary_lnz GE results[i].minlnz) AND $
                 (summary_lnz LT results[i].maxlnz)) * summary.ew2803) $
                 # (1.0/ew_limits)
      results[i].n_abs = total(nabstest GE 1.0, 1)
    endfor

    return, results
end


