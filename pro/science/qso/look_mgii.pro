
function look_mgii, this_list, spec, index, pixel_sep, $
        restwave, f_value, label, acc=acc
  
      if NOT keyword_set(acc) then acc = 0.9

      results = 0
      npix = n_elements(spec.loglam)
      xvec = findgen(npix)

      nleft = n_elements(this_list)
      work_list = this_list[sort(this_list.x)]
      linterp, xvec, double(spec.loglam), work_list.x, work_wave

;
;     Locate MgII 2796/2803 lines in line list
;
      mg1 = (where(fix(restwave) EQ 2796))[0]
      if mg1 EQ -1 then return, 0

      mg2 = (where(fix(restwave) EQ 2803))[0]
      if mg2 EQ -1 then return, 0

;
;  Now run through work_list and identify all MgII doublets
;   Leave unidentified lines is the list
;
      j = 0

      while j LE nleft - 2 do begin      
       
        indices_left = lindgen(nleft)
        mgq = where(abs(work_list[j].x - work_list.x + pixel_sep[mg2]) $
                          LT acc AND indices_left GT j, nq) 
      
        if nq EQ 0 then j = j + 1
        if nq GT 1 then message, 'Multiple hits' 
        if nq EQ 1 then begin
;
;   Don't blindly search for MgII above 7800 Ang.
;
        good = 1
        if work_wave[mgq[0]] GT alog10(7800.0) then begin
           print, 'Discarding because wavelength falls in sky lines'
           good = 0
        endif
        
;
;   Reject odd EW ratios
;
        ratio = work_list[j].y/work_list[mgq].y 
        if ratio GT 2.2 OR ratio LT 0.8 then begin
          print, 'Discarding because EW ratio= ', ratio
          good = 0
        endif

;
;   Reject zabs > zqso + 0.1
;
        if work_wave[j] GT alog10((1.1+spec.zans.z)*restwave[mg1]) then begin
           print, 'Discarding because zabs > zqso + 0.1'
           good=0
        endif


        if good EQ 0 then j=j+1  $
        else begin
          absstruct1 = create_absstruct(work_list[[j,mgq]], [index,index], $
                     [spec.zans.z, spec.zans.z], work_wave[[j,mgq]], $
                     restwave[[mg1,mg2]], f_value[[mg1,mg2]], label[[mg1,mg2]])
          print, 'found MgII    ', index, [absstruct1.wave], $
               absstruct1[0].zabs, format='(a,i5, 2f10.4, f9.6)'

          if absstruct1[0].zabs GT spec.zans.z then $
               print, 'Warning zabs is larger than z_qso'
               
          if NOT keyword_set(results) then results = absstruct1 $
          else results = struct_append(results, absstruct1)         

          now = where(indices_left NE j AND indices_left NE mgq[0],nleft)
          if nleft EQ 0 then begin
            work_list = 0
            work_wave = 0
          endif else begin
            work_list = work_list[now]
            work_wave = work_wave[now]
          endelse
        endelse
      endif 
     
     endwhile 

     if NOT keyword_set(results) then return, 0

;
;   Now loop through the remaining significant lines and look for match
;   in the rest of line list
;

     mg1_found = where(fix(results.restwave) EQ 2796, nmgii)

;
;	Return if no MgII has been found
;
     if nmgii EQ 0 then return, results

;
;	Else let's look for line with known separations from 2796.352
;
     for j=0,nleft-1 do begin

        for ii = 0, nmgii - 1 do begin

          matches = where(abs(work_list[j].x - pixel_sep - $
            results[mg1_found[ii]].pixel) LT acc, nmatch)
 
          if nmatch EQ 0 then continue
          if nmatch GT 1 then message, 'Multiple Hits'

          matches = matches[0]
          absstruct1 = create_absstruct(work_list[j], index, $
                           spec.zans.z, work_wave[j], $
                         restwave[matches], f_value[matches], label[matches])
          print, 'found ', absstruct1.label, index, absstruct1.wave, $
                 absstruct1[0].zabs, format='(a,a, i5, f10.4, 10x, f9.6)'


          results = struct_append(results, absstruct1)         

        endfor
     endfor
          

     return, results
end
