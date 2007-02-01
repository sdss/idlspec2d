function spflux_masklines, loglam, hwidth=hwidth, stellar=stellar, $
                           telluric=telluric, o2=o2, h2o=h2o, $
                           extralines=extralines, extrabands=extrabands, $
                           spans=spans
    
    if keyword_set(telluric) then begin
       o2=1
       h2o=1
    end

    o2_bands = [[6860, 7575], $
                [6970, 7710]]
    h2o_bands = [[5145, 5905, 6270, 6705, 6970, 7150., 8090., 8930.], $
                 [5202, 5950, 6290, 6725, 7075, 7350., 8360., 9250.]]
    ;if (keyword_set(telluric)) then begin
    ;   tellwave1 = [5150, 5900, 6270, 6705, 6850., 7150., 7560., 8105., 8930.]
    ;   tellwave2 = [5200, 5950, 6290, 6725, 6960., 7350., 7720., 8360., 9160.]
    ;endif

    if (NOT keyword_set(hwidth)) then $
       hwidth = 5.7e-4          ; Default is to mask +/- 5.7 pix = 400 km/sec

    mask = bytarr(size(loglam,/dimens))

    if (keyword_set(stellar)) then begin
       starwave = [ $
                  3830.0 , $    ; ? (H-7 is at 3835 Ang)
                  3889.0 , $    ; H-6
                  3933.7 , $    ; Ca_k
                  3968.5 , $    ; Ca_H (and H-5 at 3970. Ang)
                  4101.7 , $    ; H-delta
                  4300.  , $    ; G-band
                  4305.  , $    ; G-band
                  4310.  , $    ; more G-band
                  4340.5 , $    ; H-gamma
                  4861.3 , $    ; H-beta
                  5893.0 , $    ; Mg
                  6562.8 , $    ; H-alpha
                  8500.8 , $
                  8544.6 , $
                  8665.0 , $
                  8753.3 , $
                  8866.1 , $
                  9017.5 , $
                  9232.0 ]
    endif
    if keyword_set(extralines) then begin
       if keyword_set(starwave) then $
          starwave = [starwave, extralines] $
       else $
          starwave = extralines
    end
    
    if keyword_set(starwave) then begin
       airtovac, starwave
       
       for i=0L, n_elements(starwave)-1 do begin
          wmin = alog10(starwave[i])-hwidth
          wmax = alog10(starwave[i])+hwidth
          mask = mask OR (loglam GT wmin $
                          AND loglam LT wmax)

          if arg_present(spans) then begin
             if keyword_set(spanmin) then begin
                spanmin = [spanmin, wmin]
                spanmax = [spanmax, wmax]
             end else begin
                spanmin = [wmin]
                spanmax = [wmax]
             end
          endif
       endfor
    endif
    
    if (keyword_set(o2)) then begin
       tellwave1 = o2_bands[*,0]
       tellwave2 = o2_bands[*,1]
    end
    if (keyword_set(h2o)) then begin
       if keyword_set(tellwave1) then begin
          tellwave1 = [tellwave1, h2o_bands[*,0]]
          tellwave2 = [tellwave2, h2o_bands[*,1]]
       end else begin
          tellwave1 = h2o_bands[*,0]
          tellwave2 = h2o_bands[*,1]
       end
    end
    
    
    if keyword_set(extrabands) then begin
       if keyword_set(tellwave1) then begin
          tellwave1 = [tellwave1, extrabands[*,0]]
          tellwave2 = [tellwave2, extrabands[*,1]]
       end else begin
          tellwave1 = extrabands[*,0]
          tellwave2 = extrabands[*,1]
       end
       
    end
    
    if keyword_set(tellwave1) then begin
       for i=0L, n_elements(tellwave1)-1 do begin
          wmin = alog10(tellwave1[i])
          wmax = alog10(tellwave2[i])
          mask = mask OR (loglam GT wmin $
                          AND loglam LT wmax)

          if arg_present(spans) then begin
             if keyword_set(spanmin) then begin
                spanmin = [spanmin, wmin]
                spanmax = [spanmax, wmax]
             end else begin
                spanmin = [wmin]
                spanmax = [wmax]
             end
          endif
       endfor
    endif
    
    if arg_present(spans) then spans = [[spanmin],[spanmax]]
    
    return, mask
 end
