
function create_absstruct, peak, index, zqso, this_wave, $
        restwave, f_value, label 
;
;      Prepare structure template
;
      fill  = { label : '', index : 0L, $
                 pixel : 0.0, pixel_err : 0.0, $
                 wave  : 0.0d, wave_err : 0.0d, $
                 ew_pix  : 0.0, ew_r : 0.0, ew_sn : 0.0, $
                 sigma_pix  : 0.0, $
                 restwave : 0.0d, f_value : 0.0, $
                 zabs : 0.0, zabs_err : 0.0, z_qso : 0.0 }

      fill = replicate(fill,n_elements(peak))  

          pixel_width    = 1.0d-4 * alog(10)

          fill.index     = index
          fill.pixel     = peak.x
          fill.pixel_err = peak.xerr
          fill.wave      = 10^this_wave
          fill.wave_err  = pixel_width*peak.xerr*10^this_wave

          fill.ew_pix      = peak.y 
          fill.ew_sn     = peak.sn
          fill.ew_r      = pixel_width * restwave * fill.ew_pix
          fill.label     = label
          fill.restwave  = restwave
          fill.f_value   = f_value

          fill.zabs      = fill.wave / fill.restwave - 1.0
          fill.z_qso     = zqso

       return, fill
end





          
         
      
