function config_to_string, config

  if n_elements(config) gt 1 then begin
    sconfig = strarr(n_elements(config))
    for i=0, n_elements(config)-1 do sconfig[i] = config_to_string(config[i])
    return, sconfig
  endif

  if config lt 0 then config = 0
  if config lt 1000000 then $
    return, strtrim(string(config,f='(i6.6)'),2) $
  else $
    return, strtrim(config,2)

end
