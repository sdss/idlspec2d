function config_to_long, config;, subexp=subexp
  if n_elements(config) gt 1 then begin
    lconfig = lonarr(n_elements(config))
    for i=0, n_elements(config)-1 do lconfig[i] = config_to_long(config[i],subexp=subexp)
    return, lconfig
  endif
  ;conid =strsplit(config,'-',/extract)
  ;if keyword_set(subexp) then $
    ;return, long(conid[1]) $
  ;else $
  ;  return, long(conid[0])
  return, long(config)
end
