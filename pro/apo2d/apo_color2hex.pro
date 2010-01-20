;------------------------------------------------------------------------------
function apo_color2hex, colorname

   case strupcase(strtrim(colorname,2)) of
   'RED': hexname = '#FF0000'
   'YELLOW': hexname = '#FFFF00'
   else: hexname = 'black'
   endcase

   return, hexname
end

;------------------------------------------------------------------------------
