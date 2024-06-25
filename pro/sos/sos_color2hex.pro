;------------------------------------------------------------------------------
function apo_color2hex, colorname

   case strupcase(strtrim(colorname,2)) of
   'RED': hexname = '#FF0000'
   'YELLOW': hexname = '#909000'
   else: hexname = 'black'
   else: hexname = 'black'
   endcase

   return, hexname
end

;------------------------------------------------------------------------------
