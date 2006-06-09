
function gauss_l2, x, p, lvls=lvls

     xcen = p[1]
     sig_g = p[2]
     sig_e = p[3]

     d = x-xcen

     mind=min(d)
     maxd=max(d)
     nx = long(((maxd-mind)+10)/0.02) +1
     checkx = (findgen(nx) - 250)*0.02 + mind 
     
     model = p[0] * 2.5 / (1. + abs(checkx/sig_e))^6
 
     g = gauss_kernel(sig_g * 50.)
     c = convol(model, g, /edge_t)

     linterp, checkx, alog10(c), d, f
     print, total(10^f)
     frac = ((1- errorf([1.0, 2.0, 3.0, 4.0, 5.0d]/sqrt(2.0)))/2.)
     cumul = total(c, /cumul)/total(c)
     linterp, alog(cumul), checkx, alog(frac), lvls
 
;     print, lvls+xcen, xcen-lvls

    
;     plot, x, 10^f, /ylog
     return, f
end

