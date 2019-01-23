; simulate the RMS variability of a continuum plus a line

pro sim_rms_var, nepoch=nepoch, result = result

if ~keyword_set(nepoch) then nepoch=1000L


; the constant continuum
cont0 = 1.0
; the variable continuum
cfrac = 0.1
cont1 = randomn(101L, nepoch)*cfrac*cont0


; the constant line
ew = 0.5
line0 = cont0*ew
; the variable line
lfrac = 0.2
line1 = randomn(101L, nepoch)*lfrac*line0

cont = cont0 + cont1
line = line0 + line1
tot = cont + line

result={cont:cont, line:line, tot:tot}

end
