; dummy procedure to return a bunch of zeros
function sdss_astrombad, xx,yy,zz

nnn = n_elements(xx)

return, replicate(0, nnn)

end
