pro bpm_ledflat_combine2

imb1=mrdfits('pixflatave-55055-b1.old.fits.gz',0,hdrb1)
imb2=mrdfits('pixflatave-55055-b2.old.fits.gz',0,hdrb2)
imr1=mrdfits('pixflatave-55055-r1.old.fits.gz',0,hdrr1)
imr2=mrdfits('pixflatave-55055-r2.old.fits.gz',0,hdrr2)
imr2new=mrdfits('pixflatave-55300-r2.old.fits.gz',0,hdrnewr2)

imbpmb1=mrdfits('bpm-b1-count.fit.gz')
imbpmb2=mrdfits('bpm-b2-count.fit.gz')
imbpmr1=mrdfits('bpm-r1-count.fit.gz')
imbpmr2=mrdfits('bpm-r2-count.fit.gz')
imbpmr2new=mrdfits('bpm-r2-55300-count.fit.gz')


imoutb1=imb1*imbpmb1
imoutb2=imb2*imbpmb2
imoutr1=imr1*imbpmr1
imoutr2=imr2*imbpmr2
imoutr2new=imr2new*imbpmr2new

;for b1                         ;bad columns by eye


;for b2

imoutb2[334,1335:2057]=0
imoutb2[1966,1636:2057]=0
imoutb2[2493,2056:3411]=0
imoutb2[1825:1826,0:2057]=0     ;charge transfer

;for r1
imoutr1[2181:2182,0:2064]=0
imoutr1[3648:3650,0:2064]=0
imoutr1[693,2064:4127]=0
imoutr1[781,2064:4127]=0

;for r2

imoutr2[2111,1800:2065]=0  
imoutr2[551,2060:4127]=0
imoutr2[1113,2060:2505]=0
;for new r2

imoutr2new[508,2064:4127]=0
imoutr2new[2627:2628,2064:4127]=0
imoutr2new[2671:2672,2064:4127]=0
imoutr2new[2924:2926,2064:3023]=0
imoutr2new[2984:2985,2064:4127]=0
imoutr2new[3882,2064:4127]=0
imoutr2new[3955:3956,2064:4127]=0
imoutr2new[3872:3882,2612:3664]=0
imoutr2new[3947:3955,2845:3108]=0


writefits,'pixflatave-55055-b1.fits',imoutb1,hdrb1
writefits,'pixflatave-55055-b2.fits',imoutb2,hdrb2
writefits,'pixflatave-55055-r1.fits',imoutr1,hdrr1
writefits,'pixflatave-55055-r2.fits',imoutr2,hdrr2
writefits,'pixflatave-55300-r2.fits',imoutr2new,hdrnewr2

end
