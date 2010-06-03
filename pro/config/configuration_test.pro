pro configuration_test

c1=obj_new('configuration')
c2=obj_new('configuration',56000)

print, '5ai: ',c1->getNumberFibersTotal(),c2->getNumberFibersTotal()
print, '5aii: ',c1->getNumberFibersPerSpectrograph(), c2->getNumberFibersPerSpectrograph()
print, '5aiii: ',c1->getNumberFibersPerBundle(), c1->getNumberBundles(),c2->getNumberFibersPerBundle(),c2->getNumberBundles()

print, '5bi: ', (c1->getDetectorFormat('blue'))[0], (c2->getDetectorFormat('red'))[0]
print, '5bii: ', (c1->getDetectorFormat('red'))[1], (c2->getDetectorFormat('blue'))[1]

xerr=[0,0.01,0.05,0.1]
print, 'xerr=',xerr
print, 'inmask: '
temp=c1->spcalib_inmask(xerr)
temp2=c2->spcalib_inmask(xerr)
print, temp, temp2

obj_destroy,c1
obj_destroy,c2
end
