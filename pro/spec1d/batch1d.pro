pro batch1d, localfullfile

   if (NOT keyword_set(localfullfile)) then $
    localfullfile = findfile('*/spPlate-*.fits')
   localfile = fileandpath(localfullfile, path=localpath)

   platemjd = strmid(localfile, 8, 10)
   zfile = localpath + '/spZ-' + platemjd + '.fits'
   diaglog = localpath + '/spDiag1d-' + platemjd + '.log'
   diagps = localpath + '/spDiag1d-' + platemjd + '.ps'
   outfile = transpose( [ [diaglog], [diagps], [zfile] ] )

   hostfile = filepath('batch1d.par', $
    root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='examples')
   yanny_read, hostfile, pp
   hostconfig = *pp[0]
   yanny_free, pp

   endstring = 'Successful completion of SPREDUCE1D'
   sq = "\'"
   command = 'echo "cd,'+sq+localpath+sq + ' \& spreduce1d,'+sq+localfile+sq $
    + '" | idl'
   localfullfile = reform([localfullfile], 1, n_elements(localfullfile))
   batch, localfullfile, outfile, $
    hostconfig.protocol, hostconfig.remotehost, hostconfig.remotedir, $
    command, endstring

   return
end

