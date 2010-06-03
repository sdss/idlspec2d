This directory initiated by A. Bolton on 15 May 2010.

Initial contents are a quick-and-dirty attempt to get BOSS templates
up and running, and include usage of both the traditional pca_solve
and a newly implemented "bolton_mlpca" (which is checked in under
pro/spec1d).

These routines have generated the following files, checked into the
traditional (one level up) templates dir:

spMLpcaCVstar-55332.fits
spMLpcaGal-55332.fits
spMLpcaQSO-55332.fits
spMLpcaStar-55332.fits

These conform to the spEigen* convention, but are named differently
to avoid confusion.  They will need to be renamed (or the spreduce1d
filename expectations changed) in order to be fed in.

The QSO templates were generated (most recently) using pca_solve,
and the rest with bolton_mlpca.

