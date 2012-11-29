#!/usr/bin/env python

"""
Utility script to copy the spPlan* files from one production to another
while updating the RUN2D entries appropriately.

Stephen Bailey, LBL
Fall 2012
"""

import sys
import os
import os.path
import random
from glob import glob
import pyfits

#- copy spPlan file while updating the RUN2D entry
def copyplan(inplan, outplan, run2d):
    finput = open(inplan)
    foutput = open(outplan, 'w')
    for line in finput:
        if line.startswith('RUN2D'):
            xx = line.split(None, 2)    #- RUN2D VER [# Comment]
            xx[1] = run2d               #- replace RUN2D
            line = " ".join(xx) + '\n'  #- put back together with newline
            
        foutput.write(line)
        
    finput.close()
    foutput.close()

#-------------------------------------------------------------------------
import optparse
parser = optparse.OptionParser(usage = "%prog [options]",
description="""Copy spPlan files from one redux version to another while replacing RUN2D.
""")
parser.add_option("-i", "--input", type="string",  help="input directory [default $BOSS_SPECTRO_REDUX/$RUN2D/]")
parser.add_option("-o", "--output", type="string",  help="output directory")
parser.add_option("--run2d",  type="string",  help="output RUN2D version")
parser.add_option("--minmjd", type="int",  help="min MJD to include", default=0)
parser.add_option("--maxmjd", type="int",  help="max MJD to include", default=100000)
parser.add_option("-n", "--numplates", type="int",  help="number of plates to copy [default all good ones]")
parser.add_option("-R", "--randseed", type="int", default=0, help="random seed [default 0]")
### parser.add_option("--run1d", type="string",  help="output RUN1D version")
parser.add_option("-b", "--bad",   help="also copy bad quality plans, not just good ones", action="store_true")
parser.add_option("-p", "--platelist", help="override platelist location [default input/platelist.fits]")
opts, args = parser.parse_args()

#- Set random seed so that results are reproducible
random.seed(opts.randseed)

#- Default input directory $BOSS_SPECTRO_REDUX/$RUN2D/
if opts.input is None:
    opts.input = os.environ['BOSS_SPECTRO_REDUX'] + "/" + os.environ['RUN2D']

#- required options
if opts.output is None:
    print >> sys.stderr, 'ERROR: you must specify -o/--output directory'
    print >> sys.stderr, 'To see all options, run copy_spPlan.py -h'
    sys.exit(1)
    
#- choose run2d based upon output name if needed
if opts.run2d is None:
    opts.run2d = os.path.basename(opts.output)
    if opts.run2d in (None, '', '.'):
        opts.run2d = os.path.basename(os.path.dirname(opts.output))
        if opts.run2d in (None, '', '.'):        
            print "ERROR: Unable to derive RUN2D from path", opts.output
            print "ERROR: use --run2d instead"
            sys.exit(2)
    
    print "Using RUN2D="+opts.run2d
            
#- Create output directory if needed
if not os.path.isdir(opts.output):
    os.makedirs(opts.output)

#- Find platelist file
if opts.platelist is None:
    opts.platelist = opts.input + '/platelist.fits'
if not os.path.exists(opts.platelist):
    print >> sys.stderr, "ERROR: if no platelist.fits in input dir, you must specify a platelist"
    sys.exit(1)

#- Create set of plates with at least one good plugging within the
#- desired MJD range
p = pyfits.getdata(opts.platelist, 1)
goodplates = set()
for plate, mjd, quality in zip(p['PLATE'], p['MJD'], p['PLATEQUALITY']):
    if (quality.strip() == 'good' or opts.bad):
        if opts.minmjd <= mjd <= opts.maxmjd:        
            goodplates.add( plate )

#- Randomly subsample
if opts.numplates is not None:
    goodplates = set( random.sample(goodplates, opts.numplates) )

#- Make matching list of good pluggings for those plates
goodplugs = set()
for plate, mjd, quality in zip(p['PLATE'], p['MJD'], p['PLATEQUALITY']):
    if plate in goodplates and (quality.strip() == 'good' or opts.bad):
        goodplugs.add( (plate, mjd) )

#- Loop over plates, copying the plan files
ncopied = 0
for plate in sorted(goodplates):
    platedir = os.path.join(opts.input, str(plate))
    print '\rPlate', plate,
    sys.stdout.flush()
    
    #- Copy all the plan2d files since they may be needed for the coadd
    #- even if they aren't in platelist
    plan2dfiles = glob(platedir + '/spPlan2d*.par')
    for planfile in plan2dfiles:
        outdir = opts.output + "/" + str(plate)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        outplan = outdir + '/' + os.path.basename(planfile)
        copyplan(planfile, outplan, opts.run2d)

    #- Copy only the good plate-mjd plancomb files
    plancombfiles = glob(platedir + '/spPlancomb*.par')
    for planfile in plancombfiles:
        mjd = int(os.path.basename(planfile).split('-')[2][0:5])
        if (plate, mjd) in goodplugs:
            outplan = outdir + '/' + os.path.basename(planfile)
            copyplan(planfile, outplan, opts.run2d)
                    
#- final blank line print to get CR since we were being fancy with '\r...'
print

