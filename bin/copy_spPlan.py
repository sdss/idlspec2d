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
parser.add_option("--run2d", type="string",  help="output RUN2D version")
parser.add_option("-n", "--numplates", type="int",  help="number of plates to copy [default all]")
parser.add_option("-R", "--randseed", type="int", default=0, help="random seed [default 0]")
### parser.add_option("--run1d", type="string",  help="output RUN1D version")
### parser.add_option("-x", "--xxx",   help="some flag", action="store_true")
opts, args = parser.parse_args()

#- Set random seed so that results are reproducible
random.seed(opts.randseed)

#- Default input directory $BOSS_SPECTRO_REDUX/$RUN2D/
if opts.input is None:
    opts.input = os.environ['BOSS_SPECTRO_REDUX'] + "/" + os.environ['RUN2D']

#- required options
if opts.run2d is None or opts.output is None:
    print >> sys.stderr, 'ERROR: you must specify -o/--output directory and --run2d options'
    print >> sys.stderr, 'To see all options, run copy_spPlan.py -h'
    sys.exit(1)
    
#- Create output directory if needed
if not os.path.isdir(opts.output):
    os.makedirs(opts.output)

#- Get list of plates in input directory
inplates = glob(opts.input + '/[0-9][0-9][0-9][0-9]')

#- Downsample if requested
if opts.numplates is not None and opts.numplates < len(inplates):
    inplates = random.sample(inplates, opts.numplates)

#- Loop over plates, copying the plan files
for platedir in sorted(inplates):
    plate = os.path.basename(platedir)
    print '\rPlate ' + plate,
    sys.stdout.flush()
    plan2dfiles = glob(platedir + '/spPlan*.par')
    for planfile in plan2dfiles:
        outdir = opts.output + "/" + plate
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        outplan = outdir + '/' + os.path.basename(planfile)
        copyplan(planfile, outplan, opts.run2d)

#- final blank line print to get CR since we were being fancy with '\r...'
print