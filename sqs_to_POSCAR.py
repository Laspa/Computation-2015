#!/usr/bin/env python

#==============================================================================
#  Jonas Kaufman jlkaufman@hmc.edu
#  July 20, 2015
#  Script to convert a bestsqs.out file from ATAT to a POSCAR file for use  
#  with VASP
#  Uses direct coordinates
#==============================================================================
"""
Make script executable using 'chmod +x _____.py' to call as bash script
Requires Cell.py
"""
from Cell import *

# get filenames
sqs = raw_input('SQS filename (input): ')
if not sqs: sqs = 'bestsqs.out'
print sqs,'\n'
pos = raw_input('POSCAR filename (output):')
if not pos: pos = 'POSCAR'
print pos,'\n'

# load SQS
print 'Converting SQS...'
cell = Cell().loadFromSQS(sqs)
print 'Done\n'

# send to POSCAR
print 'Sending to POSCAR...'
cell.sendToPOSCAR(pos)
print 'Done\n'