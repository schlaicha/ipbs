#!/usr/bin/python

import sys

from collections import defaultdict

def nesteddict(): 
  return defaultdict(nesteddict)

if len (sys.argv) != 2:
  print "This prints the results of parallel iPBS into single files."
  print "How to use: %s # procs" % sys.argv[0]
  sys.exit()

nprocs = int(sys.argv[1])

# create a dictionary
entries = nesteddict()

# copy together the reference stuff
for i in range(0,nprocs):
  refFileName = "s%.4d:p%.4d:reference.dat" %(nprocs,i)
  refFile = open(refFileName,"r")
  solFileName = "s%.4d:p%.4d:ipbs_solution.dat" %(nprocs,i)
  solFile = open(solFileName,"r")
  refFile.readline()
  for line in refFile:
    s = line.split()
    n = len(s)
    if n == 0:
      break
    key = s[0],s[1]
    value = float(s[2])
    entries[key]['reference'] = value
  solFile.readline()
  for line in solFile:
    s = line.split()
    n = len(s)
    key = s[0],s[1]
    entries[key]['solution'] = float(s[2])

for key in entries:
  entries[key]['difference'] = (entries[key]['solution']) - (entries[key]['reference'])
  if ((entries[key]['solution'] + entries[key]['reference'])) != 0:
    entries[key]['relative'] = 100*((entries[key]['solution'] - entries[key]['reference']) 
        / (entries[key]['solution'] + entries[key]['reference']))
  else:
    entries[key]['relative'] = 0;


outfile = open("data.dat","w")
for key in entries:
  outfile.writelines("%s %s %s %s %s %s \n" %(key[0], key[1] ,entries[key]['reference'],
    entries[key]['solution'], entries[key]['difference'], entries[key]['relative']))
