#!/usr/bin/python

import sys

from pulsar import Pulsar
from population import Population

def makepulsar(line):
    a = line.split()
    p = Pulsar(gl = float(a[0]),
               gb = float(a[1]),
               period = float(a[2]),
               width_degree = float(a[3])*360.0/float(a[2]),
               dm = float(a[4]),
               dtrue = float(a[7]),
               lum_1400 = float(a[9]),
               spindex = float(a[11]),
               scindex = -3.86
               )
    return p

pop = Population()
# open file
f = open(sys.argv[1], 'r')

# make a population object from the ascii file
for line in f.readlines():
    p = makepulsar(line.strip())
    pop.population.append(p)

f.close()

# write the object to file "pypop.pop"
pop.write("pypop.pop")
