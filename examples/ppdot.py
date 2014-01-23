#!/usr/bin/python

import sys
import cPickle

import matplotlib.pyplot as plt

# open file, read in model
filename = sys.argv[1]
f = open(filename, 'rb')
pop = cPickle.load(f)
f.close()

# lists to store p/pdot
periods = [pulsar.period for pulsar in pop.population]
pdots = [pulsar.pdot for pulsar in pop.population]

# plot a scatter log-log plot of the p/pdot values
plt.loglog(periods, pdots)
plt.show()