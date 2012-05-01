#!/usr/bin/python

import sys
import argparse
import math
import random

import cPickle

import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector, Button, RadioButtons, CheckButtons

import numpy as np

from population import Population 
from pulsar import Pulsar


class Visualize:
    """
    Class for plotting different parameters of a population file
    using the matplotlib module
    """
    def __init__(self, popfile='populate.model'):
        # read in population file to population self.pop
        f = open(popfile, 'rb')
        self.pop = cPickle.load(f)
        f.close()

        self.textlabels = ['Period',
                      'DM',
                      'Gal X',
                      'Gal Y',
                      'Gal Z',
                      'W',
                      'SI',
                      'S1400',
                      'gl',
                      'gb',
                      'D']
        self.axislabels = ['Period (s)',
                      'DM (cm^-3 pc)',
                      'X (kpc)',
                      'Y (kpc)',
                      'Z (kpc)',
                      'Width (degrees)',
                      'Spectral Index',
                      'S1400 (mJy)',
                      'Galactic Longitude (degrees)',
                      'Galactic Latitude (degrees)',
                      'Distance (kpc)']

        if len(self.textlabels) != len(self.axislabels):
            print "Label list lengths not identical."
            sys.exit()

        # for each of those labels, get the data from population into lists
        # I suppose most sense is to read it into arrays
        
        # create numpy array of right size
        dataArray = np.zeros((len(self.textlabels), self.pop.size()), float)
        # loop over pulsars and fill array
        npsr = 0
        for psr in self.pop.population:
            dataArray[0][npsr] = psr.period
            dataArray[1][npsr] = psr.dm
            dataArray[2][npsr] = psr.galCoords[0]
            dataArray[3][npsr] = psr.galCoords[1]
            dataArray[4][npsr] = psr.galCoords[2]
            dataArray[5][npsr] = psr.width_degree
            dataArray[6][npsr] = psr.spindex
            dataArray[7][npsr] = psr.s_1400()
            dataArray[8][npsr] = psr.gl
            dataArray[9][npsr] = psr.gb
            dataArray[10][npsr] = psr.dtrue     
            npsr+=1

        self.dataArray = dataArray
        # delete population object, we don't need it now (I think!)
        # TBH not sure if this does anything, what with garbage collection
        del self.pop

    def display(self):
        """Method to create the plotting window and fill it with buttons"""
        # Create matplotlib window dressing
        self.fig = plt.figure(figsize=(11,9), facecolor='lightgoldenrodyellow')
        #self.fig = plt.figure(facecolor='lightgoldenrodyellow')
        self.fig.canvas.set_window_title('PyPop: Visualization')
        self.fig.clear()

        # put some buttons on the figure for the plot axes
        
        buttonW = 0.055
        buttonH = 0.045 # default button sizes

        # arguments are button left, bottom, width, height and label

        # check box to do log scale
        buttonXLog = CheckButtons(plt.axes([0.001, .3, 1.5*buttonW, 1.5*buttonH]),
                                    ['log x'],
                                    actives=[False])
        self.xlog = False

        buttonYLog = CheckButtons(plt.axes([0.001+2.5*buttonW, .3, 1.5*buttonW, 1.5*buttonH]),
                                    ['log y'],
                                    actives=[False])
        self.ylog = False


        # button to do plot
        buttonPlot = Button(plt.axes([0.001+2.0*buttonW, 0.2, 1.5*buttonW, buttonH]), 'Plot')
        
        self.fig.text(0.003, 0.98, 'PLOT SELECTION')
        self.fig.text(0.003, 0.95, 'X Axis')
        self.fig.text(0.003, 0.65, 'Y Axis')
            
        # radio buttons for x and y axis values
        radioXAxis = RadioButtons(plt.axes([0.001,
                                            0.72,
                                            len(self.textlabels)*0.02,
                                            len(self.textlabels)*0.02]), 
                                      self.textlabels,
                                      active=0,
                                      activecolor='blue')
        self.xindex = 0

        radioYAxis = RadioButtons(plt.axes([0.001,
                                            0.42,
                                            len(self.textlabels)*0.02,
                                            len(self.textlabels)*0.02]), 
                                      self.textlabels,
                                      active=1,
                                      activecolor='blue')
        self.yindex = 1

        # if radio button is clicked, change active X/Y index
        radioXAxis.on_clicked(self._radioXClick)
        radioYAxis.on_clicked(self._radioYClick)

        # add callback to the plot button(s?)
        buttonPlot.on_clicked(self._scatterplot)

        # callbacks for the log plot switches
        buttonXLog.on_clicked(self._logClick)
        buttonYLog.on_clicked(self._logClick)

        plt.show()

    def _logClick(self, label):
        """Function to switch logarithm booleans"""
        if label == 'log x':
            self.xlog = not self.xlog
        elif label == 'log y':
            self.ylog = not self.ylog
        else:
            print "Weird log label!"
            sys.exit()

    def _radioXClick(self, label):
        self.xindex = self.textlabels.index(label)

    def _radioYClick(self, label):
        self.yindex = self.textlabels.index(label)

    def _scatterplot(self, event):
        # if there's already a scatter plot, delete it!
        try:
            self.fig.delaxes(self.scatterplt)
        except:
            pass      

        # define axis position and dimensions
        self.scatterplt = self.fig.add_axes([0.3, 0.15, 0.65, 0.8])

        # axis labels
        self.scatterplt.set_xlabel(self.axislabels[self.xindex])
        self.scatterplt.set_ylabel(self.axislabels[self.yindex])

        # plot stuff. 
        self.scatterplt.plot(self.dataArray[self.xindex],
                             self.dataArray[self.yindex],
                             'r.')

        # if log switches are on, do a log plot
        if self.xlog:
            try:
                self.scatterplt.set_xscale('log')
            except ValueError:
                print "matplotlib refuses to do a log plot of the X data!"
                pass
        if self.ylog:
            try:
                self.scatterplt.set_yscale('log')
            except ValueError:
                print "matplotlib refuses to do a log plot of the Y data!"
                pass

        # finally, display the plot
        plt.show()

if __name__ == '__main__':
    """ 'Main' function for calling from command line"""
    parser = argparse.ArgumentParser(description='Visualize a population object')
    parser.add_argument('-f', metavar='fname', default='populate.model',
                          help='file containing population model (def="populate.model")')

    args = parser.parse_args()

    v = Visualize(popfile = args.f)

    v.display()
