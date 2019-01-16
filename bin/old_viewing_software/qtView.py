#!/usr/bin/env python

import sys
import os
import argparse
import random

import cPickle
import numpy as np

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib
from matplotlib.figure import Figure

from matplotlib.backends.backend_qt4agg import \
    FigureCanvasQTAgg as FigCanvas, \
    NavigationToolbar2QTAgg as NavigationToolbar

class ViewException(Exception):
    pass

class DataObj:
    def __init__(self, name, textlabels, axislabels, dataarray, size):
        self.name = name
        self.textlabels = textlabels
        self.axislabels = axislabels
        self.dataArray = dataarray
        self.count = size

def makeDataObj(frac=None, extnList=['.model', '.results']):
    """
    Function to get different parameters of a population file

    """
    textlabels = ['Period',
                    'DM',
                    'Gal X',
                    'Gal Y',
                    'Gal Z',
                    'W',
                    'alpha',
                    'rho',
                    'SI',
                    'S1400',
                    'gl',
                    'gb',
                    'D',
                    'r0',
                    'n']
    axislabels = ['Period (ms)',
                    'DM (cm^-3 pc)',
                    'X (kpc)',
                    'Y (kpc)',
                    'Z (kpc)',
                    'Width (degrees)',
                    'alpha (deg)',
                    'rho (deg)',
                    'Spectral Index',
                    'S1400 (mJy)',
                    'Galactic Longitude (degrees)',
                    'Galactic Latitude (degrees)',
                    'Distance (kpc)',
                    'GalacticRadius (kpc)',
                    'Array Index']

    if len(textlabels) != len(axislabels):
        print "Label list lengths not identical."
        sys.exit()

    dataObjList = []
    
    for filename in os.listdir(os.getcwd()):
        for extn in extnList:
            if filename.endswith(extn):
                # read in population file to population self.pop
                try:
                    f = open(filename, 'rb')
                except IOError:
                    print "Could not open file {0}.".format(filename)
                    sys.exit()

                pop = cPickle.load(f)
                f.close()

                # create numpy array of right size
                dataArray = np.zeros((len(textlabels), pop.size()), float)
                # loop over pulsars and fill array
                npsr = 0

                # factor to reduce the number of plotted pulsars
                # to improve speed
                # make this an option
                if frac is None or frac > 1.0:
                    frac = 1.0

                for psr in pop.population:
                    if random.random() < frac:
                        dataArray[0][npsr] = psr.period
                        dataArray[1][npsr] = psr.dm
                        dataArray[2][npsr] = psr.galCoords[0]
                        dataArray[3][npsr] = psr.galCoords[1]
                        dataArray[4][npsr] = psr.galCoords[2]
                        dataArray[5][npsr] = psr.width_degree
                        dataArray[6][npsr] = psr.alpha
                        dataArray[7][npsr] = psr.rho
                        dataArray[8][npsr] = psr.spindex
                        dataArray[9][npsr] = psr.s_1400()
                        dataArray[10][npsr] = psr.gl
                        dataArray[11][npsr] = psr.gb
                        dataArray[12][npsr] = psr.dtrue 
                        dataArray[13][npsr] = psr.r0
                        dataArray[14][npsr] = npsr
                        npsr+=1

                dataObjList.append(
                         DataObj(filename, textlabels, axislabels, dataArray, pop.size()))

    if len(dataObjList) == 0:
        raise ViewException('No files found matching the extensions')

    # sort list from largest -> smallest so plotting is always done that way
    dataObjList.sort(key=lambda x: x.count, reverse=True)

    return dataObjList


class VisualizeFrame(QMainWindow):
    def __init__(self, dataObjList, parent=None):

        # data object
        self.dataObjList = dataObjList

        # initialise the frame
        QMainWindow.__init__(self, parent)

        self.colour_list = ['b.', 'r.', 'g.', 'c.', 'm.', 'y.', 'k.']
        self.create_main_panel()
        self.draw_figure()


    def create_main_panel(self):
        self.main_frame = QWidget()

        self.dpi = 100
        self.fig = Figure((5., 5.), dpi = self.dpi)
        #self.fig.set_figsize_inches( ( 5., 5.0))
        self.canvas = FigCanvas(self.fig)
        self.canvas.setParent(self.main_frame)

        self.axes = self.fig.add_subplot(111)

        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

        #  The simple buttons
        #
        self.drawbutton = QPushButton("&Plot")
        self.connect(self.drawbutton, SIGNAL('clicked()'), self.draw_figure)

        self.logx = QCheckBox("&log X")
        self.logx.setChecked(False)
        self.connect(self.logx, SIGNAL('stateChanged(int)'), self.draw_figure)

        self.logy = QCheckBox("&log Y")
        self.logy.setChecked(False)
        self.connect(self.logy, SIGNAL('stateChanged(int)'), self.draw_figure)

        self.grid = QCheckBox("&Grid")
        self.grid.setChecked(False)
        self.connect(self.grid, SIGNAL('stateChanged(int)'), self.draw_figure)
        
        # make a radio button list
        self.XRadioButtonList=[]
        self.YRadioButtonList=[]
        
        self.RadioGrid = QGridLayout()

        self.XGB = QGroupBox('&X Axis')
        self.YGB = QGroupBox('&Y Axis')

        for i, l in enumerate(self.dataObjList[0].textlabels):
            item1 = QRadioButton(l)
            self.XRadioButtonList.append(item1)
            item2 = QRadioButton(l)
            self.YRadioButtonList.append(item2)
            self.connect(self.XRadioButtonList[len(self.XRadioButtonList)-1],
                            SIGNAL('stateChanged(int)'),
                            self.draw_figure)
            self.connect(self.YRadioButtonList[len(self.YRadioButtonList)-1],
                            SIGNAL('stateChanged(int)'),
                            self.draw_figure)

            if i==0:
                item1.setChecked(True)
                item2.setChecked(True)

        ## Layout for the radio button groups
        ##
        vboxX = QVBoxLayout()
        for w in self.XRadioButtonList:
            vboxX.addWidget(w)
        vboxY = QVBoxLayout()
        for w in self.YRadioButtonList:
            vboxY.addWidget(w)

        vboxX.addStretch(1)
        vboxY.addStretch(1)

        self.XGB.setLayout(vboxX)
        self.YGB.setLayout(vboxY)
        self.RadioGrid.addWidget(self.XGB, 0, 0)
        self.RadioGrid.addWidget(self.YGB, 0, 1)
        
        # make a check box list. A bit complicated!
        self.modelList = [d.name for d in self.dataObjList]
        self.modelCheckBoxList=[]
        for model in self.modelList:
            self.modelCheckBoxList.append(QCheckBox(model))
            self.modelCheckBoxList[len(self.modelCheckBoxList)-1].setChecked(False)
            self.connect(self.modelCheckBoxList[len(self.modelCheckBoxList)-1],
                            SIGNAL('stateChanged(int)'),
                            self.draw_figure)

        # position the layouts & buttons in the application
        vboxchecklist = QVBoxLayout()
        for w in self.modelCheckBoxList:
            vboxchecklist.addWidget(w)

        # line up the top row
        vboxchecklist.addStretch(1)
        hbox1 = QHBoxLayout()
        hbox1.addLayout(self.RadioGrid)
        hbox1.addWidget(self.canvas)
        hbox1.addLayout(vboxchecklist)

        # line up the buttons at the bottom
        hbox2 = QHBoxLayout()
        for w in [self.drawbutton,
                  self.logx,
                  self.logy,
                  self.grid]:
            hbox2.addWidget(w)
            hbox2.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addLayout(hbox1)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(hbox2)

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def _getCheckedRBIndex(self, buttonlist):
        for i, b in enumerate(buttonlist):
            if b.isChecked():
                return i

        # if we got here, something's broken - radio buttons 
        # should always have one checked! 
        raise ViewException('No radio button selected for plotting(!) ')


    def draw_figure(self):
        self.axes.clear()

        # get the x and y indices of selected radio buttons
        xIndex = self._getCheckedRBIndex(self.XRadioButtonList)
        yIndex = self._getCheckedRBIndex(self.YRadioButtonList)

        for dataObjIndex, val in enumerate(self.modelList):
            if self.modelCheckBoxList[dataObjIndex].isChecked():
                try:
                    self.axes.plot(self.dataObjList[dataObjIndex].dataArray[xIndex],
                                     self.dataObjList[dataObjIndex].dataArray[yIndex],
                                     self.colour_list[dataObjIndex], 
                                     label=self.dataObjList[dataObjIndex].name)
                except IndexError:
                    self.axes.plot(self.dataObjList[dataObjIndex].dataArray[xIndex],
                                   self.dataObjList[dataObjIndex].dataArray[yIndex],
                                   label=self.dataObjList[dataObjIndex].name)

        self.axes.set_xlabel(self.dataObjList[0].axislabels[xIndex], 
                                fontsize=10)
        self.axes.set_ylabel(self.dataObjList[0].axislabels[yIndex],
                                fontsize=10)

        for label in self.axes.get_xticklabels():
            label.set_fontsize(7)
        for label in self.axes.get_yticklabels():
            label.set_fontsize(7)

        if self.logx.isChecked():
            self.axes.set_xscale('log')
            for label in self.axes.get_xticklabels():
                label.set_fontsize(7)

        if self.logy.isChecked():
            self.axes.set_yscale('log')
            for label in self.axes.get_yticklabels():
                label.set_fontsize(7)

        self.axes.grid(self.grid.isChecked())
        self.canvas.draw()

if __name__ == '__main__':
    """ 'Main' function for calling from command line"""
    parser = argparse.ArgumentParser(description='Visualize a population object')
    parser.add_argument('-f', metavar='fname', default='populate.model',
                          help='file containing population model (def="populate.model")')

    parser.add_argument('-frac', nargs=1, type=float, default=None, 
                          help='plot only this fraction of pulsars')

    parser.add_argument('-extn', nargs='+', type=str,
                        default=['.results', '.model'],
                        help='extension(s) to look for when finding population models')
    args = parser.parse_args()

    dataObj = makeDataObj(frac=args.frac, extnList=args.extn)

    app = QApplication(sys.argv)
    frame = VisualizeFrame(dataObj)
    frame.show()
    app.exec_()

