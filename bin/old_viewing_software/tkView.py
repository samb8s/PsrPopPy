#!/usr/bin/env python

import sys
import os
import argparse
import random

import cPickle
import numpy as np

import Tkinter as tk
import matplotlib
matplotlib.use('TKAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import \
    FigureCanvasTkAgg as FigCanvas, \
    NavigationToolbar2TkAgg as NavigationToolbar

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

                # for each of those labels, get the data from population into lists
                # I suppose most sense is to read it into arrays
                
                # create numpy array of right size
                dataArray = np.zeros((len(textlabels), pop.size()), float)
                # loop over pulsars and fill array
                npsr = 0

                # going to throw in a factor to reduce the number of plotted pulsars
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


class VisualizeFrame(tk.Frame):
    def __init__(self, dataObjList):

        self.dataObjList = dataObjList
        tk.Frame.__init__(self, None, -1)

        self.colour_list = ['b.', 'r.', 'g.', 'c.', 'm.', 'y.', 'k.']
        self.create_main_panel()



    def create_main_panel(self):
        self.panel = tk.Panel(self)
        self.dpi = 100
        self.fig = Figure((5., 5.), dpi = self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)

        self.axes = self.fig.add_subplot(111)

        self.drawbutton = tk.Button(self.panel, -1, "Plot")
        self.Bind(tk.EVT_BUTTON, self.on_draw_button, self.drawbutton)

        self.logx = tk.CheckBox(self.panel, -1, 
            "log X",
            style=tk.ALIGN_RIGHT)
        self.Bind(tk.EVT_CHECKBOX, self.on_logx, self.logx)

        self.logy = tk.CheckBox(self.panel, -1, 
            "log Y",
            style=tk.ALIGN_RIGHT)
        self.Bind(tk.EVT_CHECKBOX, self.on_logy, self.logy)

        self.grid = tk.CheckBox(self.panel, -1, 
            "grid",
            style=tk.ALIGN_RIGHT)
        self.Bind(tk.EVT_CHECKBOX, self.on_grid, self.grid)
        
        modelList = [d.name for d in self.dataObjList]

        # create list of population models. Set all "on" by default
        self.modelCheckList = tk.CheckListBox(self.panel, -1,
                                            choices=modelList,
                                            style=tk.ALIGN_RIGHT)
        self.modelCheckList.SetChecked(range(len(modelList)))

        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas)

        #
        # Layout with box sizers
        #

        self.vbox = tk.BoxSizer(wx.VERTICAL)
        self.v_buttonbox_1 = tk.BoxSizer(wx.VERTICAL)
        self.v_buttonbox_2 = tk.BoxSizer(wx.VERTICAL)
        self.hpltbox = tk.BoxSizer(wx.HORIZONTAL) # for plot and plot selection
        self.htoolbox = tk.BoxSizer(wx.HORIZONTAL) # for drawing and log toggles

        # fill the top box first, radio buttons in a vbox, 
        # next to the canvas
        self.radioBoxX = tk.RadioBox(self.panel, 1, 'X Axis',
                                    choices = self.dataObjList[0].textlabels,
                                    majorDimension=1,
                                    style = tk.RA_SPECIFY_COLS)

        self.radioBoxY = tk.RadioBox(self.panel, 2, 'Y Axis',
                                    choices = self.dataObjList[0].textlabels,
                                    majorDimension=1,
                                    style = tk.RA_SPECIFY_COLS)

        self.xIndex=0
        self.yIndex=0

        # event for radio box with ID 1
        tk.EVT_RADIOBOX(self.panel, 1, self.onXRadioClick)
        # event for radiobox with ID 2
        tk.EVT_RADIOBOX(self.panel, 2, self.onYRadioClick)


        # top horizontal panel - canvas and radio boxes
        self.hpltbox.Add(self.radioBoxX)
        self.hpltbox.Add(self.radioBoxY)
        self.hpltbox.Add(self.canvas, 1, tk.LEFT | wx.TOP | wx.GROW)
        self.hpltbox.Add(self.modelCheckList, 0, border=3)

        # add the matplotlib toolbar
        self.vbox.Add(self.hpltbox)
        self.vbox.Add(self.toolbar, 0, tk.EXPAND)
        self.vbox.AddSpacer(10)


        # bottom horizontal panel - check buttons and plot button
        flags = tk.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.htoolbox.Add(self.drawbutton, 0, border=3, flag=flags)
        self.htoolbox.Add(self.logx, 0, border=3, flag=flags)
        self.htoolbox.Add(self.logy, 0, border=3, flag=flags)
        self.htoolbox.AddSpacer(20)
        self.htoolbox.Add(self.grid, 0, border=3, flag=flags)

        self.vbox.Add(self.htoolbox, 0, flag = tk.ALIGN_LEFT | wx.TOP)

        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)

    def onRadioClick(self, event):
        radioBox = event.GetEventObject()
        print event.GetId()
        print radioBox

    def onXRadioClick(self, event):
        # get x index of selected button
        radioBox = event.GetEventObject()
        self.xIndex=radioBox.GetSelection()

    def onYRadioClick(self, event):
        # gets yIndex of selected button
        radioBox = event.GetEventObject()
        self.yIndex=radioBox.GetSelection()

    def draw_figure(self):
        self.axes.clear()
        for dataObjIndex in self.modelCheckList.GetChecked():
            try:
                self.axes.plot(self.dataObjList[dataObjIndex].dataArray[self.xIndex],
                                 self.dataObjList[dataObjIndex].dataArray[self.yIndex],
                                 self.colour_list[dataObjIndex], 
                                 label=self.dataObjList[dataObjIndex].name)
            except IndexError:
                self.axes.plot(self.dataObjList[dataObjIndex].dataArray[self.xIndex],
                               self.dataObjList[dataObjIndex].dataArray[self.yIndex],
                               label=self.dataObjList[dataObjIndex].name)
            #self.axes.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05),
            #                   ncol=len(self.modelCheckList.GetChecked()),
            #                   )

        if len(self.modelCheckList.GetChecked())>0:
            self.axes.set_xlabel(self.dataObjList[0].axislabels[self.xIndex], 
                                    fontsize=10)
            self.axes.set_ylabel(self.dataObjList[0].axislabels[self.yIndex],
                                    fontsize=10)

            for label in self.axes.get_xticklabels():
                label.set_fontsize(7)
            for label in self.axes.get_yticklabels():
                label.set_fontsize(7)

        if self.logx.IsChecked():
            self.axes.set_xscale('log')
            for label in self.axes.get_xticklabels():
                label.set_fontsize(7)

        if self.logy.IsChecked():
            self.axes.set_yscale('log')
            for label in self.axes.get_yticklabels():
                label.set_fontsize(7)

        self.axes.grid(self.grid.IsChecked())
        self.canvas.draw()

    def on_draw_button(self, event):
        # redraw the canvas
        self.draw_figure()

    def on_logx(self, event):
        # placeholder functions in case I want to do something else
        # decided I don't want plot to redraw until the button is
        # clicke, rather than on selection of log axis
        pass

    def on_logy(self, event):
        pass

    def on_grid(self, event):
        self.draw_figure()


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

    app = tk.PySimpleApp()
    app.frame = VisualizeFrame(dataObj)
    app.frame.Show()
    app.MainLoop()

