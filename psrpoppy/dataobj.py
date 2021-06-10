#!/usr/bin/python

import os

import pickle


class DataObj:
    def __init__(self, name, labelDict, dataDict, size):
        self.name = name

        self.labelDict = labelDict
        self.dataDict = dataDict
        self.count = size


def makeDataObj(modelname=None,
                frac=None,
                extnList=['.model', '.results']):
    """
    Function to get different parameters of a population file

    """

    dataObjList = []

    if modelname is None:
        for filename in os.listdir(os.getcwd()):
            for extn in extnList:
                if filename.endswith(extn):
                    dataObjList.append(readfile_return_dataobj(filename))
    else:
        dataObjList.append(readfile_return_dataobj(modelname))

    if len(dataObjList) == 0:
        raise ViewException('No files found matching the extensions')

    # sort list from largest -> smallest so plotting is always done that way
    dataObjList.sort(key=lambda x: x.count, reverse=True)

    return dataObjList


def readfile_return_dataobj(filename):
    """Read in the given file and make a dataobj"""
    try:
        f = open(filename, 'rb')
    except IOError:
        print(("Could not open file {0}.".format(filename)))
        return

    try:
        pop = pickle.load(f)
    except pickle.UnpicklingError:
        print(("File {0} could not be unpickled!".format(filename)))
        return
    f.close()

    labelDict, dataDict, size = pop.make_plotting_dicts()

    return DataObj(filename,
                   labelDict,
                   dataDict,
                   pop.size())
