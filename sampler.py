""" 1D mesh sampler for step files
Contains functionality to generate a 1-dimensional mesh from a given path to a
step file.

A step model is made up of model faces, each with an underlying geometry. All
model face boundaries in a step model are circular. Outer boundaries are
oriented counterclockwise and inner ones clockwise.
This leads to the following representation:

model mesh (1D):
    list of 1D meshes of each model face
face mesh (1D):
    outer boundary vertex loop
    list of inner boundary vertex loops
vertex loop:
    list of vertices which are lists of 2 or 3 coordinate values
"""

import datetime
import os
import random #TODO remove

from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB #TODO remove
from OCC.Display.SimpleGui import init_display #TODO remove
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer

PATH_BOX = '../resources/step/boxWithHole.step' #TODO remove
PATH_42 = '../resources/abc/step/00090042/00090042_1cbd5c59fb8c0895c3599d4d_step_007.step' #TODO remove
PATH = PATH_BOX #TODO remove

OUTPUT_DIR = os.path.join('tmp', 'sampler')
# to make written file names unique, a timestamp prefix is used
TIMESTAMP = datetime.datetime.now()
PREFIX_SHORT = TIMESTAMP.strftime('%y%b%d_%H%M%S_')
PREFIX_LONG = TIMESTAMP.strftime('%y%b%d_%H%M%S_%f_')

def process_solid(solid):
    print(PREFIX_SHORT)
    print(PREFIX_LONG)
    pass

def process_file(path):
    solids = read_step_file(path, return_as_shapes=True)
    for solid in solids:
        process_solid(solid)

def open_file(display): #TODO remove
    solids = read_step_file(PATH, return_as_shapes=True)
    display.EraseAll()
    for solid in solids:
        color = Quantity_Color(random.random(),
                               random.random(),
                               random.random(),
                               Quantity_TOC_RGB)
        display.DisplayColoredShape(solid, color)
    display.FitAll()
    pass

if __name__ == '__main__':
    display, start_display, _, _ = init_display() #TODO remove
    open_file(display) #TODO remove
    process_file(PATH)
    # start_display() #TODO remove
