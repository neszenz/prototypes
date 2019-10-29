import glob
import os.path
import random

from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane, GeomAbs_BezierSurface, GeomAbs_BezierCurve
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.Interface import Interface_Static_SetCVal
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Display.SimpleGui import init_display
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer

OUTPUT_DIR = 'tmp/filter_v1'
OUTPUT_EXTENSION = '.step'
FILTER_TYPE = GeomAbs_Plane

def filter_all_step_files(event=None):
    inputFiles = glob.glob(os.path.join('../resources/step', '*.step'))
    for inputFile in inputFiles:
        compound = read_step_file(inputFile)
        index = 0
        for face in TopologyExplorer(compound).faces():
            if BRepAdaptor_Surface(face).GetType() == FILTER_TYPE:
                fileName = os.path.join(OUTPUT_DIR, str(index) + OUTPUT_EXTENSION)
                os.makedirs(OUTPUT_DIR, exist_ok=True)
                write_step_file(face, fileName)
                index += 1

def import_filtered_output(event=None):
    outputFiles = glob.glob(os.path.join(OUTPUT_DIR, '*' + OUTPUT_EXTENSION))
    if len(outputFiles) > 0:
        display.EraseAll()
        for outputFile in outputFiles:
            compound = read_step_file(outputFile)
            color = Quantity_Color(random.random(),
                                   random.random(),
                                   random.random(),
                                   Quantity_TOC_RGB)
            display.DisplayColoredShape(compound, color)
        display.FitAll()
    else:
        display.EraseAll()
        print("import_output_file() error - No output file(s) found")

if __name__ == '__main__':
    display, start_display, add_menu, add_function_to_menu = init_display()
    add_menu('STEP import')
    add_function_to_menu('STEP import', filter_all_step_files)
    add_function_to_menu('STEP import', import_filtered_output)
    start_display()
