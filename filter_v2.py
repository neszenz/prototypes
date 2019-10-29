import glob
import os.path
import sys

from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BezierSurface, GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution, GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface, GeomAbs_OtherSurface
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer

INPUT_DIR = os.path.join('..', 'resources', 'abc', 'step')
OUTPUT_DIR = 'tmp/filter_v2'
FILE_EXTENSION = 'step'
FILTER_TYPE = GeomAbs_BezierSurface
IGNORE_FILES = 0
FILE_COUNTER_MAX = 20
file_counter = 0
files_processed = 0
TYPE_STRINGS = {
    GeomAbs_Plane:'GeomAbs_Plane',
    GeomAbs_Cylinder:'GeomAbs_Cylinder',
    GeomAbs_Cone:'GeomAbs_Cone',
    GeomAbs_Sphere:'GeomAbs_Sphere',
    GeomAbs_Torus:'GeomAbs_Torus',
    GeomAbs_BezierSurface:'GeomAbs_BezierSurface',
    GeomAbs_BSplineSurface:'GeomAbs_BSplineSurface',
    GeomAbs_SurfaceOfRevolution:'GeomAbs_SurfaceOfRevolution',
    GeomAbs_SurfaceOfExtrusion:'GeomAbs_SurfaceOfExtrusion',
    GeomAbs_OffsetSurface:'GeomAbs_OffsetSurface',
    GeomAbs_OtherSurface:'GeomAbs_OtherSurface'
}
type_counter = {
    GeomAbs_Plane:0,
    GeomAbs_Cylinder:0,
    GeomAbs_Cone:0,
    GeomAbs_Sphere:0,
    GeomAbs_Torus:0,
    GeomAbs_BezierSurface:0,
    GeomAbs_BSplineSurface:0,
    GeomAbs_SurfaceOfRevolution:0,
    GeomAbs_SurfaceOfExtrusion:0,
    GeomAbs_OffsetSurface:0,
    GeomAbs_OtherSurface:0
}

def print_intro():
    print('+=+=+=+=+=+=+=+ start filtering... +=+=+=+=+=+=+=+')
    print('type:', TYPE_STRINGS[FILTER_TYPE])
    print('from:', INPUT_DIR)

def print_outro():
    print()
    print('+=+=+=+=+=+=+=+ filtering  done! +=+=+=+=+=+=+=+=+')
    print('model faces found in', files_processed, 'files:')
    for item in type_counter:
        print(type_counter[item], TYPE_STRINGS[item], end='')
        if item == FILTER_TYPE:
            print('    <-- filtered type')
        else:
            print('')

def compileFileList():
    file_list = glob.glob(os.path.join(INPUT_DIR, '*.'+FILE_EXTENSION))
    file_list += glob.glob(os.path.join(INPUT_DIR, '*/*.'+FILE_EXTENSION))
    return file_list

def executeOnAllFiles():
    def executeOnFile(file):
        global files_processed
        print('\n>> file', file_counter)
        compound = read_step_file(file, verbosity=False) # always 3 unknown entities
        ex = TopologyExplorer(compound)
        for face in ex.faces():
            type_counter[BRepAdaptor_Surface(face).GetType()] += 1
            if BRepAdaptor_Surface(face).GetType() == FILTER_TYPE:
                os.makedirs(OUTPUT_DIR, exist_ok=True)
                #TODO write model face to new file
                files_processed += 1
    global file_counter
    file_list = compileFileList()
    for file in file_list:
        if file_counter >= FILE_COUNTER_MAX:
            return
        if file_counter == 131: # file no. 91715
            file_counter += 1
            continue
        if file_counter >= IGNORE_FILES:
            executeOnFile(file)
        file_counter += 1

if __name__ == '__main__':
    print_intro()
    executeOnAllFiles()
    print_outro()
