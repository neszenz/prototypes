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
    tuple of...
    ...outer boundary vertex loop
    ...list of inner boundary vertex loops
vertex loop:
    list of vertices which are lists of 2 or 3 coordinate values
"""
import datetime
import os
import random #TODO remove

from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB #TODO remove
from OCC.Core.ShapeAnalysis import shapeanalysis
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopTools import TopTools_IndexedMapOfShape
from OCC.Display.SimpleGui import init_display #TODO remove
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer
## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +

PATH_BOX = '../resources/step/boxWithHole.step' #TODO remove
PATH_42 = '../resources/abc/step/00090042/00090042_1cbd5c59fb8c0895c3599d4d_step_007.step' #TODO remove
PATH = PATH_BOX #TODO remove

class SAMPLER_TYPE: # enum for calling sampler factory
    SIMPLE = 0
    NORMALIZED = 1
    string_dict = {
        SIMPLE:'SIMPLE',
        NORMALIZED:'NORMALIZED'
    }

OUTPUT_DIR = os.path.join('tmp', 'sampler')
# to make written file names unique, a timestamp prefix is used
TIMESTAMP = datetime.datetime.now()
PREFIX_SHORT = TIMESTAMP.strftime('%y%b%d_%H%M%S_')
PREFIX_LONG = TIMESTAMP.strftime('%y%b%d_%H%M%S_%f_')

## function body = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +

def generate_mesh_from_edge(edge, shape_maps):
    edge_mesh = list()
    edge_mesh += [42]
    return edge_mesh

def generate_mesh_from_wire(wire, shape_maps):
    wire_mesh = list()
    ex = TopExp_Explorer(wire, TopAbs_EDGE)
    while ex.More():
        wire_mesh += generate_mesh_from_edge(ex.Current(), shape_maps)
        ex.Next()
    return wire_mesh

def generate_mesh_from_face(face, shape_maps):
    _, wire_map, _ = shape_maps
    face_mesh = (list(), list())
    outer_loop, inner_loops = face_mesh
    ex = TopExp_Explorer(face, TopAbs_WIRE)
    outer_wire_id = wire_map.FindIndex(shapeanalysis.OuterWire(face))
    while ex.More():
        wire = ex.Current()
        wire_id = wire_map.FindIndex(wire)
        if wire_id == outer_wire_id:
            outer_loop += generate_mesh_from_wire(wire, shape_maps)
        else:
            inner_loops.append(generate_mesh_from_wire(wire, shape_maps))
        ex.Next()

    return face_mesh

def generate_shape_maps(compound):
    face_map = TopTools_IndexedMapOfShape()
    wire_map = TopTools_IndexedMapOfShape()
    edge_map = TopTools_IndexedMapOfShape()
    topexp.MapShapes(compound, TopAbs_FACE, face_map)
    topexp.MapShapes(compound, TopAbs_WIRE, wire_map)
    topexp.MapShapes(compound, TopAbs_EDGE, edge_map)
    return (face_map, wire_map, edge_map)

def generate_mesh_from_compound(compound):
    model_mesh = list()
    shape_maps = generate_shape_maps(compound)
    ex = TopologyExplorer(compound)
    print('compound has', ex.number_of_faces(), 'faces')
    for face in ex.faces():
        face_mesh = generate_mesh_from_face(face, shape_maps)
        model_mesh.append(face_mesh)
    return model_mesh

def factory(type):
    if type in SAMPLER_TYPE.string_dict:
        print('>> create sampler of type', SAMPLER_TYPE.string_dict[type])
    else:
        raise Exception('factory() error - unknown sampler type')
    def sampler(path):
        print('>> start sampling \'', path, '\'', sep='')
        compound = read_step_file(path)
        return generate_mesh_from_compound(compound)
    return sampler

def write_mesh_to_file(mesh, name):
    print('>> write to file \"', name, '\"', sep='')
    print('mesh:', mesh)

if __name__ == '__main__':
    sampler = factory(SAMPLER_TYPE.SIMPLE)
    mesh = sampler(PATH)
    write_mesh_to_file(mesh, '')
