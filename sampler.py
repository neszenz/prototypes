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

from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB #TODO remove
from OCC.Core.ShapeAnalysis import shapeanalysis
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_FORWARD, TopAbs_REVERSED
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopTools import TopTools_IndexedMapOfShape
from OCC.Display.SimpleGui import init_display #TODO remove
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer
## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +

PATH_BOX = '../resources/step/boxWithHole.step'
PATH_X = '../resources/step/surface_filling.step'
PATH_42 = '../resources/abc/step/00090042/00090042_1cbd5c59fb8c0895c3599d4d_step_007.step'
PATH_111 = '../resources/abc/step/00090111/00090111_7eac35f07183d39b4da202d3_step_000.step'
PATH_99999 = '../resources/abc/step/00099999/00099999_df6629a908dec75f8a69bda7_step_001.step'
INPUT_PATH = PATH_99999

class SAMPLER_TYPE: # enum for calling sampler factory
    SIMPLE = 0
    NORMALIZED = 1
    string_dict = {
        SIMPLE:'SIMPLE',
        NORMALIZED:'NORMALIZED'
    }

OUTPUT_DIR = os.path.join('tmp', 'sampler')
LINES_EXTENSION = '.lines'
# to make written file names unique, a timestamp prefix is used
TIMESTAMP = datetime.datetime.now()
PREFIX_SHORT = TIMESTAMP.strftime('%y%m%d_%H%M%S')
PREFIX_LONG = TIMESTAMP.strftime('%y%m%d_%H%M%S_%f')
#config
NUMBER_OF_SAMPLES = 20
INCLUDE_OUTER_WIRES = True
INCLUDE_INNER_WIRES = True

## function body = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +

def generate_shape_maps(compound):
    face_map = TopTools_IndexedMapOfShape()
    wire_map = TopTools_IndexedMapOfShape()
    edge_map = TopTools_IndexedMapOfShape()
    topexp.MapShapes(compound, TopAbs_FACE, face_map)
    topexp.MapShapes(compound, TopAbs_WIRE, wire_map)
    topexp.MapShapes(compound, TopAbs_EDGE, edge_map)
    return (face_map, wire_map, edge_map)

def sample_all_edges(compound, shape_maps):
    def generate_mesh_from_edge(edge, shape_maps):
        edge_mesh = list()
        curve = BRepAdaptor_Curve(edge)
        fp = curve.FirstParameter()
        lp = curve.LastParameter()
        p_length = lp - fp
        assert NUMBER_OF_SAMPLES > 1
        for i in range(0, NUMBER_OF_SAMPLES):
            parameter = fp + i * (p_length / (NUMBER_OF_SAMPLES-1))
            if i == NUMBER_OF_SAMPLES-1:
                parameter = lp
            point = curve.Value(parameter)
            edge_mesh.append([point.X(), point.Y(), point.Z()])
        return edge_mesh
    _, _, edge_map = shape_maps
    edge_meshes = {}
    number_of_edges = TopologyExplorer(compound).number_of_edges()
    for edge_id in range(1, number_of_edges+1):
        edge = edge_map.FindKey(edge_id)
        edge_mesh = generate_mesh_from_edge(edge, shape_maps)
        edge_meshes[edge_id] = edge_mesh
    return edge_meshes

# same structure as model mesh, but wires consist of a list of tupels containing
# infos about the edge for later matching with separately sampled edge meshes
def generate_mesh_framework(compound, shape_maps):
    def generate_face_framework(face, shape_maps):
        def generate_wire_framework(wire, face_ori, shape_maps):
            _, _, edge_map = shape_maps
            wire_framework = list()
            ex = TopExp_Explorer(wire, TopAbs_EDGE)
            while ex.More():
                edge = ex.Current()
                edge_id = edge_map.FindIndex(edge)
                need_reversed = edge.Orientation() != face_ori
                edge_infos = (edge_id, need_reversed)
                wire_framework.append(edge_infos)
                ex.Next()
            return wire_framework
        _, wire_map, _ = shape_maps
        face_framework = (list(), list())
        outer_loop, inner_loops = face_framework
        ex = TopExp_Explorer(face, TopAbs_WIRE)
        outer_wire_id = wire_map.FindIndex(shapeanalysis.OuterWire(face))
        while ex.More():
            wire = ex.Current()
            wire_id = wire_map.FindIndex(wire)
            wire_ori = wire.Orientation()
            if wire_id == outer_wire_id and INCLUDE_OUTER_WIRES:
                outer_loop += generate_wire_framework(wire, wire_ori, shape_maps)
            elif wire_id != outer_wire_id and INCLUDE_INNER_WIRES:
                inner_loops.append(generate_wire_framework(wire, wire_ori, shape_maps))
            else:
                pass
            ex.Next()

        return face_framework
    model_framework = list()
    ex = TopologyExplorer(compound)
    print('compound has', ex.number_of_faces(), 'faces')
    print('            ', ex.number_of_wires(), 'wires')
    print('            ', ex.number_of_edges(), 'edges')
    for face in ex.faces():
        face_framework = generate_face_framework(face, shape_maps)
        model_framework.append(face_framework)
    return model_framework

def combine_framework_with_edge_meshes(framework, edge_meshes):
    def process_wire_framework(framework, edge_meshes):
        wire_mesh = list()
        for edge_info in framework:
            edge_mesh = edge_meshes[edge_info[0]].copy()
            if edge_info[1]:
                edge_mesh.reverse()
            edge_mesh.pop()
            wire_mesh += edge_mesh
        framework.clear()
        framework += wire_mesh
    for face in framework:
        outer_wire = face[0]
        process_wire_framework(outer_wire, edge_meshes)
        for inner_wire in face[1]:
            process_wire_framework(inner_wire, edge_meshes)
    return framework

def factory(type):
    if type in SAMPLER_TYPE.string_dict:
        print('>> create sampler of type', SAMPLER_TYPE.string_dict[type])
    else:
        raise Exception('factory() error - unknown sampler type')
    def sampler(path):
        print('>> loading step file \'', path, '\'', sep='')
        compound = read_step_file(path, verbosity=False)
        print('>> start sampling...')
        shape_maps = generate_shape_maps(compound)
        #TODO different sampler methods
        edge_meshes = sample_all_edges(compound, shape_maps)
        framework = generate_mesh_framework(compound, shape_maps)
        return combine_framework_with_edge_meshes(framework, edge_meshes)
    return sampler

def write_mesh_to_file(mesh):
    def generate_output_file_path():
        name = PREFIX_SHORT + LINES_EXTENSION
        path = os.path.join(OUTPUT_DIR, name)
        if os.path.exists(path):
            name = PREFIX_LONG + LINES_EXTENSION
        assert not os.path.exists(path)
        return path
    def write_wire_mesh(wire_mesh, output):
        for iVertex in range(0, len(wire_mesh)):
            v0 = wire_mesh[iVertex]
            v1 = wire_mesh[(iVertex + 1) % len(wire_mesh)]
            output.write('%s %s %s' % (v0[0], v0[1], v0[2]))
            output.write(' ')
            output.write('%s %s %s' % (v1[0], v1[1], v1[2]))
            output.write('\n')
    path = generate_output_file_path()
    print('>> write to file \"', path, '\"', sep='')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    with open(path, 'a') as output:
        for face_mesh in mesh:
            outer_wire_mesh = face_mesh[0]
            write_wire_mesh(outer_wire_mesh, output)
            for inner_wire_mesh in face_mesh[1]:
                write_wire_mesh(inner_wire_mesh, output)

if __name__ == '__main__':
    sampler = factory(SAMPLER_TYPE.SIMPLE)
    mesh = sampler(INPUT_PATH)
    write_mesh_to_file(mesh)
