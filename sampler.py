""" 1D mesh sampler for step files
Contains functionality to generate a 1-dimensional mesh from a given path to a
step file.

A step model is made up of model faces, each with an underlying geometry. All
model face boundaries in a step model are circular. Outer boundaries are
oriented counterclockwise and inner ones clockwise.
This leads to the following representation:

model mesh (1D):
    tuple of...
    ...list of vertices which are lists of 3 floats
    ...list of indexed 1-dim. meshes of each model face
    ...surface type string
face mesh (1D):
    tuple of...
    ...outer boundary vertex loop
    ...list of inner boundary vertex loops
vertex loop:
    list of vertex indices
"""
import datetime
import numpy as np
import os
import pickle

from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BezierSurface, GeomAbs_BSplineSurface, GeomAbs_SurfaceOfRevolution, GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface, GeomAbs_OtherSurface
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Core.ShapeAnalysis import shapeanalysis
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_FORWARD, TopAbs_REVERSED
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopTools import TopTools_IndexedMapOfShape
from OCC.Display.SimpleGui import init_display
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer

from mesh1D import SuperVertex, Mesh1D, FILE_EXTENSION

SURFACE_TYPE_STRINGS = {
    GeomAbs_Plane :               'Plane',
    GeomAbs_Cylinder :            'Cylinder',
    GeomAbs_Cone :                'Cone',
    GeomAbs_Sphere :              'Sphere',
    GeomAbs_Torus :               'Torus',
    GeomAbs_BezierSurface :       'BezierSurface',
    GeomAbs_BSplineSurface :      'BSplineSurface',
    GeomAbs_SurfaceOfRevolution : 'SurfaceOfRevolution',
    GeomAbs_SurfaceOfExtrusion :  'SurfaceOfExtrusion',
    GeomAbs_OffsetSurface :       'OffsetSurface',
    GeomAbs_OtherSurface :        'OtherSurface'
}

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
PATH_BOX = '../resources/step/boxWithHole.step'
PATH_SURFFIL = '../resources/step/surface_filling.step'
PATH_FILLET = '../resources/step/fillet.step'
PATH_42 = '../resources/abc/step/00090042/00090042_1cbd5c59fb8c0895c3599d4d_step_007.step'
PATH_111 = '../resources/abc/step/00090111/00090111_7eac35f07183d39b4da202d3_step_000.step'
PATH_99999 = '../resources/abc/step/00099999/00099999_df6629a908dec75f8a69bda7_step_001.step'
PATH_98613 = '../resources/abc/step/00098613/00098613_56d3ec39e4b0747e94b812ee_step_007.step'
INPUT_PATH = PATH_98613

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
PREFIX_SHORT = TIMESTAMP.strftime('%y%m%d_%H%M%S')
PREFIX_LONG = TIMESTAMP.strftime('%y%m%d_%H%M%S_%f')
#config
NUMBER_OF_SAMPLES = 20
INCLUDE_OUTER_WIRES = True
INCLUDE_INNER_WIRES = True
SIMPLIFY_VERTEX_LIST = False # removes doubly vertices; performance impact

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
            if i == NUMBER_OF_SAMPLES-1:
                parameter = lp
            else:
                parameter = fp + i * (p_length / (NUMBER_OF_SAMPLES-1))
            point = curve.Value(parameter)
            sv = SuperVertex(x=point.X(), y=point.Y(), z=point.Z())
            edge_mesh.append(sv)
        return edge_mesh
    _, _, edge_map = shape_maps
    edge_meshes = {}
    number_of_edges = TopologyExplorer(compound).number_of_edges()
    for edge_id in range(1, number_of_edges+1):
        edge = edge_map.FindKey(edge_id)
        edge_mesh = generate_mesh_from_edge(edge, shape_maps)
        edge_meshes[edge_id] = edge_mesh
    return edge_meshes

def indexed_edge_meshes_from(edge_meshes):
    def get_vertex_index(vertices, vertex):
        for iVertex in range(0, len(vertices)):
            v = vertices[iVertex]
            if np.allclose(v, vertex):
                return iVertex
        return -1
    def removeDoubles(vertices):
        remove_map = []
        number_of_vertices = len(vertices)
        i = 0
        while i < number_of_vertices:
            vi = vertices[i]
            j = i+1
            while j < number_of_vertices:
                vj = vertices[j]
                if np.allclose(vi, vj):
                    remove_map.append((j, i))
                    del vertices[j]
                    number_of_vertices -= 1
                else:
                    j += 1
            i += 1
        return remove_map
    def updateIndeces(edge_meshes, remove_map):
        for edge_id in edge_meshes:
            edge_mesh = edge_meshes[edge_id]
            for iIndex in range(0, len(edge_mesh)):
                for pair in remove_map:
                    i_removed, i_mapped = pair
                    if edge_mesh[iIndex] == i_removed:
                        edge_mesh[iIndex] = i_mapped
                    elif edge_mesh[iIndex] > i_removed:
                        edge_mesh[iIndex] -= 1
    vertices = []
    indexed_edge_meshes = {}
    for edge_id in edge_meshes:
        edge_mesh = edge_meshes[edge_id]
        indexed_edge_mesh = []
        for vertex in edge_mesh:
            vertex_index = len(vertices)
            vertices.append(vertex)
            assert vertex_index >= 0
            indexed_edge_mesh.append(vertex_index)
        assert len(edge_mesh) == len(indexed_edge_mesh)
        indexed_edge_meshes[edge_id] = indexed_edge_mesh
    assert len(edge_meshes) == len(indexed_edge_meshes)
    if SIMPLIFY_VERTEX_LIST:
        old_num_of_vertices = len(vertices)
        remove_map = removeDoubles(vertices)
        assert old_num_of_vertices == len(vertices) + len(remove_map)
        print('simplifyed vertex list by', round(100*len(remove_map)/old_num_of_vertices), '%')
        updateIndeces(indexed_edge_meshes, remove_map)
    else:
        print('vertex list not simplifyed, can contain doubly vertices')
    return vertices, indexed_edge_meshes

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
        face_type = BRepAdaptor_Surface(face).GetType()
        face_framework = (list(), list(), SURFACE_TYPE_STRINGS[face_type])
        outer_loop, inner_loops, _ = face_framework
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
    print('model compound contains', ex.number_of_faces(), 'faces', end='')
    print(',', ex.number_of_wires(), 'wires', end='')
    print(',', ex.number_of_edges(), 'edges')
    for face in ex.faces():
        face_framework = generate_face_framework(face, shape_maps)
        model_framework.append(face_framework)
    return model_framework

def combine_framework_with_edge_meshes(framework, edge_meshes):
    def process_wire_framework(framework, edge_meshes):
        wire_mesh = list()
        for edge_info in framework:
            edge_id, need_reversed = edge_info
            edge_mesh = edge_meshes[edge_id].copy()
            if need_reversed:
                edge_mesh.reverse()
            edge_mesh.pop()
            wire_mesh += edge_mesh
        return wire_mesh
    face_meshes = list()
    for face in framework:
        outer_wire, inner_wirees, face_type = face
        new_outer = process_wire_framework(outer_wire, edge_meshes)
        new_inner = list()
        for inner_wire in inner_wirees:
            new_inner.append(process_wire_framework(inner_wire, edge_meshes))
        face_meshes.append((new_outer, new_inner, face_type))
    return face_meshes

def factory(type):
    if type in SAMPLER_TYPE.string_dict:
        print('>> create sampler of type', SAMPLER_TYPE.string_dict[type])
    else:
        raise Exception('factory() error - unknown sampler type')
    #TODO different sampler methods
    def sampler(path):
        print('>> loading step file \'', path, '\'', sep='')
        compound = read_step_file(path, verbosity=False)
        shape_maps = generate_shape_maps(compound)
        print('>> start sampling...')
        framework = generate_mesh_framework(compound, shape_maps)
        edge_meshes = sample_all_edges(compound, shape_maps)
        vertices, indexed_edge_meshes = indexed_edge_meshes_from(edge_meshes)
        face_meshes = combine_framework_with_edge_meshes(framework, indexed_edge_meshes)
        return Mesh1D(path, vertices, face_meshes)
    return sampler

def noDoublyLoopInsertions(mesh1D):
    def checkWire(wire_mesh):
        ok = True
        for i in range(0, len(wire_mesh)):
            iIndex = wire_mesh[i]
            vi = mesh1D.vertices[iIndex]
            for j in range(i+1, len(wire_mesh)):
                jIndex = wire_mesh[j]
                vj = mesh1D.vertices[jIndex]
                if np.allclose(vi, vj):
                    print('vertex', vi, 'occures both at', i, 'and', j)
                    ok = False
        return ok
    ok = True
    for i in range(0, mesh1D.number_of_faces()):
        face_mesh = mesh1D.face_meshes[i]
        print('\nface', i)
        outer_wire_mesh, inner_wire_meshes = face_mesh
        if not checkWire(outer_wire_mesh):
            ok = False
        for inner_wire_mesh in inner_wire_meshes:
            if not checkWire(inner_wire_mesh):
                ok = False
    return ok

def write_mesh_to_file(mesh1D):
    def generate_output_file_path():
        name = PREFIX_SHORT + FILE_EXTENSION
        path = os.path.join(OUTPUT_DIR, name)
        if os.path.exists(path):
            name = PREFIX_LONG + FILE_EXTENSION
        assert not os.path.exists(path)
        return path
    path = generate_output_file_path()
    print('>> write to file \"', path, '\"', sep='')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pickle.dump(mesh1D, open(path, 'wb'))

if __name__ == '__main__':
    sampler = factory(SAMPLER_TYPE.SIMPLE)
    mesh1D = sampler(INPUT_PATH)
    # print('no doubly vertices per loop:', noDoublyLoopInsertions(mesh1D))
    write_mesh_to_file(mesh1D)
