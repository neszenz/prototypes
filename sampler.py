""" 1D mesh sampler for STEP files
Contains functionality to generate a 1-dimensional mesh from a given path to a
STEP file.

A STEP model is made up of model faces, each with an underlying geometry. All
model face boundaries in a step model are circular. Outer boundaries are
oriented counterclockwise and inner ones clockwise.
This leads to the representation as described in meshkD.py.
"""
import datetime
import numpy as np
import os
import pickle

from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Core.ShapeAnalysis import shapeanalysis
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_FORWARD, TopAbs_REVERSED
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopTools import TopTools_IndexedMapOfShape
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer

import paths
from meshkD import SuperVertex, MeshkD

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_PATH = paths.PATH_BOX

class SAMPLER_TYPE: # enum for calling sampler factory
    SIMPLE = 0
    NORMALIZED = 1
    string_dict = {
        SIMPLE:'SIMPLE',
        NORMALIZED:'NORMALIZED'
    }

OUTPUT_DIR = 'tmp'
# to make written file names unique, a timestamp prefix is used
TIMESTAMP = datetime.datetime.now()
OUTPUT_PREFIX = TIMESTAMP.strftime('%y%m%d_%H%M%S')

NUMBER_OF_SAMPLES = 10 # only for SIMPLE sampling method
INCLUDE_OUTER_WIRES = True
INCLUDE_INNER_WIRES = True
SIMPLIFY_VERTEX_LIST = False # removes doubly vertices (cylinders); big performance impact

## functions = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
def generate_shape_maps(compound):
    face_map = TopTools_IndexedMapOfShape()
    wire_map = TopTools_IndexedMapOfShape()
    edge_map = TopTools_IndexedMapOfShape()
    topexp.MapShapes(compound, TopAbs_FACE, face_map)
    topexp.MapShapes(compound, TopAbs_WIRE, wire_map)
    topexp.MapShapes(compound, TopAbs_EDGE, edge_map)
    return (face_map, wire_map, edge_map)

# same structure as model mesh, but wires consist of a list of edge_info tuples
# for later replacement during actual sampling:
# edge_info = (wire, edge)
def generate_mesh_framework(compound, shape_maps):
    def generate_face_framework(face, shape_maps):
        def generate_wire_framework(wire):
            wire_framework = []
            ex = TopExp_Explorer(wire, TopAbs_EDGE)
            while ex.More():
                edge = ex.Current()
                edge_info = (wire, edge)
                wire_framework.append(edge_info)
                ex.Next()
            return wire_framework
        _, wire_map, _ = shape_maps

        wire_frameworks = []
        outer_loop, inner_loops = [], []

        ex = TopExp_Explorer(face, TopAbs_WIRE)
        outer_wire_id = wire_map.FindIndex(shapeanalysis.OuterWire(face))
        while ex.More():
            wire = ex.Current()
            wire_id = wire_map.FindIndex(wire)
            if wire_id == outer_wire_id and INCLUDE_OUTER_WIRES:
                outer_loop += generate_wire_framework(wire)
            elif wire_id != outer_wire_id and INCLUDE_INNER_WIRES:
                inner_loops.append(generate_wire_framework(wire))
            else:
                pass
            ex.Next()

        wire_frameworks.append(outer_loop)
        wire_frameworks += inner_loops

        return wire_frameworks, face
    model_framework = []
    ex = TopologyExplorer(compound)
    print('model compound contains', ex.number_of_faces(), 'faces', end='')
    print(',', ex.number_of_wires(), 'wires', end='')
    print(',', ex.number_of_edges(), 'edges')
    cnt = 1 #TODO used for testing
    for face in ex.faces():
        if cnt != 1: #TODO used for testing
            cnt += 1 #TODO used for testing
            # continue #TODO used for testing
        face_framework = generate_face_framework(face, shape_maps)
        model_framework.append(face_framework)
        cnt += 1 #TODO used for testing
    return model_framework

def edge_sampler_simple(edge_info, face, shape_maps):
    wire, edge = edge_info
    face_map, _, _ = shape_maps
    edge_mesh = []
    if NUMBER_OF_SAMPLES < 2:
        return edge_mesh
    surface = BRepAdaptor_Surface(face)
    curve, fp, lp = BRep_Tool.CurveOnSurface(edge, face)
    curve = BRepAdaptor_Curve2d(edge, face) # should be the same, but for consistency w/ rest of code
    p_length = lp - fp
    for i in range(0, NUMBER_OF_SAMPLES):
        if i == NUMBER_OF_SAMPLES-1:
            parameter = lp
        else:
            parameter = fp + i*(p_length / (NUMBER_OF_SAMPLES-1))
        p2d = curve.Value(parameter)
        p3d = surface.Value(p2d.X(), p2d.Y())
        sv = SuperVertex(x=p3d.X(), y=p3d.Y(), z=p3d.Z(), u=p2d.X(), v=p2d.Y())
        sv.face_id = face_map.FindIndex(face)
        sv.face = face
        sv.edges_with_p = [(edge, parameter)]
        edge_mesh.append(sv)
    # here, the edges orientation are made consistent with that of the wire
    if edge.Orientation() != wire.Orientation(): # corrects ori to keep edges consistent in wire
        edge_mesh.reverse()
    return edge_mesh
def sample_edges_in_framework(framework, shape_maps, sampler_type):
    def add_wire(face_mesh, wire_framework, shape_maps, sampler_type):
        vertices, wire_meshes, _, face = face_mesh

        wire_vertices = []
        wire_mesh = []

        edge_mesh_loop = []

        # collect all edge meshes
        for edge_info in wire_framework:
            if sampler_type == SAMPLER_TYPE.SIMPLE:
                edge_mesh = edge_sampler_simple(edge_info, face, shape_maps)
            else:
                pass #TODO different sampler methods

            edge_mesh_loop.append(edge_mesh)

        # merge edge meshes together into a wire vertex loop
        for i in range(len(edge_mesh_loop)):
            curr_em = edge_mesh_loop[i]
            next_em = edge_mesh_loop[(i+1) % len(edge_mesh_loop)]
            if curr_em[-1] == next_em[0]:
                assert not curr_em[-1].edges_with_p is None
                assert not next_em[0].edges_with_p is None
                assert len(curr_em[-1].edges_with_p) == 1
                assert len(next_em[0].edges_with_p) == 1

                # merge edge list in connection and insert current into wire_mesh
                next_em[0].edges_with_p = curr_em[-1].edges_with_p + next_em[0].edges_with_p
                curr_em.pop()

                wire_vertices += curr_em
            else:
                print('sample_edges_in_framework() error - wire connectivity degenerated')
                wire_vertices += curr_em
                return

        # construct indexed wire mesh
        offset = len(vertices)
        wire_mesh = [offset+i for i in range(len(wire_vertices))]

        # here, the model face orientation gets adapted by the mesh
        if face.Orientation() == TopAbs_REVERSED:
            wire_mesh.reverse()

        # insert result into face_mesh
        vertices += wire_vertices
        wire_meshes.append(wire_mesh)

        return
    def reverse_u_of_all_vertices(face_mesh):
        vertices, _, _, _ = face_mesh

        for sv in vertices:
            sv.reverse_u()

        return
    face_meshes = []

    for face_framework in framework:
        wire_frameworks, face = face_framework
        vertices = []
        wire_meshes = []
        face_mesh = (vertices, wire_meshes, [], face)

        for wire_framework in wire_frameworks:
            add_wire(face_mesh, wire_framework, shape_maps, sampler_type)

        if face.Orientation() == TopAbs_REVERSED:
            reverse_u_of_all_vertices(face_mesh)

        face_meshes.append(face_mesh)

    return face_meshes

def factory(sampler_type):
    if sampler_type in SAMPLER_TYPE.string_dict:
        print('>> creating sampler of type', SAMPLER_TYPE.string_dict[sampler_type])
    else:
        raise Exception('factory() error - unknown sampler type')
    def sampler(path):
        print('>> loading step file \'', path, '\'', sep='')
        compound = read_step_file(path, verbosity=False)
        shape_maps = generate_shape_maps(compound)
        print('>> generating framework and sampling...')
        framework = generate_mesh_framework(compound, shape_maps)
        face_meshes = sample_edges_in_framework(framework, shape_maps, sampler_type)
        return MeshkD(path, face_meshes)
    return sampler

def noDoublyLoopInsertions(mesh1D):
    def checkWire(vertices, wire_mesh):
        ok = True
        for i in range(0, len(wire_mesh)):
            i_index = wire_mesh[i]
            svi = vertices[i_index]
            for j in range(i+1, len(wire_mesh)):
                j_index = wire_mesh[j]
                svj = vertices[j_index]
                if svi == svj:
                    print('vertex', svi, 'occures both at', i, 'and', j, end='')
                    print(' (i_index:', i_index, ', j_index:', j_index, ')', sep='')
                    ok = False
        return ok
    ok = True

    for i in range(mesh1D.number_of_faces()):
        print('face', i)
        vertices, wire_meshes, _, _ = mesh1D.face_meshes[i]

        for wire_mesh in wire_meshes:
            if not checkWire(vertices, wire_mesh):
                ok = False

    return ok

def write_mesh_to_file(mesh1D):
    def generate_output_file_path():
        name = OUTPUT_PREFIX + MeshkD.FILE_EXTENSION
        path = os.path.join(OUTPUT_DIR, name)
        if os.path.exists(path):
            raise Exception('output file name already exists:', path)
        assert not os.path.exists(path)
        return path
    path = generate_output_file_path()
    print('>> write to file \"', path, '\"', sep='')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pickle.dump(mesh1D, open(path, 'wb'))

if __name__ == '__main__':
    sampler = factory(SAMPLER_TYPE.SIMPLE)
    mesh1D = sampler(INPUT_PATH)
    print('no doubly vertices per loop:', noDoublyLoopInsertions(mesh1D))
    write_mesh_to_file(mesh1D)
