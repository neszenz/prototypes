""" 1D mesh sampler for STEP files
Contains functionality to generate a 1-dimensional mesh from a given path to a
STEP file.

A STEP model is made up of model faces, each with an underlying geometry. All
model face boundaries in a step model are circular. Outer boundaries are
oriented counterclockwise and inner ones clockwise.
This leads to the representation as described in meshkD.py.
"""
import numpy as np

from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCC.Core.gp import gp_Pnt, gp_Vec
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Core.ShapeAnalysis import shapeanalysis
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_FORWARD, TopAbs_REVERSED
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopTools import TopTools_IndexedMapOfShape
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer

import paths
from meshkD import SuperVertex, MeshkD, write_to_file

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
class SAMPLER_TYPES: # enum for calling sampler factory
    SIMPLE = 0
    string_dict = {
        SIMPLE:'SIMPLE',
    }

SAMPLER_TYPE = SAMPLER_TYPES.SIMPLE
NUMBER_OF_SAMPLES = 10 # only for SIMPLE sampling method
MIN_NUMBER_OF_SAMPLES = 3 # prevent multiple shared edges when edes get simplified
INCLUDE_OUTER_WIRES = True
INCLUDE_INNER_WIRES = True
REMOVE_SINGULARITIES = True
SIMPLIFY_LINEAR_EDGES = False
REORDER_WIRES_FOR_CLOSEST_ENDPOINTS = True

# __main__ config
INPUT_PATH = paths.STEP_42
OUTPUT_DIR = paths.DIR_TMP

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
    def collect_interface(edge, face):
        surface = BRepAdaptor_Surface(face)
        curve2d = BRepAdaptor_Curve2d(edge, face)
        curve3d = BRepAdaptor_Curve(edge)

        fp = curve2d.FirstParameter()
        lp = curve2d.LastParameter()
        assert fp == curve3d.FirstParameter()
        assert lp == curve3d.LastParameter()
        p_length = lp - fp

        return surface, curve2d, curve3d, fp, lp, p_length
    def sample_derivative(curve3d, parameter):
        _ = gp_Pnt()
        d1 = gp_Vec()
        curve3d.D1(parameter, _, d1)

        return np.array((d1.X(), d1.Y(), d1.Z()))
    def sample_supervertex(surface, curve2d, parameter):
        p2d = curve2d.Value(parameter)
        p3d = surface.Value(p2d.X(), p2d.Y())

        sv = SuperVertex(x=p3d.X(), y=p3d.Y(), z=p3d.Z(), u=p2d.X(), v=p2d.Y())

        return sv
    def remove_singularities(edge_mesh, derivatives):
        sv_index = 1
        while sv_index < len(edge_mesh):
            sv_last = edge_mesh[sv_index-1]
            sv_curr = edge_mesh[sv_index]

            if np.allclose(sv_last.XYZ_vec3(), sv_curr.XYZ_vec3()):
                del edge_mesh[sv_index]
                del derivatives[sv_index]
            else:
                sv_index += 1

        return
    def simplify_linear_segments(edge_mesh, derivatives):
        assert len(edge_mesh) == len(derivatives)

        i = 1
        while i+1 < len(edge_mesh) and len(edge_mesh) >= MIN_NUMBER_OF_SAMPLES+1:
            d_last = derivatives[i-1]
            d_curr = derivatives[i]
            d_next = derivatives[i+1]

            if np.allclose(d_last, d_curr) and np.allclose(d_curr, d_next):
                del edge_mesh[i]
                del derivatives[i]
            else:
                i += 1

        return
    wire, edge = edge_info
    face_map, _, _ = shape_maps

    edge_mesh = []
    derivatives = []

    if NUMBER_OF_SAMPLES < MIN_NUMBER_OF_SAMPLES:
        return edge_mesh

    surface, curve2d, curve3d, fp, lp, p_length = collect_interface(edge, face)

    for i in range(0, NUMBER_OF_SAMPLES):
        if i == NUMBER_OF_SAMPLES-1:
            parameter = lp
        else:
            parameter = fp + i*(p_length / (NUMBER_OF_SAMPLES-1))

        d1 = sample_derivative(curve3d, parameter)
        derivatives.append(d1)

        sv = sample_supervertex(surface, curve2d, parameter)
        sv.face_id = face_map.FindIndex(face)
        sv.face = face
        sv.edges_with_p = [(edge, parameter)]
        edge_mesh.append(sv)

    if REMOVE_SINGULARITIES:
        remove_singularities(edge_mesh, derivatives)

    if SIMPLIFY_LINEAR_EDGES:
        simplify_linear_segments(edge_mesh, derivatives)

    # here, the edges orientation are made consistent with that of the wire
    if edge.Orientation() != wire.Orientation(): # corrects ori to keep edges consistent in wire
        edge_mesh.reverse()

    return edge_mesh
def sample_edges_in_framework(framework, shape_maps):
    def add_wire(face_mesh, wire_framework, shape_maps):
        def order_edge_meshes(edge_mesh_loop):
            new_edge_mesh_loop = []

            new_edge_mesh_loop.append(edge_mesh_loop[0])
            del edge_mesh_loop[0]

            # always append edge mesh with closest start point
            while len(edge_mesh_loop) > 0:
                curr_em = new_edge_mesh_loop[-1]
                clos_i = -1
                clos_dist = float('inf')

                for i in range(len(edge_mesh_loop)):
                    tmp_em = edge_mesh_loop[i]
                    tmp_dist = np.linalg.norm(curr_em[-1].XYZ_vec3() - tmp_em[0].XYZ_vec3())
                    if tmp_dist < clos_dist:
                        clos_i = i
                        clos_dist = tmp_dist

                new_edge_mesh_loop.append(edge_mesh_loop[clos_i])
                del edge_mesh_loop[clos_i]

            return new_edge_mesh_loop
        vertices, wire_meshes, _, face = face_mesh

        wire_vertices = []
        wire_mesh = []

        edge_mesh_loop = []

        # collect all edge meshes
        for edge_info in wire_framework:
            if SAMPLER_TYPE == SAMPLER_TYPES.SIMPLE:
                edge_mesh = edge_sampler_simple(edge_info, face, shape_maps)
            else:
                pass #TODO different sampler methods

            if len(edge_mesh) >= MIN_NUMBER_OF_SAMPLES:
                edge_mesh_loop.append(edge_mesh)

        # repair certain cases of degenerated order and geometry of wires
        if REORDER_WIRES_FOR_CLOSEST_ENDPOINTS:
            edge_mesh_loop = order_edge_meshes(edge_mesh_loop)

        # merge edge meshes together into a wire vertex loop
        for i in range(len(edge_mesh_loop)):
            curr_em = edge_mesh_loop[i]
            next_em = edge_mesh_loop[(i+1) % len(edge_mesh_loop)]

            assert not curr_em[-1].edges_with_p is None
            assert not next_em[0].edges_with_p is None
            assert len(curr_em[-1].edges_with_p) == 1
            assert len(next_em[0].edges_with_p) == 1

            if not np.allclose(curr_em[-1].XYZ_vec3(), next_em[0].XYZ_vec3()):
                print('sample_edges_in_framework() warning - wire connectivity may be degenerated')

            # merge edge list in connection and insert current into wire_mesh
            next_em[0].edges_with_p = curr_em[-1].edges_with_p + next_em[0].edges_with_p
            curr_em.pop()

            wire_vertices += curr_em

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
            add_wire(face_mesh, wire_framework, shape_maps)

        if face.Orientation() == TopAbs_REVERSED:
            reverse_u_of_all_vertices(face_mesh)

        face_meshes.append(face_mesh)

    return face_meshes

def sample(path):
    print('>> start sampler of type', SAMPLER_TYPES.string_dict[SAMPLER_TYPE])
    print('>> loading step file \'', path, '\'', sep='')
    compound = read_step_file(path, verbosity=False)
    shape_maps = generate_shape_maps(compound)

    print('>> generating framework and sampling...')
    framework = generate_mesh_framework(compound, shape_maps)
    face_meshes = sample_edges_in_framework(framework, shape_maps)

    return MeshkD(path, face_meshes)

def noDoublyLoopInsertions(mesh):
    def checkWire(vertices, wire_mesh):
        ok = True
        for i in range(0, len(wire_mesh)):
            i_index = wire_mesh[i]
            svi = vertices[i_index]
            for j in range(i+1, len(wire_mesh)):
                j_index = wire_mesh[j]
                svj = vertices[j_index]
                if svi == svj:
                    # print('vertex', svi, 'occures both at', i, 'and', j, end='')
                    # print(' (i_index:', i_index, ', j_index:', j_index, ')', sep='')
                    ok = False
        return ok
    ok = True

    for i in range(mesh.number_of_faces()):
        # print('face', i)
        vertices, wire_meshes, _, _ = mesh.face_meshes[i]

        for wire_mesh in wire_meshes:
            if not checkWire(vertices, wire_mesh):
                ok = False

    return ok

if __name__ == '__main__':
    mesh = sample(INPUT_PATH)
    print('no doubly vertices per loop:', noDoublyLoopInsertions(mesh))
    write_to_file(mesh, OUTPUT_DIR)
