""" Chew93_Surface mesher for STEP files using sampler.py w/ OpenCASCADE
Contains functionality to generate a 2-dimensional mesh from a given path to a
STEP file.

To gain a foundation for the meshing process, sampler.py is used to generate a
1-dimensional mesh from the STEP file which is than iterated upon. Here, the
first step is to compute constrained Delaunay Triangulations w/ pytriangle and
then apply my implementation of Chew's Surface Delaunay Refinement algorithm.
"""
from gerobust.predicates import clockwise, counter_clockwise
import numpy as np
import triangle as shewchuk_triangle
import openmesh

import paths
import sampler
from meshkD import SuperVertex, MeshkD, write_to_file
from util import *

from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Lin
from OCC.Core.IntCurvesFace import IntCurvesFace_Intersector
from OCC.Core.TopAbs import TopAbs_ON, TopAbs_IN

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
sampler.NUMBER_OF_SAMPLES = 10
sampler.INCLUDE_INNER_WIRES = True
SMALLEST_ANGLE = np.deg2rad(30)
SIZE_THRESHOLD = 0.01#float('inf')

# __main__ config
INPUT_PATH = paths.STEP_SPHERE
OUTPUT_DIR = paths.DIR_TMP

def triangulate_cdt(face_mesh):
    def segments_from_wires(vertices, wire_meshes):
        segments = []

        for wire_index in range(len(wire_meshes)):
            wire_mesh = wire_meshes[wire_index]
            wire_segments = []

            for i in range(len(wire_mesh)):
                i0 = wire_mesh[i]
                i1 = wire_mesh[(i+1) % len(wire_mesh)]

                if wire_index == 0:
                    wire_segments.append((i0, i1)) # keep ccw order
                else:
                    wire_segments.insert(0, (i1, i0)) # change cw to ccw for pytriangle

            segments += wire_segments

        return segments
    def trim_triangulation(triangles, wire_meshes):
        t_index = 0
        while t_index < len(triangles):
            i0, i1, i2 = triangles[t_index]

            for inner_wire_mesh in wire_meshes[1:]:
                if i0 in inner_wire_mesh and\
                   i1 in inner_wire_mesh and\
                   i2 in inner_wire_mesh:
                    del triangles[t_index]
                    t_index -= 1
                    break

            t_index += 1

        return
    vertices, wire_meshes, triangles, _ = face_mesh

    triangles.clear()

    if len(wire_meshes) == 0:
        return

    t = shewchuk_triangle.Triangle()

    points = [(sv.u, sv.v) for sv in vertices]
    boundary_offset = len(wire_meshes[0])
    markersBoundary = [1] * boundary_offset
    markersInner = [0] * (len(points)-boundary_offset)
    markers = markersBoundary + markersInner
    segments = segments_from_wires(vertices, wire_meshes)

    t.set_points(points, markers=markers)
    t.set_segments(segments)

    t.triangulate(mode='pzQ')

    for triNode in t.get_triangles():
        ([i0, i1, i2], neighbors, attri) = triNode
        triangles.append((i0, i1, i2))

    trim_triangulation(triangles, wire_meshes)

    return

def parse_into_openmesh(face_mesh):
    vertices, _, triangles, _ = face_mesh

    omesh = openmesh.TriMesh()

    vertex_handles = []
    for sv in vertices:
        vertex_handles.append(omesh.add_vertex(sv.XYZ_vec3()))

    for triangle in triangles:
        i0, i1, i2 = triangle
        vh0 = vertex_handles[i0]
        vh1 = vertex_handles[i1]
        vh2 = vertex_handles[i2]
        omesh.add_face(vh0, vh1, vh2)

    return omesh

def flip_until_scdt(omesh, vertices):
    def flip_one_non_scdt_edge(omesh, vertices):
        def collect_quadrilateral_vertices(omesh, eh):
            heh0 = omesh.halfedge_handle(eh, 0)
            heh1 = omesh.halfedge_handle(eh, 1)

            heh_tmp = heh0
            vh0 = omesh.from_vertex_handle(heh_tmp)
            heh_tmp = omesh.next_halfedge_handle(heh_tmp)
            vh1 = omesh.from_vertex_handle(heh_tmp)
            heh_tmp = omesh.next_halfedge_handle(heh_tmp)
            vh2 = omesh.from_vertex_handle(heh_tmp)
            assert omesh.to_vertex_handle(heh_tmp) == vh0
            heh_tmp = heh1
            assert omesh.to_vertex_handle(heh_tmp) == vh0
            assert omesh.from_vertex_handle(heh_tmp) == vh1
            heh_tmp = omesh.next_halfedge_handle(heh_tmp)
            vh3 = omesh.to_vertex_handle(heh_tmp)

            return vh0, vh1, vh2, vh3
        def would_flip_triangle_in_2D(omesh, eh, vertices):
            vh0, vh1, vh2, vh3 = vhs

            sv0 = vertices[vh0.idx()]
            sv1 = vertices[vh1.idx()]
            sv2 = vertices[vh2.idx()]
            sv3 = vertices[vh3.idx()]
            assert np.allclose(sv0.XYZ_vec3(), omesh.point(vh0))
            assert np.allclose(sv1.XYZ_vec3(), omesh.point(vh1))
            assert np.allclose(sv2.XYZ_vec3(), omesh.point(vh2))
            assert np.allclose(sv3.XYZ_vec3(), omesh.point(vh3))

            p0 = sv0.UV_vec2()
            p1 = sv1.UV_vec2()
            p2 = sv2.UV_vec2()
            p3 = sv3.UV_vec2()

            # situation before flip
            t0 = (tuple(p0), tuple(p1), tuple(p2))
            t1 = (tuple(p0), tuple(p3), tuple(p1))
            assert counter_clockwise(*t0) == True
            assert counter_clockwise(*t1) == True

            # situation after flip
            t2 = (tuple(p3), tuple(p2), tuple(p0))
            t3 = (tuple(p3), tuple(p1), tuple(p2))

            return counter_clockwise(*t2) != True or counter_clockwise(*t3) != True
        def flip_maximizes_minimum_angle(omesh, vhs):
            def get_min_face_angle(triangle):
                p0, p1, p2 = triangle

                alpha = calculate_angle_in_corner(p0, p1, p2)
                beta = calculate_angle_in_corner(p1, p2, p0)
                gamma = calculate_angle_in_corner(p2, p0, p1)
                assert np.allclose((alpha+beta+gamma), np.pi)

                min_angle = min(alpha, beta, gamma)

                return min_angle
            vh0, vh1, vh2, vh3 = vhs

            p0 = omesh.point(vh0)
            p1 = omesh.point(vh1)
            p2 = omesh.point(vh2)
            p3 = omesh.point(vh3)

            # situation before flip
            t0 = (p0, p1, p2)
            t1 = (p0, p3, p1)
            min_angle_before = min(get_min_face_angle(t0), get_min_face_angle(t1))

            # situation after flip
            t2 = (p3, p2, p0)
            t3 = (p3, p1, p2)
            min_angle_after = min(get_min_face_angle(t2), get_min_face_angle(t3))

            return min_angle_after > min_angle_before
        for eh in omesh.edges():
            if omesh.is_boundary(eh):
                continue

            vhs = collect_quadrilateral_vertices(omesh, eh)

            if would_flip_triangle_in_2D(omesh, vhs, vertices):
                continue

            if flip_maximizes_minimum_angle(omesh, vhs):
                assert omesh.is_flip_ok(eh)
                omesh.flip(eh)
                return True
            else:
                continue

        return False
    while flip_one_non_scdt_edge(omesh, vertices):
        continue

    return

#TODO cog stub currently
def find_largest_failing_triangle(omesh):
    def vertices_of_face_handle(omesh, fh):
        hh = omesh.halfedge_handle(fh)
        vh0 = omesh.from_vertex_handle(hh)
        hh = omesh.next_halfedge_handle(hh)
        vh1 = omesh.from_vertex_handle(hh)
        hh = omesh.next_halfedge_handle(hh)
        vh2 = omesh.from_vertex_handle(hh)

        p0 = omesh.point(vh0)
        p1 = omesh.point(vh1)
        p2 = omesh.point(vh2)

        return p0, p1, p2
    delta = openmesh.FaceHandle(-1)
    delta_size = float('-inf')

    for fh in omesh.faces():
        p0, p1, p2 = vertices_of_face_handle(omesh, fh)
        size = calculate_area(p0, p1, p2)
        if size < SIZE_THRESHOLD:
            continue
        if size > delta_size:
            delta = fh
            delta_size = size

    print(fh.idx(), delta_size)

    return delta

def insert_cog(omesh, vertices, delta):
    hh = omesh.halfedge_handle(delta)
    vh0 = omesh.from_vertex_handle(hh)
    hh = omesh.next_halfedge_handle(hh)
    vh1 = omesh.from_vertex_handle(hh)
    hh = omesh.next_halfedge_handle(hh)
    vh2 = omesh.from_vertex_handle(hh)

    p0 = vertices[vh0.idx()].UV_vec2()
    p1 = vertices[vh1.idx()].UV_vec2()
    p2 = vertices[vh2.idx()].UV_vec2()
    cog = (p0+p1+p2) / 3

    sv_cog = SuperVertex(u=cog[0], v=cog[1])
    sv_cog.face_id = vertices[0].face_id
    sv_cog.face = vertices[0].face
    sv_cog.project_to_XYZ()
    vertices.append(sv_cog)

    vh_cog = omesh.add_vertex(sv_cog.XYZ_vec3())
    openmesh.TriMesh.split(omesh, delta, vh_cog)
    omesh.garbage_collection()

    return

def split_segment(omesh, vertices, eh):
    hh = omesh.halfedge_handle(eh, 0)
    vh0 = omesh.from_vertex_handle(hh)
    vh1 = omesh.to_vertex_handle(hh)

    sv0 = vertices[vh0.idx()]
    sv1 = vertices[vh1.idx()]
    svhw = SuperVertex.compute_halfway_on_shared_edge(sv0, sv1)
    vertices.append(svhw)

    vhhw = omesh.add_vertex(svhw.XYZ_vec3())
    omesh.split_edge(eh, vhhw)
    omesh.garbage_collection()

    return

def are_consistent(omesh, vertices):
    for vh in omesh.vertices():
        if not all(vertices[vh.idx()].XYZ_vec3() == omesh.point(vh)):
            return False

    return True

def parse_triangles_back(omesh, face_mesh):
    vertices, _, triangles, _ = face_mesh

    assert are_consistent(omesh, vertices)

    triangles.clear()
    for fh in omesh.faces():
        heh0 = omesh.halfedge_handle(fh)
        heh1 = omesh.next_halfedge_handle(heh0)
        heh2 = omesh.next_halfedge_handle(heh1)
        vh0 = omesh.from_vertex_handle(heh0)
        vh1 = omesh.from_vertex_handle(heh1)
        vh2 = omesh.from_vertex_handle(heh2)
        triangles.append((vh0.idx(), vh1.idx(), vh2.idx()))

    return

def chew93_Surface(face_mesh):
    vertices, wire_meshes, triangles, _ = face_mesh

    # step 1: compute initial CDT and iteratively transform into SCDT
    triangulate_cdt(face_mesh)
    omesh = parse_into_openmesh(face_mesh)
    flip_until_scdt(omesh, vertices)

    delta = find_largest_failing_triangle(omesh)

    cnt = 0
    while delta.is_valid():
        # if cnt >= 0:
            # break
        print('iteration', cnt)
        cnt += 1
        insert_cog(omesh, vertices, delta)

        #TODO

        flip_until_scdt(omesh, vertices)
        delta = find_largest_failing_triangle(omesh)

    parse_triangles_back(omesh, face_mesh)

    return

def triangulate(path):
    mesh = sampler.sample(path)

    print('>> mesh samples with chew93_Surface')
    for f_index in range(mesh.number_of_faces()):
        print('face', f_index+1, 'of', mesh.number_of_faces(), '. . .')
        face_mesh = mesh.face_meshes[f_index]
        chew93_Surface(face_mesh)

    mesh.reset_bounding_boxes()

    return mesh

if __name__ == '__main__':
    TMP2 = False
    for TMP in range(1):
        mesh = triangulate(INPUT_PATH)
        TMP2 = False
        write_to_file(mesh, OUTPUT_DIR)
