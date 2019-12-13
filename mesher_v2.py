""" Chew93_Surface mesher for STEP files using sampler.py w/ OpenCASCADE
Contains functionality to generate a 2-dimensional mesh from a given path to a
STEP file.

To gain a foundation for the meshing process, sampler.py is used to generate a
1-dimensional mesh from the STEP file which is than iterated upon. Here, the
first step is to compute constrained Delaunay Triangulations w/ pytriangle and
then apply my implementation of Chew's Surface Delaunay Refinement algorithm.
"""
from gerobust.predicates import clockwise, counter_clockwise
from gerobust import wrapper
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
SIZE_THRESHOLD = float('inf')
MAX_ITERATIONS = 0 # -1 for unlimited

# __main__ config
INPUT_PATH = paths.STEP_BOX
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

def collect_quadrilateral_vertices(omesh, eh):
    hh0 = omesh.halfedge_handle(eh, 0)
    hh1 = omesh.halfedge_handle(eh, 1)

    heh_tmp = hh0
    vh0 = omesh.from_vertex_handle(heh_tmp)
    heh_tmp = omesh.next_halfedge_handle(heh_tmp)
    vh1 = omesh.from_vertex_handle(heh_tmp)
    heh_tmp = omesh.next_halfedge_handle(heh_tmp)
    vh2 = omesh.from_vertex_handle(heh_tmp)
    assert omesh.to_vertex_handle(heh_tmp) == vh0
    heh_tmp = hh1
    assert omesh.to_vertex_handle(heh_tmp) == vh0
    assert omesh.from_vertex_handle(heh_tmp) == vh1
    heh_tmp = omesh.next_halfedge_handle(heh_tmp)
    vh3 = omesh.to_vertex_handle(heh_tmp)

    return vh0, vh1, vh2, vh3

def collect_triangle_vertices(omesh, fh):
    hh = omesh.halfedge_handle(fh)

    vh0 = omesh.from_vertex_handle(hh)
    hh = omesh.next_halfedge_handle(hh)
    vh1 = omesh.from_vertex_handle(hh)
    hh = omesh.next_halfedge_handle(hh)
    vh2 = omesh.from_vertex_handle(hh)
    assert omesh.to_vertex_handle(hh) == vh0

    return vh0, vh1, vh2

def flip_until_scdt(omesh, vertices):
    def flip_one_non_scdt_edge(omesh, vertices):
        def would_flip_triangle_in_2D(omesh, vhs, vertices):
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

#TODO change largest area to largest circumradius
def find_largest_failing_triangle(omesh):
    def shape_test(p0, p1, p2):
        alpha = calculate_angle_in_corner(p0, p1, p2)
        beta = calculate_angle_in_corner(p1, p2, p0)
        gamma = calculate_angle_in_corner(p2, p0, p1)

        return min(alpha, beta, gamma) > SMALLEST_ANGLE
    def size_test(p0, p1, p2):
        t_size = calculate_area(p0, p1, p2)

        return t_size <= SIZE_THRESHOLD, t_size
    delta = openmesh.FaceHandle(-1)
    delta_size = float('-inf')

    for fh in omesh.faces():
        vh0, vh1, vh2 = collect_triangle_vertices(omesh, fh)
        p0 = omesh.point(vh0)
        p1 = omesh.point(vh1)
        p2 = omesh.point(vh2)

        well_shaped = shape_test(p0, p1, p2)
        well_sized, t_size = size_test(p0, p1, p2)

        if well_shaped and well_sized:
            continue

        if t_size > delta_size:
            delta_size = t_size
            delta = fh
        else:
            assert delta.idx() >= 0

    return delta

def calculate_surface_circumcenter(omesh, vertices, delta):
    def barycentric_circumcenter(pa, pb, pc):
        a = np.linalg.norm(pc-pb)
        b = np.linalg.norm(pa-pc)
        c = np.linalg.norm(pb-pa)

        a_squared = a*a
        b_squared = b*b
        c_squared = c*c
        bx = a_squared * (b_squared + c_squared - a_squared)
        by = b_squared * (c_squared + a_squared - b_squared)
        bz = c_squared * (a_squared + b_squared - c_squared)

        return np.array((bx, by, bz))
    vh0, vh1, vh2 = collect_triangle_vertices(omesh, delta)

    sv0 = vertices[vh0.idx()]
    sv1 = vertices[vh1.idx()]
    sv2 = vertices[vh2.idx()]

    bx, by, bz = barycentric_circumcenter(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())

    u, v    = (bx*sv0.UV_vec2()  + by*sv1.UV_vec2()  + bz*sv2.UV_vec2())  / (bx+by+bz)
    # x, y, z = (bx*sv0.XYZ_vec3() + by*sv1.XYZ_vec3() + bz*sv2.XYZ_vec3()) / (bx+by+bz)

    # scc = SuperVertex(x, y, z, u, v)
    scc = SuperVertex(u=u, v=v)
    scc.face_id = sv0.face_id
    scc.face = sv0.face
    scc.project_to_XYZ()

    return scc

#TODO
def travel(omesh, vertices, delta, scc):
    segment = openmesh.EdgeHandle(-1)
    triangle = openmesh.FaceHandle(-1)

    #TODO

    return segment, triangle

def insert_inner_vertex(omesh, vertices, fh, sv_new):
    # vertex insert (vertex list and omesh)
    new_index = len(vertices)
    vertices.append(sv_new)

    vh_new = omesh.add_vertex(sv_new.XYZ_vec3())
    assert vh_new.idx() == new_index

    # integrating into triangulation of omesh
    hh0 = omesh.halfedge_handle(fh)
    hh1 = omesh.next_halfedge_handle(hh0)
    hh2 = omesh.next_halfedge_handle(hh1)
    vh0 = omesh.from_vertex_handle(hh0)
    vh1 = omesh.from_vertex_handle(hh1)
    vh2 = omesh.from_vertex_handle(hh2)

    p0 = vertices[vh0.idx()].UV_vec2()
    p1 = vertices[vh1.idx()].UV_vec2()
    p2 = vertices[vh2.idx()].UV_vec2()
    p_new = sv_new.UV_vec2()

    orient0 = wrapper.orientation_fast(tuple(p0), tuple(p1), tuple(p_new))
    orient1 = wrapper.orientation_fast(tuple(p1), tuple(p2), tuple(p_new))
    orient2 = wrapper.orientation_fast(tuple(p2), tuple(p0), tuple(p_new))

    if np.allclose(orient0, 0.0): # p0, p1 and p_new are colinear
        assert orient1 > 0.0
        assert orient2 > 0.0
        omesh.split_edge(omesh.edge_handle(hh0), vh_new)
    elif np.allclose(orient1, 0.0): # p1, p2 and p_new are colinear
        assert orient2 > 0.0
        omesh.split_edge(omesh.edge_handle(hh1), vh_new)
    elif np.allclose(orient2, 0.0): # p2, p0 and p_new are colinear
        omesh.split_edge(omesh.edge_handle(hh2), vh_new)
    else:
        assert orient0 > 0.0
        assert orient1 > 0.0
        assert orient2 > 0.0
        openmesh.TriMesh.split(omesh, fh, vh_new)

    omesh.garbage_collection()

    return

#TODO
def remove_inner_vertex(omesh, vertices, vh):
    assert vertices[vh.idx()].edges_with_p is None
    assert not omesh.is_boundary(vh)

    #TODO

    return

#TODO
def split_segment(omesh, vertices, eh):
    def remove_encroaching_vertices(omesh, vertices, vh, h):
        def remove_one_encroaching_vertex(omesh, vertices, phw, h):
            for vh in omesh.vertices():
                pcurr = np.array((omesh.point(vh)))

                #TODO consider inf dist if segment in line of sight
                if np.linalg.norm(phw - pcurr) < h:
                    remove_vertex(omesh, vertices, vh)
                    return True

            return False
        phw = np.array((omesh.point(vhhw)))

        while remove_one_encroaching_vertex(omesh, phw, h):
            continue

        return
    assert omesh.is_boundary(eh)

    hh = omesh.halfedge_handle(eh, 0)
    vh0 = omesh.from_vertex_handle(hh)
    vh1 = omesh.to_vertex_handle(hh)
    sv0 = vertices[vh0.idx()]
    sv1 = vertices[vh1.idx()]

    # calculate halfway vertex
    svhw = SuperVertex.compute_halfway_on_shared_edge(sv0, sv1)
    vertices.append(svhw)

    # split segment
    vhhw = omesh.add_vertex(svhw.XYZ_vec3())
    omesh.split_edge(eh, vhhw)
    omesh.garbage_collection()

    # remove encroaching vertices
    h = np.linalg.norm(svhw.XYZ_vec3(), sv0.XYZ_vec3())
    remove_encroaching_vertices(omesh, vertices, vhhw, h)

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

    # step 2+3: find largest triangle that fails shape ans size criteria
    delta = find_largest_failing_triangle(omesh)

    #TODO remove
    scc = calculate_surface_circumcenter(omesh, vertices, delta)
    insert_inner_vertex(omesh, vertices, delta, scc)
    flip_until_scdt(omesh, vertices)
    #TODO remove

    iter_counter = 0
    while delta.is_valid() and iter_counter != MAX_ITERATIONS:
        print('iteration', iter_counter)

        # step 4: travel from any triangle vertex to c and return hit segment id
        scc = calculate_surface_circumcenter(omesh, vertices, delta)
        segment, triangle = travel(omesh, vertices, delta, scc)

        if triangle.is_valid():
            # step 5: no segment was hit, insert scc into triangle
            insert_inner_vertex(omesh, vertices, triangle, scc)
        else:
            # step 6: segment was hit, split
            assert segment.is_valid()
            split_segment(omesh, vertices, segment)

        # after vertex insertion and possible deletion, restore SCDT criteria
        flip_until_scdt(omesh, vertices)

        # update for next iteration
        delta = find_largest_failing_triangle(omesh)
        iter_counter += 1

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
