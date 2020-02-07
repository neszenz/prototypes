""" Chew93_Surface mesher for STEP files using sampler.py w/ OpenCASCADE
Contains functionality to generate a 2-dimensional mesh from a given path to a
STEP file.

To gain a foundation for the meshing process, sampler.py is used to generate a
1-dimensional mesh from the STEP file which is than iterated upon. Here, the
first step is to compute constrained Delaunay Triangulations w/ pytriangle and
then apply my implementation of Chew's Surface Delaunay Refinement algorithm.
"""
import copy
import ctypes
from gerobust.predicates import clockwise, counter_clockwise
from gerobust import wrapper
import numpy as np
import triangle as shewchuk_triangle
import openmesh as om
import sys

import paths
import sampler
from meshkD import SuperVertex, MeshkD, write_to_file
from util import *

from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Lin
from OCC.Core.IntCurvesFace import IntCurvesFace_Intersector
from OCC.Core.TopAbs import TopAbs_ON, TopAbs_IN
from OCC.Core.TopAbs import TopAbs_REVERSED

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
sampler.NUMBER_OF_SAMPLES = 10
sampler.MIN_NUMBER_OF_SAMPLES = 3
sampler.INCLUDE_INNER_WIRES = True
sampler.SIMPLIFY_LINEAR_EDGES = False

MAX_ITERATIONS = -1 # -1 for unlimited

# shape and size test options
SMALLEST_ANGLE = np.deg2rad(30)
USE_SIZE_TEST = False
DISTANCE_THRESHOLD = 0.05

PRIORITIZE_AREA         = 0
PRIORITIZE_CIRCUMRADIUS = 1
PRIORITIZE_DISTANCE     = 2
PRIORITY_FACTOR = PRIORITIZE_CIRCUMRADIUS

# options to avoid 'ping-pong encroachment'
SKIP_ALL_DOMAIN_CORNERS = False # skip all triangles in domain corners
SKIP_PPE_DOMAIN_CORNERS = True # skip those with angle smaller than threshold
PPE_THRESHOLD = SMALLEST_ANGLE#np.deg2rad(60)

# __main__ config
INPUT_PATH = paths.custom(2)
OUTPUT_DIR = paths.TMP_DIR

DEBUG_VERTICES = []

def triangulate_dt(vertices):
    triangles = []

    t = shewchuk_triangle.Triangle()

    points = [(sv.u, sv.v) for sv in vertices]
    markers = [1] * len(points)
    segments = [(i, (i+1)%len(vertices)) for i in range(len(vertices))]

    t.set_points(points, markers=markers)
    t.set_segments(segments)

    t.triangulate(mode='pzQ')

    for triNode in t.get_triangles():
        ([i0, i1, i2], neighbors, attri) = triNode
        triangles.append((i0, i1, i2))

    return triangles

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

def set_supervertex_property(omesh, vh, sv):
    omesh.set_vertex_property('supervertex', vh, id(sv))
    return

def sv_from_vh(omesh, vh):
    int_id = int(omesh.vertex_property_array('supervertex')[vh.idx()])
    return ctypes.cast(int_id, ctypes.py_object).value

def vh_from_sv(omesh, sv):
    supervertices = omesh.vertex_property_array('supervertex')
    result = np.where(supervertices == id(sv))
    assert np.shape(result) == (1, 1)
    return om.VertexHandle(result[0][0])

def parse_into_openmesh(face_mesh):
    vertices, _, triangles, _ = face_mesh

    omesh = om.TriMesh()

    vertex_handles = []
    for sv in vertices:
        vh = omesh.add_vertex(sv.XYZ_vec3())
        set_supervertex_property(omesh, vh, sv)
        vertex_handles.append(vh)

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

def collect_triangle_vertex_handles(omesh, fh):
    hh = omesh.halfedge_handle(fh)

    vh0 = omesh.from_vertex_handle(hh)
    hh = omesh.next_halfedge_handle(hh)
    vh1 = omesh.from_vertex_handle(hh)
    hh = omesh.next_halfedge_handle(hh)
    vh2 = omesh.from_vertex_handle(hh)
    assert omesh.to_vertex_handle(hh) == vh0

    return vh0, vh1, vh2

def collect_triangle_supervertices(omesh, fh):
    vh0, vh1, vh2 = collect_triangle_vertex_handles(omesh, fh)

    sv0 = sv_from_vh(omesh, vh0)
    sv1 = sv_from_vh(omesh, vh1)
    sv2 = sv_from_vh(omesh, vh2)

    return sv0, sv1, sv2

# creates scdt criteria by edge flips via repeated linear search
def flip_until_scdt(omesh):
    def flip_one_non_scdt_edge(omesh):
        def would_flip_triangle_in_2D(omesh, vhs):
            vh0, vh1, vh2, vh3 = vhs

            sv0 = sv_from_vh(omesh, vh0)
            sv1 = sv_from_vh(omesh, vh1)
            sv2 = sv_from_vh(omesh, vh2)
            sv3 = sv_from_vh(omesh, vh3)
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

            if would_flip_triangle_in_2D(omesh, vhs):
                continue

            if flip_maximizes_minimum_angle(omesh, vhs):
                assert omesh.is_flip_ok(eh)
                omesh.flip(eh)
                return True

        return False
    while flip_one_non_scdt_edge(omesh):
        continue

    return

# restores scdt criteria by edge flips via depth-first-search (dfs)
def restore_scdt(omesh, vh_start):
    def create_backlog_from_neighbors(omesh, vh):
        backlog = []

        for vvh in omesh.vv(vh):
            backlog.append(vvh)

        return backlog
    def process_vertex(omesh, vh):
        def would_flip_triangle_in_2D(omesh, vhs):
            vh0, vh1, vh2, vh3 = vhs

            sv0 = sv_from_vh(omesh, vh0)
            sv1 = sv_from_vh(omesh, vh1)
            sv2 = sv_from_vh(omesh, vh2)
            sv3 = sv_from_vh(omesh, vh3)
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
        backlog_extension = []

        incident_edges = []
        for hh in omesh.voh(vh):
            incident_edges.append(omesh.edge_handle(hh))

        for eh in incident_edges:
            if omesh.is_boundary(eh):
                continue

            vhs = collect_quadrilateral_vertices(omesh, eh)

            if would_flip_triangle_in_2D(omesh, vhs):
                continue

            if flip_maximizes_minimum_angle(omesh, vhs):
                assert omesh.is_flip_ok(eh)
                backlog_extension += vhs
                omesh.flip(eh)

        return backlog_extension
    backlog = create_backlog_from_neighbors(omesh, vh_start)

    while len(backlog) > 0:
        vh = backlog.pop()

        backlog_extension = process_vertex(omesh, vh)

        backlog += backlog_extension

    return

def find_largest_failing_triangle(omesh):
    def is_domain_corner(omesh, sv0, sv1, sv2):
        if not sv0.edges_with_p is None:
            if len(sv0.edges_with_p) == 2:
                return True
        if not sv1.edges_with_p is None:
            if len(sv1.edges_with_p) == 2:
                return True
        if not sv2.edges_with_p is None:
            if len(sv2.edges_with_p) == 2:
                return True

        return False
    def skip_to_avoid_ping_pong_encroachment(omesh, fh):
        def in_danger_of_encroachment(omesh, hh0, hh1):
            assert hh1 == omesh.next_halfedge_handle(hh0)
            eh0 = omesh.edge_handle(hh0)
            eh1 = omesh.edge_handle(hh1)

            if omesh.is_boundary(eh0) and omesh.is_boundary(eh1):
                p0 = omesh.point(omesh.from_vertex_handle(hh0))
                p1 = omesh.point(omesh.to_vertex_handle(hh0))
                p2 = omesh.point(omesh.to_vertex_handle(hh1))
                angle = calculate_angle_between_vectors(normalize(p0-p1), normalize(p2-p1))

                return angle < PPE_THRESHOLD

            return False
        hh0 = omesh.halfedge_handle(fh)
        hh1 = omesh.next_halfedge_handle(hh0)
        hh2 = omesh.next_halfedge_handle(hh1)

        if in_danger_of_encroachment(omesh, hh0, hh1) or\
           in_danger_of_encroachment(omesh, hh1, hh2) or\
           in_danger_of_encroachment(omesh, hh2, hh0):
            return True

        return False
    def shape_test(sv0, sv1, sv2):
        p0 = sv0.XYZ_vec3()
        p1 = sv1.XYZ_vec3()
        p2 = sv2.XYZ_vec3()

        alpha = calculate_angle_in_corner(p0, p1, p2)
        beta = calculate_angle_in_corner(p1, p2, p0)
        gamma = calculate_angle_in_corner(p2, p0, p1)

        return min(alpha, beta, gamma) > SMALLEST_ANGLE
    def size_test(sv0, sv1, sv2):
        def approximation_distance(sv0, sv1, sv2):
            cog_2d = (sv0.UV_vec2() + sv1.UV_vec2() + sv2.UV_vec2()) / 3
            sv_cog = SuperVertex(u=cog_2d[0], v=cog_2d[1])
            sv_cog.set_same_face_as(sv0)
            sv_cog.project_to_XYZ()

            cog_surface = sv_cog.XYZ_vec3()
            cog_3d = (sv0.XYZ_vec3() + sv1.XYZ_vec3() + sv2.XYZ_vec3()) / 3
            distance = np.linalg.norm(cog_surface - cog_3d)

            return distance
        t_distance = approximation_distance(sv0, sv1, sv2)

        if USE_SIZE_TEST:
            return t_distance < DISTANCE_THRESHOLD, t_distance
        else:
            return True, t_distance
    delta = om.FaceHandle(-1)
    delta_size = float('-inf')

    for fh in omesh.faces():
        sv0, sv1, sv2 = collect_triangle_supervertices(omesh, fh)
        p0 = sv0.XYZ_vec3()
        p1 = sv1.XYZ_vec3()
        p2 = sv2.XYZ_vec3()

        if SKIP_ALL_DOMAIN_CORNERS and is_domain_corner(omesh, sv0, sv1, sv2):
            continue
        if SKIP_PPE_DOMAIN_CORNERS and skip_to_avoid_ping_pong_encroachment(omesh, fh):
            continue

        well_shaped = shape_test(sv0, sv1, sv2)
        well_sized, t_distance = size_test(sv0, sv1, sv2)

        if well_shaped and well_sized:
            continue

        if PRIORITY_FACTOR == PRIORITIZE_AREA:
            t_size = calculate_area(p0, p1, p2)
        elif PRIORITY_FACTOR == PRIORITIZE_CIRCUMRADIUS:
            t_size = calculate_circumradius_v2(p0, p1, p2)
        elif PRIORITY_FACTOR == PRIORITIZE_DISTANCE:
            t_size = t_distance
        else:
            raise Exception('PRIORITY_FACTOR unknown')

        if t_size > delta_size:
            delta_size = t_size
            delta = fh

    return delta

# for triangle p0-p1-p2 returns orientation test results of all three edges with px
def orientations(p0, p1, p2, px):
    ori0 = wrapper.orientation_fast(tuple(p0), tuple(p1), tuple(px))
    ori1 = wrapper.orientation_fast(tuple(p1), tuple(p2), tuple(px))
    ori2 = wrapper.orientation_fast(tuple(p2), tuple(p0), tuple(px))

    return ori0, ori1, ori2

def calculate_refinement(omesh, delta, meta_block):
    def calculate_circumcenter_3d(A, B, C):
        # by 'Oscar Lanzi III' from
        # https://sci.math.narkive.com/nJMeroLe/circumcenter-of-a-3d-triangle
        a = np.linalg.norm(B-C)
        b = np.linalg.norm(A-C)
        c = np.linalg.norm(A-B)
        a_square = a**2
        b_square = b**2
        c_square = c**2
        A_comp = A*(a_square*(b_square+c_square-a_square))
        B_comp = B*(b_square*(a_square+c_square-b_square))
        C_comp = C*(c_square*(a_square+b_square-c_square))
        divisor = 2*(a_square*b_square+a_square*c_square+b_square*c_square)-(a**4+b**4+c**4)
        circumcenter = (A_comp+B_comp+C_comp) / divisor

        return circumcenter
    def calculate_circumcenter_2d_candidates(c_3d, n, face):
        def array_from_inter(inter, i):
            u = inter.UParameter(i)
            v = inter.VParameter(i)
            return np.array((u, v))
        center_line = gp_Lin(gp_Pnt(c_3d[0], c_3d[1], c_3d[2]), gp_Dir(n[0], n[1], n[2]))

        inter = IntCurvesFace_Intersector(face, 0.0001)
        inter.Perform(center_line, -float('inf'), float('inf'))

        if inter.IsDone():
            occ_offset = 1 # occ indices start at 1
            c_2d_candidates = [array_from_inter(inter, i+occ_offset) for i in range(inter.NbPnt())]
            if face.Orientation() == TopAbs_REVERSED:
                for c_2d in c_2d_candidates:
                    # reverse u axis for reversed faces
                    c_2d[0] = reverse_u(c_2d[0], face)
            return c_2d_candidates
        else:
            raise Exception('calculate_surface_circumcenter() error - intersector not done')
    def find_point_inside(omesh, fh, points):
        def lies_inside(omesh, fh, p):
            sv0, sv1, sv2 = collect_triangle_supervertices(omesh, fh)
            p0 = tuple(sv0.UV_vec2())
            p1 = tuple(sv1.UV_vec2())
            p2 = tuple(sv2.UV_vec2())

            ori0, ori1, ori2 = orientations(p0, p1, p2, p)

            if np.allclose(ori0, 0.0) or np.allclose(ori1, 0.0) or np.allclose(ori2, 0.0):
                return True

            if ori0 > 0.0 and ori1 > 0.0 and ori2 > 0.0:
                return True

            return False
        for i in range(len(points)):
            if lies_inside(omesh, fh, points[i]):
                return i

        return -1
    def scc_from_c_2d(c_2d, other_sv):
        scc = SuperVertex(u=c_2d[0], v=c_2d[1])
        scc.face = other_sv.face
        scc.face_id = other_sv.face_id
        scc.project_to_XYZ()

        return scc
    def travel(omesh, delta, hh, c_3d, c_2d_candidates, normal, meta_block):
        def traveled_too_far(omesh, hh, p_orig, ray_dir):
            vh_opposite = omesh.to_vertex_handle(hh)
            p_opposite = sv_from_vh(omesh, vh_opposite).XYZ_vec3()
            p_opposite -= p_orig # move system to origin

            p_projected = np.abs(np.dot(p_opposite, ray_dir))
            if p_projected > 1.0:
                print('traveled too far:', p_projected)

            return p_projected > 1.0
        def halfedge_crossed_by_ray_shadow(omesh, hh, ray_ori, ray_dir, normal):
            # print('halfedge_crossed_by_ray_shadow: ', end='')
            assert np.allclose(np.linalg.norm(ray_dir), 1.0)
            assert np.allclose(np.linalg.norm(normal), 1.0)
            p_from = sv_from_vh(omesh, omesh.from_vertex_handle(hh)).XYZ_vec3()
            p_to = sv_from_vh(omesh, omesh.to_vertex_handle(hh)).XYZ_vec3()

            v_from = shortest_vector_between_two_lines(p_from, normal, ray_ori, ray_dir)
            v_to = shortest_vector_between_two_lines(p_to, normal, ray_ori, ray_dir)

            dp = np.dot(v_from, v_to)
            # print(dp < 0)

            return dp < 0
        meta_block[MeshkD.NT_INVK] += 1

        svh0 = sv_from_vh(omesh, omesh.from_vertex_handle(hh))
        svh1 = sv_from_vh(omesh, omesh.to_vertex_handle(hh))
        sv_orig = SuperVertex.compute_halfway(svh0, svh1)
        p_orig = sv_orig.XYZ_vec3()
        ray_dir = normalize(c_3d - p_orig)

        # halt if we encounter a boundary edge (split case)
        while not omesh.is_boundary(omesh.edge_handle(hh)):
            meta_block[MeshkD.NT_LOOP] += 1
            # get inside neighboring triangle
            hh = omesh.opposite_halfedge_handle(hh)
            fh = omesh.face_handle(hh)

            # halt if surface circumcenter is found to be insides fh (insert case)
            inside_index = find_point_inside(omesh, fh, c_2d_candidates)
            if inside_index >= 0:
                return fh, scc_from_c_2d(c_2d_candidates[inside_index], sv_orig)

            hh = omesh.next_halfedge_handle(hh) # go to first of the opposite halfedges

            if traveled_too_far(omesh, hh, p_orig, ray_dir):
                pass

            # which of the remaining edges is crossed to travel further in that direction
            if not halfedge_crossed_by_ray_shadow(omesh, hh, p_orig, ray_dir, normal):
                # shadow apparently crossed the other edge
                hh = omesh.next_halfedge_handle(hh)
                assert halfedge_crossed_by_ray_shadow(omesh, hh, p_orig, ray_dir, normal)
            meta_block[MeshkD.NT_SHAD] += 1

        return omesh.edge_handle(hh), None
    hh = omesh.halfedge_handle(delta)
    vh0 = omesh.from_vertex_handle(hh)
    vh1 = omesh.to_vertex_handle(hh)
    hh = omesh.next_halfedge_handle(hh) # hh is now opposite of vh0
    vh2 = omesh.to_vertex_handle(hh)

    sv0 = sv_from_vh(omesh, vh0)
    sv1 = sv_from_vh(omesh, vh1)
    sv2 = sv_from_vh(omesh, vh2)

    normal = calculate_normal_normalized(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())
    c_3d = calculate_circumcenter_3d(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())
    #TODO remove
    # if iter_counter == 20:
    #     sv_c_3d = SuperVertex(*c_3d, same_as=sv0)
    #     sv_c_3d.project_to_UV()
    #     DEBUG_VERTICES.append(sv_c_3d)
    #TODO remove
    c_2d_candidates = calculate_circumcenter_2d_candidates(c_3d, normal, sv0.face)
    #TODO remove
    # if iter_counter == 20:
    #     for (u, v) in c_2d_candidates:
    #         sv_c_2d = SuperVertex(u=u, v=v, same_as=sv0)
    #         # DEBUG_VERTICES.append(sv_c_2d)
    #TODO remove

    # test if surface circumcenter is in delta
    inside_index = find_point_inside(omesh, delta, c_2d_candidates)
    if inside_index >= 0:
        return delta, scc_from_c_2d(c_2d_candidates[inside_index], sv0)

    # it is outside -> one angle needs to be greater than the other
    angle0 = calculate_angle_in_corner(sv2.XYZ_vec3(), sv0.XYZ_vec3(), sv1.XYZ_vec3())
    angle1 = calculate_angle_in_corner(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3())
    angle2 = calculate_angle_in_corner(sv1.XYZ_vec3(), sv2.XYZ_vec3(), sv0.XYZ_vec3())

    if angle0 > angle1 and angle0 > angle2:
        hh_selected = hh
    elif angle1 > angle0 and angle1 > angle2:
        hh_selected = omesh.next_halfedge_handle(hh)
    elif angle2 > angle0 and angle2 > angle1:
        hh_selected = omesh.prev_halfedge_handle(hh)
    else:
        raise Exception('calculate_refinement() error - neither inside nor outside')

    handle, scc = travel(omesh, delta, hh_selected, c_3d, c_2d_candidates, normal, meta_block)

    return handle, scc

def insert_inner_vertex(omesh, vertices, fh, sv_new, meta_block):
    # vertex insert (vertex list and omesh)
    vertices.append(sv_new)

    vh_new = omesh.add_vertex(sv_new.XYZ_vec3())
    set_supervertex_property(omesh, vh_new, sv_new)

    # integrating into triangulation of omesh
    hh0 = omesh.halfedge_handle(fh)
    hh1 = omesh.next_halfedge_handle(hh0)
    hh2 = omesh.next_halfedge_handle(hh1)
    vh0 = omesh.from_vertex_handle(hh0)
    vh1 = omesh.from_vertex_handle(hh1)
    vh2 = omesh.from_vertex_handle(hh2)

    p0 = sv_from_vh(omesh, vh0).UV_vec2()
    p1 = sv_from_vh(omesh, vh1).UV_vec2()
    p2 = sv_from_vh(omesh, vh2).UV_vec2()
    p_new = sv_new.UV_vec2()

    ori0, ori1, ori2 = orientations(p0, p1, p2, p_new)

    if np.allclose(ori0, 0.0): # p0, p1 and p_new are colinear
        assert ori1 > 0.0
        assert ori2 > 0.0
        eh = omesh.edge_handle(hh0)
        omesh.split_edge(eh, vh_new)
    elif np.allclose(ori1, 0.0): # p1, p2 and p_new are colinear
        assert ori2 > 0.0
        eh = omesh.edge_handle(hh1)
        omesh.split_edge(eh, vh_new)
    elif np.allclose(ori2, 0.0): # p2, p0 and p_new are colinear
        eh = omesh.edge_handle(hh2)
        omesh.split_edge(eh, vh_new)
    else:
        assert ori0 > 0.0
        assert ori1 > 0.0
        assert ori2 > 0.0
        om.TriMesh.split(omesh, fh, vh_new)

    omesh.garbage_collection()
    meta_block[MeshkD.NV_INST] += 1

    assert all(omesh.point(vh_new) == sv_new.XYZ_vec3())
    return vh_new

def remove_inner_vertex(omesh, vertices, vh):
    def remove_incident_faces(omesh, vh):
        incident_faces = []
        for fh in omesh.vf(vh):
            incident_faces.append(fh)

        while len(incident_faces) > 0:
            fh = incident_faces.pop()
            omesh.delete_face(fh, False)

        omesh.garbage_collection()

        return
    # print('remove_inner_vertex()')
    assert sv_from_vh(omesh, vh).edges_with_p is None
    assert not omesh.is_boundary(vh)

    # collect 1-ring neigbors for computing new triangulation
    neighbor_handles = []
    neighbor_vertices = []
    for vvh in omesh.vv(vh):
        neighbor_handles.append(vvh)
        neighbor_vertices.append(sv_from_vh(omesh, vvh))

    triangles = triangulate_dt(neighbor_vertices)

    # replace old by new triangulation
    remove_incident_faces(omesh, vh)
    for (i0, i1, i2) in triangles:
        vh0 = neighbor_handles[i0]
        vh1 = neighbor_handles[i1]
        vh2 = neighbor_handles[i2]
        omesh.add_face(vh0, vh1, vh2)

    # only now remove vertex to avoid vertex handles to become invalid
    vertices.remove(sv_from_vh(omesh, vh))
    omesh.delete_vertex(vh, False)
    omesh.garbage_collection()

    return

#TODO consider inf dist if segment in line of sight
#TODO imporove encroaching vertices removal performance
def split_edge(omesh, vertices, eh, meta_block):
    def remove_encroaching_vertices(omesh, vertices, vh, h):
        def remove_one_encroaching_vertex(omesh, vertices, phw, h):
            for vh in omesh.vertices():
                if omesh.is_boundary(vh):
                    continue

                pcurr = np.array((omesh.point(vh)))

                if np.linalg.norm(phw - pcurr) < h:
                    remove_inner_vertex(omesh, vertices, vh)
                    return True

            return False
        phw = np.array((omesh.point(vh)))

        n_deleted = 0
        while remove_one_encroaching_vertex(omesh, vertices, phw, h):
            n_deleted += 1

        return n_deleted
    hh = omesh.halfedge_handle(eh, 0)
    vh0 = omesh.from_vertex_handle(hh)
    vh1 = omesh.to_vertex_handle(hh)
    sv0 = sv_from_vh(omesh, vh0)
    sv1 = sv_from_vh(omesh, vh1)

    if omesh.is_boundary(eh):
        sv_new = SuperVertex.compute_halfway_on_shared_edge(sv0, sv1)
    else:
        sv_new = SuperVertex.compute_halfway(sv0, sv1)

    vertices.append(sv_new)

    # split segment
    vh_new = omesh.add_vertex(sv_new.XYZ_vec3())
    set_supervertex_property(omesh, vh_new, sv_new)
    omesh.split_edge(eh, vh_new)
    omesh.garbage_collection()
    meta_block[MeshkD.NV_SPLT] += 1

    # remove encroaching vertices
    if omesh.is_boundary(eh):
        h = np.linalg.norm(sv_new.XYZ_vec3() - sv0.XYZ_vec3())
        n_deleted = remove_encroaching_vertices(omesh, vertices, vh_new, h)
        meta_block[MeshkD.NV_DELT] += n_deleted

    return vh_from_sv(omesh, sv_new)

def parse_back(omesh, face_mesh):
    vertices, _, triangles, _ = face_mesh
    assert len(vertices) == len(omesh.vertices())

    # Supervertices are stored only once in the vertices list and openmesh
    # changes index order eventually during vertex deletion. Therefor the list's
    # order needs to be updated without overwriting (to avoid possibly deleting
    # a supervertex)
    number_of_vertices = len(vertices)
    for vh in omesh.vertices():
        vertices.append(sv_from_vh(omesh, vh))
    del vertices[:number_of_vertices]

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

def chew93_Surface(face_mesh, meta_block):
    vertices, wire_meshes, triangles, _ = face_mesh

    # step 1: compute initial CDT and iteratively transform into SCDT
    print('generating initial mesh...', end='')
    sys.stdout.flush()
    triangulate_cdt(face_mesh)
    omesh = parse_into_openmesh(face_mesh)
    flip_until_scdt(omesh)

    # step 2+3: find largest triangle that fails shape ans size criteria
    print('\b\b\b - done!\nrefining mesh...')
    delta = find_largest_failing_triangle(omesh)

    global iter_counter
    iter_counter = 0
    while delta.is_valid() and iter_counter != MAX_ITERATIONS:
        print('iteration', iter_counter)
        #TODO remove
        # if iter_counter == 20:
        #     sv0, sv1, sv2 = collect_triangle_supervertices(omesh, delta)
        #     (x, y, z), (u, v) = calculate_cog(sv0, sv1, sv2)
        #     sv_cog = SuperVertex(x, y, z, u, v, same_as=sv0)
        #     DEBUG_VERTICES.append(sv_cog)
        #TODO remove

        handle, scc = calculate_refinement(omesh, delta, meta_block)

        if isinstance(handle, om.FaceHandle):
            vh_new = insert_inner_vertex(omesh, vertices, handle, scc, meta_block)
        else:
            assert isinstance(handle, om.EdgeHandle)
            vh_new = split_edge(omesh, vertices, handle, meta_block)

        meta_block[MeshkD.NV_REFI] += 1

        # after vertex insertion and possible deletion, restore SCDT criteria
        restore_scdt(omesh, vh_new)

        # update for next iteration
        delta = find_largest_failing_triangle(omesh)
        iter_counter += 1

    parse_back(omesh, face_mesh)
    # vertices += DEBUG_VERTICES

    return

def triangulate(path):
    mesh = sampler.sample(path)

    print('>> mesh samples with chew93_Surface')
    for f_index in range(mesh.number_of_faces()):
        print('face', f_index+1, 'of', mesh.number_of_faces(), '. . .')
        face_mesh = mesh.face_meshes[f_index]
        meta_block = mesh.meta_blocks[f_index]
        chew93_Surface(face_mesh, meta_block)
        print('done after', meta_block[MeshkD.NV_REFI], 'refinements')

    mesh.reset_bounding_boxes()

    return mesh

if __name__ == '__main__':
    for TMP in range(1):
        # MAX_ITERATIONS = TMP
        mesh = triangulate(INPUT_PATH)
        write_to_file(mesh, OUTPUT_DIR)
