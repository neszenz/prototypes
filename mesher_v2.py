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
sampler.SIMPLIFY_LINEAR_EDGES = False
SMALLEST_ANGLE = np.deg2rad(30)
SIZE_THRESHOLD = float('inf')
MAX_ITERATIONS = 12 # -1 for unlimited
USE_TRAVEL_TEST = True

# __main__ config
INPUT_PATH = paths.STEP_SURFFIL
OUTPUT_DIR = paths.DIR_TMP

DEBUG_VERTICES = [] #TODO remove

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
            else:
                continue

        return False
    while flip_one_non_scdt_edge(omesh):
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
    delta = om.FaceHandle(-1)
    delta_size = float('-inf')

    for fh in omesh.faces():
        vh0, vh1, vh2 = collect_triangle_vertex_handles(omesh, fh)
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

# calculate barycentric coordinates of 3-dimensional circumcenter
def calculate_bcc(A, B, C):
    a = np.linalg.norm(C-B)
    b = np.linalg.norm(A-C)
    c = np.linalg.norm(B-A)
    a_squared = a*a
    b_squared = b*b
    c_squared = c*c

    bx = a_squared * (b_squared + c_squared - a_squared)
    by = b_squared * (c_squared + a_squared - b_squared)
    bz = c_squared * (a_squared + b_squared - c_squared)

    return np.array((bx, by, bz))

def cartesian_from_barycentric(p0, p1, p2, bx, by, bz):
    cartesian = (bx*p0 + by*p1 + bz*p2) / (bx+by+bz)
    return np.array(cartesian)

#TODO correct approxcimation by normal convergence
# calculate surface circumcenter based on circumcenter in barycentric coordinates
def scc_from_bcc(sv0, sv1, sv2, bcc):
    u, v    = cartesian_from_barycentric(sv0.UV_vec2(),  sv1.UV_vec2(),  sv2.UV_vec2(),  *bcc)
    # x, y, z = cartesian_from_barycentric(sv0.XYZ_vec3(), sv1.XYZ_vec3(), sv2.XYZ_vec3(), *bcc)

    # scc = SuperVertex(x, y, z, u, v)
    scc = SuperVertex(u=u, v=v)
    scc.face_id = sv0.face_id
    scc.face = sv0.face
    scc.project_to_XYZ()

    return scc

# while the cc in 3d can be calc. exactly, for the 2d ones only possible candidates can be given
def circumcenters_from_bcc(svA, svB, svC, bcc):
    def array_from_inter(inter, i):
        u = inter.UParameter(i)
        v = inter.VParameter(i)
        return np.array((u, v))
    a = svA.XYZ_vec3()
    b = svB.XYZ_vec3()
    c = svC.XYZ_vec3()

    c_3d = cartesian_from_barycentric(svA.XYZ_vec3(), svB.XYZ_vec3(), svC.XYZ_vec3(), *bcc)
    n = normalize(np.cross(b-a, c-a))

    center_line = gp_Lin(gp_Pnt(c_3d[0], c_3d[1], c_3d[2]), gp_Dir(n[0], n[1], n[2]))

    inter = IntCurvesFace_Intersector(svA.face, 0.0001)
    inter.Perform(center_line, -float('inf'), float('inf'))

    if inter.IsDone():
        occ_offset = 1 # occ indices start at 1
        c_2d_candidates =  [array_from_inter(inter, i+occ_offset) for i in range(inter.NbPnt())]
        return c_2d_candidates, c_3d, n
    else:
        raise Exception('calculate_surface_circumcenter() error - intersector not done')

# for triangle p0-p1-p2 returns orientation test results of all three edges with px
def orientations(p0, p1, p2, px):
    ori0 = wrapper.orientation_fast(tuple(p0), tuple(p1), tuple(px))
    ori1 = wrapper.orientation_fast(tuple(p1), tuple(p2), tuple(px))
    ori2 = wrapper.orientation_fast(tuple(p2), tuple(p0), tuple(px))

    return ori0, ori1, ori2

def find_point_inside(omesh, fh, points):
    def lies_inside(omesh, fh, p):
        sv0, sv1, sv2 = collect_triangle_supervertices(omesh, fh)
        p0 = tuple(sv0.UV_vec2())
        p1 = tuple(sv1.UV_vec2())
        p2 = tuple(sv2.UV_vec2())

        ori0, ori1, ori2 = orientations(p0, p1, p2, p)

        if np.allclose(ori0, 0.0) or np.allclose(ori1, 0.0) or np.allclose(ori2, 0.0):
            return False

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

#TODO
def travel(omesh, delta, hh, sv_orig, c_3d, c_2d_candidates, normal):
    def halfedge_intersects_ray_shadow(omesh, hh, p_orig, ray_dir, normal): #TODO test/ fix
        print('halfedge_intersects_ray_shadow()')
        assert np.allclose(np.linalg.norm(ray_dir), 1.0)
        assert np.allclose(np.linalg.norm(normal), 1.0)
        p_from = sv_from_vh(omesh, omesh.from_vertex_handle(hh)).XYZ_vec3()
        p_to = sv_from_vh(omesh, omesh.to_vertex_handle(hh)).XYZ_vec3()

        # to simplify, move ray origin to system origin
        p0 = p_from - p_orig
        p1 = p_to - p_orig

        # run binary search until we converge to machine percision
        while not np.allclose(np.linalg.norm(p1-p0), 0.0):
            ph = p0 + ((p1-p0) * 0.5)
            #TODO remove
            global DEBUG_VERTICES
            x, y, z = ph + p_orig
            DEBUG_VERTICES.append(SuperVertex(x, y, z))
            #TODO remove

            d0 = project_point_onto_normalized_vector(p0, ray_dir)
            d1 = project_point_onto_normalized_vector(p1, ray_dir)
            dh = project_point_onto_normalized_vector(ph, ray_dir)
            #TODO remove
            x, y, z = d0 + p_orig
            DEBUG_VERTICES.append(SuperVertex(x, y, z))
            x, y, z = d1 + p_orig
            DEBUG_VERTICES.append(SuperVertex(x, y, z))
            x, y, z = dh + p_orig
            DEBUG_VERTICES.append(SuperVertex(x, y, z))
            #TODO remove
            dp0 = p0 - d0
            dp1 = p1 - d1
            dph = ph - dh

            # angle0 = calculate_angle_between_vectors(dp0, normal)
            # angle1 = calculate_angle_between_vectors(dp1, normal)
            # angleh = calculate_angle_between_vectors(dph, normal)
            error0h = np.abs(np.linalg.norm(dph) - np.linalg.norm(dp0))
            error1h = np.abs(np.linalg.norm(dph) - np.linalg.norm(dp1))

            # error0h = np.abs(angleh - angle0)
            # error1h = np.abs(angleh - angle1)

            if error0h > error1h:
                p0 = ph
            else:
                p1 = ph

        p_result = p0 + p_orig

        assert not np.allclose(np.linalg.norm(p_from - p_result), 0.0)

        if np.allclose(np.linalg.norm(p_to - p_result), 0.0):
            print('False')
            return False
        else:
            print('True')
            return True
    p_orig = sv_orig.XYZ_vec3()
    ray_dir = normalize(c_3d - p_orig)

    # halt if we encounter a boundary edge (split case)
    while not omesh.is_boundary(omesh.edge_handle(hh)):
        hh = omesh.opposite_halfedge_handle(hh)
        fh = omesh.face_handle(hh)

        # halt if surface circumcenter is found to be insides fh (insert case)
        inside_index = find_point_inside(omesh, fh, c_2d_candidates)
        if inside_index >= 0:
            return fh, scc_from_c_2d(c_2d_candidates[inside_index], sv_orig)

        # which of the remaining edges is crossed to travel further in that direction
        hh = omesh.next_halfedge_handle(hh)
        if not halfedge_intersects_ray_shadow(omesh, hh, p_orig, ray_dir, normal):
            # shadow apparently crossed the other edge
            hh = omesh.next_halfedge_handle(hh)
            assert halfedge_intersects_ray_shadow(omesh, hh, p_orig, ray_dir, normal) #TODO remove after testing

    return omesh.edge_handle(hh), None

def calculate_refinement(omesh, delta):
    hhc = omesh.halfedge_handle(delta)
    hha = omesh.next_halfedge_handle(hhc)
    hhb = omesh.next_halfedge_handle(hha)
    svA = sv_from_vh(omesh, omesh.from_vertex_handle(hhc))
    svB = sv_from_vh(omesh, omesh.from_vertex_handle(hha))
    svC = sv_from_vh(omesh, omesh.from_vertex_handle(hhb))

    bcc = calculate_bcc(svA.XYZ_vec3(), svB.XYZ_vec3(), svC.XYZ_vec3())
    bx, by, bz = bcc

    # based on how barycentric coordinates work, we can first test for simple cases
    if np.allclose(bx, 0.0):
        return omesh.edge_handle(hha), None
    if np.allclose(by, 0.0):
        return omesh.edge_handle(hhb), None
    if np.allclose(bz, 0.0):
        return omesh.edge_handle(hhc), None

    # scc_approx = scc_from_bcc(svA, svB, svC, bcc)
    c_2d_candidates, c_3d, normal = circumcenters_from_bcc(svA, svB, svC, bcc)

    inside_index = find_point_inside(omesh, delta, c_2d_candidates)
    if inside_index >= 0:
        return delta, scc_from_c_2d(c_2d_candidates[inside_index], svA)
    if all(bcc > 0.0):
        for c_2d in c_2d_candidates:
            tmp = SuperVertex(u=c_2d[0], v=c_2d[1])
            tmp.face = svA.face
            tmp.face_id = svA.face_id
            tmp.project_to_XYZ()
            global DEBUG_VERTICES
            DEBUG_VERTICES.append(tmp)
        return delta, scc_from_bcc(svA, svB, svC, bcc)
    assert not all(bcc > 0.0)
    # if all(bcc > 0.0):
        # return delta, scc_approx

    # none of the simple cases are met; we need to invoke the travel algorithm
    if bx < 0.0:
        if USE_TRAVEL_TEST:
            handle, scc = travel(omesh, delta, hha, svA, c_3d, c_2d_candidates, normal)
        else:
            handle, scc = omesh.edge_handle(hha), scc_approx
    elif by < 0.0:
        if USE_TRAVEL_TEST:
            handle, scc = travel(omesh, delta, hhb, svB, c_3d, c_2d_candidates, normal)
        else:
            handle, scc = omesh.edge_handle(hhb), scc_approx
    elif bz < 0.0:
        if USE_TRAVEL_TEST:
            handle, scc = travel(omesh, delta, hhc, svC, c_3d, c_2d_candidates, normal)
        else:
            handle, scc = omesh.edge_handle(hhc), scc_approx
    else:
        raise Exception('calculate_refinement() error - bcc meets none of the cases')

    return handle, scc

def insert_inner_vertex(omesh, vertices, fh, sv_new):
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

    return

def remove_inner_vertex(omesh, vertices, vh):
    assert sv_from_vh(omesh, vh).edges_with_p is None
    assert not omesh.is_boundary(vh)

    neighbors = []
    for vvh in omesh.vv(vh):
        neighbors.append(sv_from_vh(omesh, vvh))

    triangles = triangulate_dt(neighbors)

    vertices.remove(sv_from_vh(omesh, vh))
    omesh.delete_vertex(vh, False)
    omesh.garbage_collection()

    for (i0, i1, i2) in triangles:
        vh0 = vh_from_sv(omesh, neighbors[i0])
        vh1 = vh_from_sv(omesh, neighbors[i1])
        vh2 = vh_from_sv(omesh, neighbors[i2])
        omesh.add_face(vh0, vh1, vh2)

    return

#TODO consider inf dist if segment in line of sight
#TODO imporove encroaching vertices removal performance
def split_edge(omesh, vertices, eh):
    def remove_encroaching_vertices(omesh, vertices, vh, h):
        def remove_one_encroaching_vertex(omesh, vertices, phw, h):
            for vh in omesh.vertices():
                if omesh.is_boundary(vh):
                    continue

                pcurr = np.array((omesh.point(vh)))

                #TODO consider inf dist if segment in line of sight
                if np.linalg.norm(phw - pcurr) < h:
                    remove_inner_vertex(omesh, vertices, vh)
                    return True

            return False
        phw = np.array((omesh.point(vh)))

        while remove_one_encroaching_vertex(omesh, vertices, phw, h):
            continue

        return
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

    # remove encroaching vertices
    if omesh.is_boundary(eh):
        h = np.linalg.norm(sv_new.XYZ_vec3() - sv0.XYZ_vec3())
        remove_encroaching_vertices(omesh, vertices, vh_new, h)

    return

#TODO make less stupid
def parse_back(omesh, face_mesh):
    vertices, _, triangles, _ = face_mesh

    vertices_tmp = []
    assert len(vertices) == len(omesh.vertices())
    for i in range(len(omesh.vertices())):
        sv = sv_from_vh(omesh, om.VertexHandle(i))
        vertices_tmp.append(sv)
    for i in range(len(vertices_tmp)):
        vertices[i] = vertices_tmp[i]

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
    flip_until_scdt(omesh)

    # step 2+3: find largest triangle that fails shape ans size criteria
    delta = find_largest_failing_triangle(omesh)

    global iter_counter
    iter_counter = 0
    while delta.is_valid() and iter_counter != MAX_ITERATIONS:
        print('iteration', iter_counter)

        handle, scc = calculate_refinement(omesh, delta)

        if isinstance(handle, om.FaceHandle):
            insert_inner_vertex(omesh, vertices, handle, scc)
        else:
            assert isinstance(handle, om.EdgeHandle)
            split_edge(omesh, vertices, handle)

        # after vertex insertion and possible deletion, restore SCDT criteria
        flip_until_scdt(omesh)

        # update for next iteration
        delta = find_largest_failing_triangle(omesh)
        iter_counter += 1

    parse_back(omesh, face_mesh)
    # vertices += DEBUG_VERTICES #TODO remove

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
        # MAX_ITERATIONS = TMP
        mesh = triangulate(INPUT_PATH)
        TMP2 = False
        write_to_file(mesh, OUTPUT_DIR)
