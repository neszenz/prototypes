""" Chew93_Surface mesher for STEP files using sampler.py w/ OpenCASCADE
Contains functionality to generate a 2-dimensional mesh from a given path to a
STEP file.

To gain a foundation for the meshing process, sampler.py is used to generate a
1-dimensional mesh from the STEP file which is than iterated upon. Here, the
first step is to compute constrained Delaunay Triangulations w/ pytriangle and
then apply my implementation of Chew's Surface Delaunay Refinement algorithm.
"""
import datetime
from gerobust.predicates import clockwise, counter_clockwise
import numpy as np
import os
import pickle
import triangle as shewchuk_triangle
import openmesh

import paths
import sampler
from meshkD import SuperVertex, MeshkD
from util import *

from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Lin
from OCC.Core.IntCurvesFace import IntCurvesFace_Intersector
from OCC.Core.TopAbs import TopAbs_ON, TopAbs_IN

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_PATH = paths.STEP_SURFFIL
sampler.NUMBER_OF_SAMPLES = 10
sampler.INCLUDE_INNER_WIRES = True
SMALLEST_ANGLE = np.deg2rad(30)
SIZE_THRESHOLD = float('inf')

OUTPUT_DIR = paths.DIR_TMP

## BEGIN OF CHEW93_SURFACE + = + = + = + = + = + = + = + = + = + = + = + = + = +
# 1st: builds CDT meeting normal criteria; 2nd: filps edges until surface CDT
def triangulate(vertices, segment_loops):
    def triangulate_cdt(vertices, segment_loops):
        assert len(segment_loops) > 0
        triangles = []

        t = shewchuk_triangle.Triangle()

        boundary_offset = len(segment_loops[0])
        points = [(sv.u, sv.v) for sv in vertices]
        markersBoundary = [1] * boundary_offset
        markersInner = [0] * (len(points)-boundary_offset)
        markers = markersBoundary + markersInner
        segments = []
        for segment_loop in segment_loops:
            segments += segment_loop

        t.set_points(points, markers=markers)
        t.set_segments(segments)

        t.triangulate(mode='pzQ')

        for triNode in t.get_triangles():
            ([i0, i1, i2], neighbors, attri) = triNode
            triangles.append((i0, i1, i2))

        return triangles
    def trim_triangulation(triangles, segment_loops):
        def is_on_segment(i, segment_loop):
            for segment in segment_loop:
                if i in segment:
                    return True

            return False
        t_index = 0
        while t_index < len(triangles):
            i0, i1, i2 = triangles[t_index]

            for inner_segment_loop in segment_loops[1:]:
                if is_on_segment(i0, inner_segment_loop) and \
                   is_on_segment(i1, inner_segment_loop) and \
                   is_on_segment(i2, inner_segment_loop):
                    del triangles[t_index]
                    t_index -= 1
                    break

            t_index += 1

        return
    def scdt_from_cdt(vertices, cdt_triangles):
        def openmesh_from_cdt(vertices, triangles):
            omesh = openmesh.TriMesh()

            vertex_handles = []
            for sv in vertices:
                vertex_handles.append(omesh.add_vertex(sv.XYZ_vec3()))

            for triangle in cdt_triangles:
                i0, i1, i2 = triangle
                vh0 = vertex_handles[i0]
                vh1 = vertex_handles[i1]
                vh2 = vertex_handles[i2]
                omesh.add_face(vh0, vh1, vh2)

            return omesh
        def are_consistent(omesh, vertices):
            for vh in omesh.vertices():
                if not all(vertices[vh.idx()].XYZ_vec3() == omesh.point(vh)):
                    return False

            return True
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
                assert counter_clockwise(*t0) == True #TODO ran into errors maybe
                assert counter_clockwise(*t1) == True #TODO b/o degenerating flips

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
        def triangles_from_omesh(omesh):
            triangles = []

            for fh in omesh.faces():
                heh0 = omesh.halfedge_handle(fh)
                heh1 = omesh.next_halfedge_handle(heh0)
                heh2 = omesh.next_halfedge_handle(heh1)
                vh0 = omesh.from_vertex_handle(heh0)
                vh1 = omesh.from_vertex_handle(heh1)
                vh2 = omesh.from_vertex_handle(heh2)
                triangles.append((vh0.idx(), vh1.idx(), vh2.idx()))

            return triangles
        omesh = openmesh_from_cdt(vertices, cdt_triangles)
        assert are_consistent(omesh, vertices)

        # if TMP:
            # for i in range(TMP2):
                # flip_one_non_scdt_edge(omesh, vertices)
        # else:
        while flip_one_non_scdt_edge(omesh, vertices):
            continue

        scdt_triangles = triangles_from_omesh(omesh)

        return scdt_triangles

    triangles = triangulate_cdt(vertices, segment_loops)
    trim_triangulation(triangles, segment_loops)

    #TODO check/ enforce face-normal-criteria

    triangles = scdt_from_cdt(vertices, triangles)

    return triangles

def find_largest_failing_triangle(scdt):
    def shape_test(scdt, t_index):
        vertices, _, triangles = scdt
        i0, i1, i2 = triangles[t_index]
        p0 = vertices[i0].XYZ_vec3()
        p1 = vertices[i1].XYZ_vec3()
        p2 = vertices[i2].XYZ_vec3()

        alpha = calculate_angle_in_corner(p0, p1, p2)
        beta = calculate_angle_in_corner(p1, p2, p0)
        gamma = calculate_angle_in_corner(p2, p0, p1)

        return min(alpha, beta, gamma) > SMALLEST_ANGLE
    def size_test(scdt, t_index):
        vertices, _, triangles = scdt
        i0, i1, i2 = triangles[t_index]
        p0 = vertices[i0].XYZ_vec3()
        p1 = vertices[i1].XYZ_vec3()
        p2 = vertices[i2].XYZ_vec3()

        t_size = calculate_area(p0, p1, p2)

        return t_size <= SIZE_THRESHOLD, t_size
    _, _, triangles = scdt

    delta_index = -1
    delta_size = 0.0

    for t_index in range(len(triangles)):
        well_shaped = shape_test(scdt, t_index)
        well_sized, t_size = size_test(scdt, t_index)

        if well_shaped and well_sized:
            continue
        else:
            if t_size > delta_size:
                delta_size = t_size
                delta_index = t_index
            else:
                assert delta_index >= 0

    return delta_index

def calculate_circumcenter(scdt, delta_index):
    vertices, _, triangles = scdt
    i0, i1, i2 = triangles[delta_index]
    A = vertices[i0].XYZ_vec3()
    B = vertices[i1].XYZ_vec3()
    C = vertices[i2].XYZ_vec3()

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

# computes surface circumcenter in ambient space
def calculate_surface_circumcenter(scdt, delta_index):
    def calculate_center_line(scdt, delta_index):
        vertices, _, triangles = scdt
        i0, i1, i2 = triangles[delta_index]
        a = vertices[i0].XYZ_vec3()
        b = vertices[i1].XYZ_vec3()
        c = vertices[i2].XYZ_vec3()

        cc = calculate_circumcenter(scdt, delta_index)
        n = normalize(np.cross(b-a, c-a))

        return cc, gp_Lin(gp_Pnt(cc[0], cc[1], cc[2]), gp_Dir(n[0], n[1], n[2]))
    def array_from_pnt(pnt):
        return np.array((pnt.X(), pnt.Y(), pnt.Z()))
    vertices, _, triangles = scdt
    face = vertices[triangles[delta_index][0]].face

    cc, center_line = calculate_center_line(scdt, delta_index)
    intersector = IntCurvesFace_Intersector(face, 0.0001)
    intersector.Perform(center_line, -float('inf'), float('inf'))
    if intersector.IsDone():
        if intersector.NbPnt() < 1:
            return None
        elif intersector.NbPnt() == 1:
            pnt = intersector.Pnt(1)
            assert intersector.State(1) == TopAbs_IN

            scc = np.array((pnt.X(), pnt.Y(), pnt.Z()))

            return scc
        else:
            points = [array_from_pnt(intersector.Pnt(i+1)) for i in range(intersector.NbPnt())]

            closest_index = 0
            closest_dist = np.linalg.norm(points[closest_index] - cc)

            for p_index in range(1, len(points)):
                p = points[p_index]
                current_dist = np.linalg.norm(p - cc)

                if current_dist < closest_dist:
                    closest_index = p_index
                    closest_dist = current_dist

            return points[closest_index]
    else:
        raise Exception('calculate_surface_circumcenter() error - intersector not done')

def longest_edge_vertex_indices(scdt, delta_index): #TODO check whether this works correctly
    vertices, _, triangles = scdt
    i0, i1, i2 = triangles[delta_index]
    p0 = vertices[i0].XYZ_vec3()
    p1 = vertices[i1].XYZ_vec3()
    p2 = vertices[i2].XYZ_vec3()

    p01_length = np.linalg.norm(p1-p0)
    p12_length = np.linalg.norm(p2-p1)
    p20_length = np.linalg.norm(p0-p2)

    if p01_length > p12_length:
        if p01_length > p20_length:
            return i0, i1
        else:
            assert p20_length > p01_length

            return i2, i0
    else:
        if p12_length > p20_length:
            return i1, i2
        else:
            assert p20_length > p12_length

            return i2, i0

    return (-1, -1)

#TODO check for guaranteed termination
#TODO avoid problem regions? (p.40)
def segment_index_of_longest_edge(scdt, edge_vertex_indices):
    _, segment_loops, _ = scdt
    ev0, ev1 = edge_vertex_indices

    for l_index in range(len(segment_loops)):
        segment_loop = segment_loops[l_index]
        for s_index in range(len(segment_loop)):
            s0, s1 = segment_loop[s_index]
            if (ev0 == s0 and ev1 == s1) or\
               (ev0 == s1 and ev1 == s0):
                return l_index, s_index

    return None

def insert_halfway_vertex_of_edge(scdt, edge_vertex_indices):
    vertices, _, _ = scdt
    assert len(vertices) > 0

    ev0, ev1 = edge_vertex_indices
    sv0 = vertices[ev0]
    sv1 = vertices[ev1]

    p0 = sv0.UV_vec2()
    p1 = sv1.UV_vec2()
    p01 = p1 - p0
    hw = p0 + (p01/2)

    svhw = SuperVertex(u=hw[0], v=hw[1])
    svhw.face = sv0.face
    svhw.face_id = sv0.face_id
    svhw.project_to_XYZ()

    vertices.append(svhw)

    return

def split_segment(scdt, segment_index):
    def count_segments(segment_loops):
        nos = 0

        for segment_loop in segment_loops:
            nos += len(segment_loop)

        return nos
    vertices, segment_loops, _ = scdt
    l_index, s_index = segment_index
    print('split_segment()')
    print(segment_index, '->', segment_loops[l_index][s_index])

    segment_loop = segment_loops[l_index]
    vi0, vi1 = segment_loop[s_index]

    # compute and insert halfway vertex
    sv0 = vertices[vi0]
    sv1 = vertices[vi1]
    sv_halfway = SuperVertex.compute_halfway_on_shared_edge(sv0, sv1)
    # segment vertices are never deleted -> insert before inner vertices to avoid index shift
    hw_index = count_segments(segment_loops) # == num of boundary vertices
    vertices.insert(hw_index, sv_halfway)

    # delete old and insert new segments
    print('del', segment_loop[s_index])
    del segment_loop[s_index]
    new_s0 = (vi0, hw_index)
    new_s1 = (hw_index, vi1)
    print('insert', new_s1)
    print('insert', new_s0)
    segment_loop.insert(s_index, new_s1)
    segment_loop.insert(s_index, new_s0)

    # delete encroaching vertices
    new_segments_length = np.linalg.norm(sv0.XYZ_vec3() - sv_halfway.XYZ_vec3())
    inner_vertices_offset = hw_index+1
    sv_index = inner_vertices_offset
    while sv_index < len(vertices):
        sv = vertices[sv_index]
        dist_to_hw = np.linalg.norm(sv.XYZ_vec3() - sv_halfway.XYZ_vec3())

        if dist_to_hw < new_segments_length:
            del vertices[sv_index]
        else:
            sv_index += 1

    return

def insert_inner_vertex(scdt, c):
    vertices, _, _ = scdt
    assert len(vertices) > 0

    svc = SuperVertex(x=c[0], y=c[1], z=c[2])
    svc.face_id = vertices[0].face_id
    svc.face = vertices[0].face
    svc.project_to_UV()

    vertices.append(svc)

    return

def chew93_Surface(vertices, wire_meshes):
    def segments_from_wires(vertices, wire_meshes):
        segment_loops = []

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

            segment_loops.append(wire_segments)

        return segment_loops
    # step 1: initial surface Delaunay triangulation
    segment_loops = segments_from_wires(vertices, wire_meshes)
    triangles = triangulate(vertices, segment_loops)
    scdt = vertices, segment_loops, triangles

    # step 2+3: find largest triangle that fails shape ans size criteria
    delta_index = find_largest_failing_triangle(scdt)
    #TODO remove
    #cc = calculate_circumcenter(scdt, delta_index)
    #p0 = vertices[triangles[delta_index][0]].XYZ_vec3()
    #p1 = vertices[triangles[delta_index][1]].XYZ_vec3()
    #p2 = vertices[triangles[delta_index][2]].XYZ_vec3()
    #n = normalize(calculate_normal(p0, p1, p2))
    #insert_inner_vertex(scdt, cc+n)
    #TODO remove

    iteration_counter = 0
    while delta_index >= 0:
        print('iteration', iteration_counter)
        iteration_counter += 1
        # step 4: travel across the from any of delta's corners to c
        try:
            c = calculate_surface_circumcenter(scdt, delta_index)
        except Exception as e:
            print(e)
            return triangles

        if c is None:
            edge_vertex_indices = longest_edge_vertex_indices(scdt, delta_index)
            s_index = segment_index_of_longest_edge(scdt, edge_vertex_indices)

            if s_index is None:
                # custom step for 'no intersection, no segment' case
                insert_halfway_vertex_of_edge(scdt, edge_vertex_indices)
            else:
                # step 6: split segment; remove encroaching inner vertices
                split_segment(scdt, s_index)
        else:
            # step 5: no segment was hit; insert c
            insert_inner_vertex(scdt, c)

        # update for next loop
        triangles = triangulate(vertices, segment_loops)
        scdt = vertices, segment_loops, triangles
        delta_index = find_largest_failing_triangle(scdt)

    return triangles
## END OF CHEW93_SURFACE = + = + = + = + = + = + = + = + = + = + = + = + = + = +

def calculate_triangulation(mesh):
    for face_index in range(mesh.number_of_faces()):
        print('face', face_index+1, 'of', mesh.number_of_faces(), '. . .')
        vertices, wire_meshes, triangles, _ = mesh.face_meshes[face_index]

        triangles.clear()
        triangles += chew93_Surface(vertices, wire_meshes)

    return

def mesher(path, write_mesh1D=True):
    simple_sampler = sampler.factory(sampler.SAMPLER_TYPE.SIMPLE)
    mesh = simple_sampler(INPUT_PATH)

    print('>> mesh samples with chew93_Surface')
    calculate_triangulation(mesh)
    mesh.reset_bounding_boxes()

    return mesh

def write_mesh_to_file(mesh2D):
    def generate_output_file_path():
        timestamp = datetime.datetime.now()
        name = timestamp.strftime('%y%m%d_%H%M%S') + MeshkD.FILE_EXTENSION
        path = os.path.join(OUTPUT_DIR, name)

        if os.path.exists(path):
            name = timestamp.strftime('%y%m%d_%H%M%S%f') + MeshkD.FILE_EXTENSION
            path = os.path.join(OUTPUT_DIR, name)

        assert not os.path.exists(path)

        return path
    path = generate_output_file_path()

    print('>> write to file \"', path, '\"', sep='')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pickle.dump(mesh2D, open(path, 'wb'))

    return

if __name__ == '__main__':
    TMP = False
    for i in range(1):
        TMP2 = i
        print(TMP2)
        mesh = mesher(INPUT_PATH)
        write_mesh_to_file(mesh)
        TMP = False
