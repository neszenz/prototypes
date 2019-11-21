""" Chew93_Surface mesher for STEP files using sampler.py w/ OpenCASCADE
Contains functionality to generate a 2-dimensional mesh from a given path to a
STEP file.

To gain a foundation for the meshing process, sampler.py is used to generate a
1-dimensional mesh from the STEP file which is than iterated upon. Here, the
first step is to compute constrained Delaunay Triangulations w/ pytriangle and
then apply my implementation of Chew's Surface Delaunay Refinement algorithm.
"""
import numpy as np
import os
import pickle
import triangle as shewchuk_triangle
import openmesh

import paths
import sampler
from meshkD import SuperVertex, MeshkD

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_PATH = paths.PATH_BOX
sampler.NUMBER_OF_SAMPLES = 10

OUTPUT_DIR = 'tmp'
# unitize timestamp prefix w/ sampler, so that output files have the same name
OUTPUT_PREFIX = sampler.OUTPUT_PREFIX

def normalize(vector):
    vNorm = np.linalg.norm(vector)
    if vNorm == 0.0:
        return vector
    else:
        return vector / vNorm

def calculate_angle(p0, p1, p2):
    p10 = normalize(p0 - p1)
    p12 = normalize(p2 - p1)

    dp = np.dot(p10, p12)
    max(-1.0, min(1.0, dp))

    return np.arccos(dp)

## BEGINE OF CHEW93_SURFACE  = + = + = + = + = + = + = + = + = + = + = + = + = +
# 1st: builds CDT meeting normal criteria; 2nd: filps edges until surface CDT
def triangulate(vertices, pslg):
    def triangulate_cdt(vertices, pslg):
        segments, holes, boundary_offset = pslg
        triangles = []

        t = shewchuk_triangle.Triangle()

        points = [(sv.u, sv.v) for sv in vertices]
        markersBoundary = [1] * boundary_offset
        markersInner = [0] * (len(points)-boundary_offset)
        markers = markersBoundary + markersInner

        t.set_points(points, markers=markers)
        t.set_segments(segments)
        t.set_holes(holes)

        t.triangulate(mode='pzQ')

        for triNode in t.get_triangles():
            ([i0, i1, i2], neighbors, attri) = triNode
            triangles.append((i0, i1, i2))

        return triangles
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
        def flip_one_non_scdt_edge(omesh):
            def get_min_angle(omesh, eh):
                def get_min_face_angle(omesh, heh):
                    p0 = omesh.point(omesh.from_vertex_handle(heh))
                    heh = omesh.next_halfedge_handle(heh)
                    p1 = omesh.point(omesh.from_vertex_handle(heh))
                    heh = omesh.next_halfedge_handle(heh)
                    p2 = omesh.point(omesh.from_vertex_handle(heh))

                    alpha = calculate_angle(p0, p1, p2)
                    beta = calculate_angle(p1, p2, p0)
                    gamma = calculate_angle(p2, p0, p1)

                    assert np.allclose((alpha+beta+gamma), np.pi)

                    return min(alpha, beta, gamma)
                heh0 = omesh.halfedge_handle(eh, 0)
                heh1 = omesh.halfedge_handle(eh, 1)

                min_angle0 = get_min_face_angle(omesh, heh0)
                min_angle1 = get_min_face_angle(omesh, heh0)

                return min(min_angle0, min_angle1)
            for eh in omesh.edges():
                if omesh.is_boundary(eh):
                    continue
                # assert omesh.is_flip_ok(eh)
                if not omesh.is_flip_ok(eh): #TODO remove
                    continue

                min_angle_before = get_min_angle(omesh, eh)
                omesh.flip(eh)
                min_angle_after = get_min_angle(omesh, eh)

                print(min_angle_before, min_angle_after)
                if min_angle_after <= min_angle_before:
                    print('nope')
                    omesh.flip(eh) # redo if no improvment
                else:
                    print('flip')
                    return True

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

        # while flip_one_non_scdt_edge(omesh):
            # continue
        for i in range(100):
            flip_one_non_scdt_edge(omesh)

        scdt_triangles = triangles_from_omesh(omesh)

        return scdt_triangles
    triangles = triangulate_cdt(vertices, pslg)

    #TODO check/ enforce face-normal-criteria

    triangles = scdt_from_cdt(vertices, triangles)

    return triangles

def find_largest_failing_triangle(scdt):
    delta = -1

    #TODO

    return delta

# computes surface circumcenter in ambient space
def surface_circumcenter(scdt, delta):
    #TODO

    return np.array((0.0, 0.0, 0.0)) #TODO  stub

#TODO figure out how to do this
def travel(scdt, delta, c):
    #TODO

    return -1 #TODO stub

def split_segment(scdt, segment_index):
    #TODO

    return

def insert_circumcenter(scdt, c):
    svc = SuperVertex(x=c[0], y=c[1], z=c[2])

    #TODO

    return

def chew93_Surface(vertices, wire_meshes):
    def pslg_from_wires(wire_meshes):
        def calculate_hole(wire_mesh):
            #TODO

            return []
        segments = []
        holes = []
        boundary_offset = 0

        if len(wire_meshes) > 0:
            boundary_offset = len(wire_meshes[0])

        for wire_index in range(len(wire_meshes)):
            wire_mesh = wire_meshes[wire_index]
            wire_segments = []

            for i in range(len(wire_mesh)):
                i0 = wire_mesh[i]
                i1 = wire_mesh[(i+1) % len(wire_mesh)]

                if wire_index == 0:
                    assert max(i0, i1) < boundary_offset
                    wire_segments.append((i0, i1)) # keep ccw order
                else:
                    wire_segments.insert(0, (i1, i0)) # change cw to ccw for pytriangle
                    holes += calculate_hole(wire_mesh)

            segments += wire_segments

        return segments, holes, boundary_offset
    # step 1: initial surface Delaunay triangulation
    pslg = pslg_from_wires(wire_meshes)
    triangles = triangulate(vertices, pslg)
    scdt = vertices, pslg, triangles

    # step 2+3: find largest triangle that fails shape ans size criteria
    #delta = find_largest_failing_triangle(scdt)

    #while delta >= 0:
        # step 4: travel across the from any of delta's corners to c
        #c = surface_circumcenter(scdt, delta)
        #hit_segment_index = travel(scdt, delta, c)

        #if hit_segment_index >= 0:
            # step 6: split segment; remove encroaching inner vertices
            #split_segment(scdt, hit_segment_index)
        #else:
            # step 5: no segment was hit; insert c
            #insert_circumcenter(scdt, c)

        # update for next loop
        #scdt = vertices, pslg, triangulate_scdt(vertices, pslg)
        #delta = find_largest_failing_triangle(scdt)

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

    # print('>> mesh samples w/ chew93_2D')
    print('>> mesh samples as CDT per face')
    calculate_triangulation(mesh)

    return mesh

def write_mesh_to_file(mesh2D):
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
    pickle.dump(mesh2D, open(path, 'wb'))

    return

if __name__ == '__main__':
    mesh = mesher(INPUT_PATH)
    # write_mesh_to_file(mesh)
