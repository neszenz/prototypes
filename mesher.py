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
import triangle as shewchuck_triangle

import paths
import sampler
from meshkD import SuperVertex, MeshkD

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_PATH = paths.PATH_98613
sampler.NUMBER_OF_SAMPLES = 20

OUTPUT_DIR = 'tmp'
# unitize timestamp prefix w/ sampler, so that output files have the same name
OUTPUT_PREFIX = sampler.OUTPUT_PREFIX

## BEGINE OF CHEW93_SURFACE  = + = + = + = + = + = + = + = + = + = + = + = + = +
# 1st: builds CDT meeting normal criteria; 2nd: filps edges until surface CDT
def triangulate(vertices, segments):
    triangles = []

    t = shewchuck_triangle.Triangle()

    points = [(sv.u, sv.v) for sv in vertices]
    markers = [0] * len(points)
    t.set_points(points, markers=markers)
    t.set_segments(segments)

    t.triangulate(mode='pzQ')

    for triNode in t.get_triangles():
        ([i0, i1, i2], neighbors, attri) = triNode
        triangles.append((i0, i1, i2))

    #TODO trim triangles and make scdt

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
    def segments_from_wires(wire_meshes):
        segments = []

        for wire_mesh in wire_meshes:
            for i in range(len(wire_mesh)):
                i0 = i
                i1 = (i+1) % len(wire_mesh)
                segments.append((i0, i1))

        return segments
    # step 1: initial surface Delaunay triangulation
    segments = segments_from_wires(wire_meshes)
    triangles = triangulate(vertices, segments)
    scdt = vertices, segments, triangles

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
        #scdt = vertices, segments, triangulate_scdt(vertices, segments)
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
    write_mesh_to_file(mesh)
