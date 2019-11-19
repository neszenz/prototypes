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
import dr2d #TODO remove
from meshkD import SuperVertex, Mesh2D

## config and enum + = + = + = + = + = + = + = + = + = + = + = + = + = + = + = +
INPUT_PATH = paths.PATH_FILLET
sampler.NUMBER_OF_SAMPLES = 10
dr2d.SIZE_THRESHOLD = 10.0

OUTPUT_DIR = os.path.join('tmp', 'mesher')
# unitize timestamp prefix w/ sampler, so that output files have the same name
OUTPUT_PREFIX = sampler.OUTPUT_PREFIX

## BEGINE OF CHEW93_SURFACE  = + = + = + = + = + = + = + = + = + = + = + = + = +
# 1st: builds CDT meeting normal criteria; 2nd: filps edges until surface CDT
def triangulate(svertices, segments):
    triangles = []

    #TODO

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
    new_sv = SuperVertex(x=c[0], y=c[1], z=c[2])
    #TODO
    return

def chew93_Surface(pslg):
    svertices, segments = pslg

    # step 1: initial surface Delaunay triangulation
    scdt = svertices, segments, triangulate(svertices, segments)

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
        #scdt = svertices, segments, triangulate(svertices, segments)
        #delta = find_largest_failing_triangle(scdt)

    return scdt
## END OF CHEW93_SURFACE = + = + = + = + = + = + = + = + = + = + = + = + = + = +

def generate_mesh(mesh1D):
    def pslg_from_face_mesh(mesh1D, face_index):
        svertices = []
        segments = []
        pslg = (svertices, segments)

        #TODO

        return pslg
    def apply_offset(triangles, offset):
        #TODO
        return
    svertices = []
    face_meshes = []

    for face_index in range(mesh1D.number_of_faces()):
        print('face', face_index+1, 'of', mesh1D.number_of_faces(), '. . .')
        _, _, face_type = mesh1D.face_meshes[face_index]
        pslg = pslg_from_face_mesh(mesh1D, face_index)
        face_svertices, _, face_triangles = chew93_Surface(pslg)
        apply_offset(face_triangles, len(svertices))
        print() #TODO remove
        for sv in face_svertices: #TODO remove
            print(sv)
        for t in face_triangles: #TODO remove
            print(t)
        face_meshes.append((face_triangles, face_type))
        svertices += face_svertices

    return svertices, face_meshes

def mesher(path, write_mesh1D=True):
    simple_sampler = sampler.factory(sampler.SAMPLER_TYPE.SIMPLE)
    mesh1D = simple_sampler(INPUT_PATH)
    if write_mesh1D:
        sampler.write_mesh_to_file(mesh1D)
    print('>> mesh samples w/ chew93_2D')
    vertices, face_meshes = generate_mesh(mesh1D)
    return Mesh2D(path, vertices, face_meshes)

def write_mesh_to_file(mesh2D):
    def generate_output_file_path():
        name = OUTPUT_PREFIX + Mesh2D.FILE_EXTENSION
        path = os.path.join(OUTPUT_DIR, name)
        if os.path.exists(path):
            raise Exception('output file name already exists:', path)
        assert not os.path.exists(path)
        return path
    path = generate_output_file_path()
    print('>> write to file \"', path, '\"', sep='')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pickle.dump(mesh2D, open(path, 'wb'))

if __name__ == '__main__':
    mesh2D = mesher(INPUT_PATH, write_mesh1D=False)
    # write_mesh_to_file(mesh2D)
