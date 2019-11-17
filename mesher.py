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
INPUT_PATH = paths.PATH_BOX
sampler.NUMBER_OF_SAMPLES = 10
dr2d.SIZE_THRESHOLD = 0.5

OUTPUT_DIR = os.path.join('tmp', 'mesher')
# unitize timestamp prefix w/ sampler, so that output files have the same name
OUTPUT_PREFIX = sampler.OUTPUT_PREFIX

def triangulate_v1(mesh1D):
    def collect_vertex_indices_of_face(face_mesh):
        vertices = []
        outer_wire, inner_wires, _ = face_mesh
        vertices += outer_wire
        for inner_wire in inner_wires:
            vertices += inner_wire
        return vertices
    def vertices_from_vertex_indices(mesh1D, vertex_indices):
        vertices = []
        for vertex_index in vertex_indices:
            vertices.append(mesh1D.vertices[vertex_index])
        return vertices
    def triangulate_face(vertices, vertex_offset):
        triangles = []
        t = shewchuck_triangle.Triangle()
        pointsBoundary = [[sv.u, sv.v] for sv in vertices]
        markersBoundary = [1] * len(pointsBoundary)
        pointsInner = []
        markersInner = [0] * len(pointsInner)
        points = pointsBoundary + pointsInner
        markers = markersBoundary + markersInner
        t.set_points(points, markers=markers)
        t.triangulate(mode='zQ')
        for triNode in t.get_triangles():
            ([i0, i1, i2], neighbors, attri) = triNode
            i0 += vertex_offset
            i1 += vertex_offset
            i2 += vertex_offset
            triangles.append((i0, i1, i2))
        return triangles
    vertices = []
    triangles = []
    for face_mesh in mesh1D.face_meshes:
        face_vertex_indices = collect_vertex_indices_of_face(face_mesh)
        face_vertices = vertices_from_vertex_indices(mesh1D, face_vertex_indices)
        triangles+= triangulate_face(face_vertices, len(vertices))
        vertices += face_vertices
    return vertices, triangles
def triangulate_v2(mesh1D):
    def pslg_from_face_mesh(mesh1D, face_index):
        def process_wire(vertices, wire_mesh):
            wire_vertices = []
            wire_segments = []
            for iIndex in range(0, len(wire_mesh)):
                i0 = iIndex
                i1 = (iIndex+1) % len(wire_mesh)
                wire_vertices.append(vertices[wire_mesh[i0]])
                wire_segments.append((i0, i1))
            return wire_vertices, wire_segments
        def apply_offset(segments, vertex_offset):
            offset_segments = []
            for segment in segments:
                offset_segment = (segment[0]+vertex_offset, segment[1]+vertex_offset)
                offset_segments.append((offset_segment))
            return offset_segments
        def make2D(super_vertices):
            vertices2D = []
            for sv in super_vertices:
                vertices2D.append(tuple(sv.UV_vec2()))
            return vertices2D
        face_vertices = []
        face_segments = []
        outer_wire, inner_wires, _ = mesh1D.face_meshes[face_index]
        outer_vertices, outer_segments = process_wire(mesh1D.vertices, outer_wire)
        face_vertices += outer_vertices
        face_segments += outer_segments
        for inner_wire in inner_wires:
            inner_vertices, inner_segments = process_wire(mesh1D.vertices, inner_wire)
            face_segments += apply_offset(inner_segments, len(face_vertices))
            face_vertices += inner_vertices
        face_vertices = make2D(face_vertices)
        return face_vertices, face_segments
    def apply_offset(triangles, vertex_offset):
        offset_triangles = []
        for triangle in triangles:
            i0, i1, i2 = triangle
            offset_triangle = (i0+vertex_offset, i1+vertex_offset, i2+vertex_offset)
            offset_triangles.append(offset_triangle)
        return offset_triangles
    def make3D(vertices, index):
        super_vertices = []
        for v in vertices:
            sv = SuperVertex(u=v[0], v=v[1], face_id=index+1)
            sampler.project_to_XYZ(sv, INPUT_PATH) #TODO super not final; temp!!1!
            super_vertices.append(sv)
        return super_vertices
    vertices = []
    face_meshes = []
    for face_index in range(0, mesh1D.number_of_faces()):
        _, _, face_type = mesh1D.face_meshes[face_index]
        print('face', face_index, 'of', mesh1D.number_of_faces(), '. . .')
        face_pslg = pslg_from_face_mesh(mesh1D, face_index)
        face_vertices, _, face_triangles = dr2d.chew93(face_pslg)
        face_vertices = make3D(face_vertices, face_index)
        face_triangles = apply_offset(face_triangles, len(vertices))
        face_meshes.append((face_triangles, face_type))
        vertices += face_vertices
    return vertices, face_meshes

def mesher(path, write_mesh1D=True):
    simple_sampler = sampler.factory(sampler.SAMPLER_TYPE.SIMPLE)
    mesh1D = simple_sampler(INPUT_PATH)
    if write_mesh1D:
        sampler.write_mesh_to_file(mesh1D)
    print('>> mesh samples w/ chew93_2D')
    vertices, face_meshes = triangulate_v2(mesh1D)
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
    mesh2D = mesher(INPUT_PATH, write_mesh1D=True)
    write_mesh_to_file(mesh2D)
